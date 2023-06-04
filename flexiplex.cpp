// Copyright 2022 Nadia Davidson 
// This program is distributed under the MIT License.
// We also ask that you cite this software in publications
// where you made use of it for any part of the data analysis.

#include <iostream>
#include <fstream>
#include <istream>
#include <string>
#include <sstream>
#include <unordered_set>
#include <unordered_map>
#include <stdlib.h>
#include <unistd.h>
#include <vector>
#include <algorithm>
#include <thread>
#include <numeric>
#include <Rcpp.h>
#include <ctime>

#include "edlib.h"

using namespace std;

const static string VERSION="0.96.2";


// compliment nucleotides - used to reverse compliment string
char compliment(char& c){
  switch(c){
  case 'A' : return 'T';
  case 'T' : return 'A';
  case 'G' : return 'C';
  case 'C' : return 'G';
  default: return 'N';
  }
}

// Reverse complement
std::string reverse_compliment(const std::string& seq) {
  std::string rev_seq(seq.rbegin(), seq.rend());
  std::transform(rev_seq.begin(), rev_seq.end(), rev_seq.begin(),
                 [](char c) { return compliment(c); });
  return rev_seq;
}

//Holds the search string patterns
struct SearchSeq {
  string primer;
  string polyA;
  string umi_seq;
  string temp_barcode;
} search_pattern ;

//Holds the found barcode and associated information 
struct Barcode {
  string barcode;
  string umi;
  int editd;
  int flank_editd;
  int flank_start;
  int flank_end;
  bool unambiguous;
} ;

struct SearchResult {
  string read_id;
  string qual_scores;
  string line;
  string rev_line;
  vector<Barcode> vec_bc_for;
  vector<Barcode> vec_bc_rev;
};

// Code for fast edit distance calculation for short sequences modified from
// https://en.wikibooks.org/wiki/Algorithm_Implementation/Strings/Levenshtein_distance#C++
// s2 is always assumned to be the shorter string (barcode)
unsigned int edit_distance(const std::string& s1, const std::string& s2, unsigned int &end, int max_editd){

  const std::string_view s1_view(s1);
  const std::string_view s2_view(s2);

  const std::size_t len1 = s1_view.size() + 1;
  const std::size_t len2 = s2_view.size() + 1;

  std::vector<unsigned int> dist_holder(len1 * len2);
  //initialise the edit distance matrix.
  //penalise for gaps at the start and end of the shorter sequence (j)
  //but not for shifting the start/end of the longer sequence (i,0)
  dist_holder[0] = 0; //[0][0]
  for (std::size_t j = 1; j < len2; ++j)
    dist_holder[j] = j; //[0][j];
  for (std::size_t i = 1; i < len1; ++i)
    dist_holder[i * len2] = 0; //[i][0];

  unsigned int best = len2;
  end = len1 - 1;

  // loop over the distance matrix elements and calculate running distance
  for (std::size_t j = 1; j < len2; ++j) {
    bool any_below_threshold = false; // flag used for early exit
    for (std::size_t i = 1; i < len1; ++i) {
      unsigned int sub = (s1_view[i - 1] == s2_view[j - 1]) ? 0 : 1; // match / mismatch score

      const unsigned int &top_left = dist_holder[(i - 1) * len2 + (j - 1)];
      const unsigned int &left = dist_holder[i * len2 + (j - 1)];
      const unsigned int &top = dist_holder[(i - 1) * len2 + j];

      unsigned int min_value = std::min({top + 1, left + 1, top_left + sub});
      dist_holder[i * len2 + j] = min_value;

      if (min_value <= max_editd)
        any_below_threshold = true;
      if (j == (len2 - 1) && min_value < best) {
        // if this is the last row in j
        // check if this is the best running score
        // update the end position of alignment
        best = min_value;
        end = i;
      }
    }
    if(!any_below_threshold){ //early exit to save time.
      return(100);
    }
  }
  return best; //return edit distance
}


//Given a string 'seq' search for substring with primer and polyT sequence followed by
//a targeted search in the region for barcode
//Seaquence seearch is performed using edlib

Barcode get_barcode(const string & seq,
		    const unordered_set<string> &known_barcodes,
		    int flank_max_editd,
		    int barcode_max_editd){
  
  const int OFFSET=5; //wiggle room in bases of the expected barcode start site to search.

  //initialise struct variables for return:
  Barcode barcode;
  barcode.editd=100; barcode.flank_editd=100; barcode.unambiguous = false;

  //initialise edlib configuration
  EdlibEqualityPair additionalEqualities[5] = {{'?','A'},{'?','C'},{'?','G'},{'?','T'},{'?','N'}};
  EdlibAlignConfig edlibConf = {flank_max_editd, EDLIB_MODE_HW, EDLIB_TASK_PATH, additionalEqualities, 5};

  //search for primer and ployT (barcode and umi as wildcards)
  string search_string=
    search_pattern.primer+
    search_pattern.temp_barcode+
    search_pattern.umi_seq+
    search_pattern.polyA;
  EdlibAlignResult result = edlibAlign(search_string.c_str(), search_string.length(), seq.c_str(), seq.length(), edlibConf);
  if(result.status != EDLIB_STATUS_OK || result.numLocations==0 ){
    edlibFreeAlignResult(result);
    return(barcode); // no match found - return
  } //fill in info about found primer and polyT location
  barcode.flank_editd=result.editDistance;
  barcode.flank_start=result.startLocations[0];
  barcode.flank_end=result.endLocations[0];

  // Extract sub-patterns from aligment directly
  vector<long unsigned int> subpattern_lengths = {
    search_pattern.primer.length(),
    search_pattern.temp_barcode.length(),
    search_pattern.umi_seq.length(),
    search_pattern.polyA.length()
  };

  std::vector<long unsigned int> subpattern_ends;
  subpattern_ends.resize(subpattern_lengths.size());
  std::inclusive_scan(subpattern_lengths.begin(), subpattern_lengths.end(), subpattern_ends.begin());

  vector<int> read_to_subpatterns;
  read_to_subpatterns.reserve(subpattern_ends.size());

  // initialise pointers
  int i_read = barcode.flank_start;
  int i_pattern = 0;
  int i_subpattern = 0;

  // walk through edlib aligment
  // 0 for match
  // 1 for insertion to target
  // 2 for insertion to query
  // 3 for mismatch 
  std::vector<unsigned char> alignment_vector(result.alignment, result.alignment + result.alignmentLength);
  for (const auto& value : alignment_vector) {
    if (value != 1) {
      i_read++;
    }
    if (value != 2) {
      i_pattern++;
    }
    if (i_pattern >= subpattern_ends[i_subpattern]) {
      read_to_subpatterns.emplace_back(i_read);
      i_subpattern++;
    }
  }

  edlibFreeAlignResult(result);
  
  //if not checking against known list of barcodes, return sequence after the primer
  //also check for a perfect match straight up as this will save computer later.
  string exact_bc=seq.substr(read_to_subpatterns[0], read_to_subpatterns[1] - read_to_subpatterns[0]);
  if(known_barcodes.empty() || (known_barcodes.find(exact_bc) != known_barcodes.end())){ 
    barcode.barcode=exact_bc;
    barcode.editd=0;
    barcode.unambiguous=true;
    barcode.umi=seq.substr(read_to_subpatterns[1], search_pattern.umi_seq.length());
    return(barcode);
  }
  
  // otherwise widen our search space and the look for matches with errors
  string barcode_seq=seq.substr(read_to_subpatterns[0]-OFFSET,search_pattern.temp_barcode.length()+2*OFFSET);
  
  //iterate over all the known barcodes, checking each sequentially
  unsigned int editDistance; unsigned int endDistance;
  for (const auto &known_bc : known_barcodes) {
    editDistance = edit_distance(barcode_seq,known_bc,endDistance,barcode_max_editd);
    if (editDistance == barcode.editd) {
      barcode.unambiguous=false;    
    } else if (editDistance < barcode.editd && editDistance <= barcode_max_editd) { //if best so far, update
      barcode.unambiguous=true;
      barcode.editd=editDistance;
      barcode.barcode=known_bc;
      barcode.umi=seq.substr(read_to_subpatterns[0]-OFFSET+endDistance,search_pattern.umi_seq.length());//assumes no error in UMI seq.
      if(editDistance==0) { //if perfect match is found we're done.
        return(barcode);
      }
    }
  }
  return(barcode); //return the best matched barcode and associated information
}

//search a read for one or more barcodes (parent function that calls get_barcode)
vector<Barcode> big_barcode_search(const string & sequence,const unordered_set<string> & known_barcodes,
				   int max_flank_editd, int max_editd)
  vector<Barcode> return_vec; //vector of all the barcodes found

  //search for barcode
  Barcode result=get_barcode(sequence,known_barcodes,max_flank_editd,max_editd); //,ss);
  if(result.editd<=max_editd && result.unambiguous) //add to return vector if edit distance small enough
    return_vec.emplace_back(result);
  
  //if a result was found, mask out the flanking sequence and search again in case there are more.
  if (!return_vec.empty()) {
    string masked_sequence = sequence;
    for(const auto &barcode : return_vec){
      int flank_length = barcode.flank_end - barcode.flank_start;
      masked_sequence.replace(barcode.flank_start, flank_length,string(flank_length,'X'));
    } //recursively call this function until no more barcodes are found
    vector<Barcode> masked_res = big_barcode_search(masked_sequence,known_barcodes,max_flank_editd,max_editd);
    return_vec.insert(return_vec.end(),masked_res.begin(),masked_res.end()); //add to result
  }
    
  return(return_vec);
}

// utility function to check true/false input options
bool get_bool_opt_arg(string value){
  transform(value.begin(), value.end(), value.begin(), ::tolower);
  if( value.compare("true")==0 | value.compare("t")==0 | value.compare("1")==0){
    return true;
  } else if (value.compare("false")!=0 | value.compare("f")!=0 | value.compare("0")!=0){
    return false;
  } else {
    Rcpp::Rcout << "Unknown argument to boolean option\n";
    print_usage();
    exit(1);
  } 
}

// print information about barcodes
void print_stats(const string& read_id, const vector<Barcode>& vec_bc, ostream& out_stream) {
  for (const auto& barcode : vec_bc) {
    out_stream << read_id << '\t'
               << barcode.barcode << "\t"
	             << barcode.flank_editd << "\t"
	             << barcode.editd << "\t"
	             << barcode.umi << '\n';
  }
}

void print_line(const string& id, const string& read, const string& quals, ostream& out_stream) {
  const char delimiter = quals.empty() ? '>' : '@';
  
  //output to the new read file
  out_stream << delimiter << id << '\n'
             << read << '\n';
  
  if (!quals.empty())
    out_stream << '+' << id << '\n'
               << quals << '\n';
}

//print fastq or fasta lines..
void print_read(string read_id,string read, string qual,
		vector<Barcode> & vec_bc, string prefix,
		bool split, unordered_set<string> & found_barcodes,
		bool trim_barcodes){
  //loop over the barcodes found... usually will just be one
  for(int b=0; b<vec_bc.size() ; b++){
    
    //format the new read id. Using FLAMES format.
    stringstream ss;
    ss << (b+1) << "of" << vec_bc.size() ;
    string barcode=vec_bc.at(b).barcode;
    string new_read_id=barcode+"_"+vec_bc.at(b).umi+"#"+read_id+ss.str();
    
    // work out the start and end base in case multiple barcodes
    int read_start=vec_bc.at(b).flank_end;
    int read_length=read.length()-read_start;
    for(int f=0; f<vec_bc.size(); f++){
      int temp_read_length=vec_bc.at(f).flank_start-read_start;
      if(temp_read_length>0 && temp_read_length<read_length)
	read_length=temp_read_length;
    }
    string qual_new=""; //don't trim the quality scores if it's a fasta file
    if(qual!="") qual_new=qual.substr(read_start,read_length);
    string read_new=read.substr(read_start,read_length);

    if(b==0 && !trim_barcodes){ //override if read shouldn't be cut
      new_read_id=read_id;
      read_new=read;
      qual_new=qual;
      b=vec_bc.size(); //force loop to exit after this iteration
    }
    
    if(split){ //to a file if spliting by barcode
      string outname=prefix+"_"+barcode+".";
      if(qual=="") outname+="fasta"; else outname+="fastq";
      ofstream outstream;
      if(found_barcodes.insert(barcode).second)
	outstream.open(outname); //override file if this is the first read for the barcode
      else
	outstream.open(outname,ofstream::app); //remove file if this is the first read for the barcode
      print_line(new_read_id,read_new,qual_new,outstream);
      outstream.close();
    } else {
      print_line(new_read_id,read_new,qual_new,std::cout);
    }
  }
}

// separated out from main so that this can be run with threads
void search_read(vector<SearchResult> & reads, unordered_set<string> & known_barcodes,
			 int flank_edit_distance, int edit_distance){
  for (auto& read : reads) {
    //forward search
    read.vec_bc_for = big_barcode_search(read.line, known_barcodes, flank_edit_distance, edit_distance);

    //Check the reverse compliment of the read
    read.rev_line = reverse_compliment(read.line);
    read.vec_bc_rev = big_barcode_search(read.rev_line, known_barcodes, flank_edit_distance, edit_distance);
  }
}

// [[Rcpp::export]]
int flexiplex(){
  std::ios_base::sync_with_stdio(false);

  Rcpp::Rcout << "FLEXIPLEX " << VERSION << "\n";

  //Variables to store user options
  //Set these to their defaults
  int expected_cells=0; //(d)
  int edit_distance=2; //(e)
  int flank_edit_distance=8; //(f)

  //set the output filenames
  string out_stat_filename="reads_barcodes.txt";
  string out_bc_filename="barcodes_counts.txt";
  string out_filename_prefix="flexiplex"; //(n)

  bool split_file_by_barcode=false; //(s)
  bool remove_barcodes=true; //(r)
  
  search_pattern.primer = "CTACACGACGCTCTTCCGATCT"; //(p)
  search_pattern.polyA = string(9,'T'); //(T)
  search_pattern.umi_seq = string(12,'?'); //(length u)
  search_pattern.temp_barcode = string(16,'?'); //(length b)
  
  //Set of known barcodes 
  unordered_set<string> known_barcodes;
  unordered_set<string> found_barcodes;

  // threads
  int n_threads=1;
  
  /*** Pass command line option *****/
  int c;
  int params=1;
  ifstream file;
  string line;

  while((c =  getopt(argc, argv, "k:i:l:r:b:u:e:f:n:s:h:p:")) != EOF){
    switch(c){
    case 'k': { //k=list of known barcodes
      string file_name(optarg);
      string bc;
      /**** READ BARCODE LIST FROM FILE ******/
      file.open(file_name);
      Rcpp::Rcout << "Setting known barcodes from "<< file_name << "\n";
      if(!(file.good())){ //if the string given isn't a file
	stringstream bc_list(file_name); string s;
	while (getline(bc_list, bc, ',')) //tokenize
	  known_barcodes.insert(bc);
      } else {
	// otherwise get the barcodes from the file..
	while ( getline (file,line) ){
	  istringstream line_stream(line);
	  line_stream >> bc;
	  known_barcodes.insert(bc); 
	}
	file.close();
      }
      Rcpp::Rcout << "Number of known barcodes: " << known_barcodes.size() << "\n";
      if(known_barcodes.size()==0){
	print_usage();
	exit(1); //case barcode file is empty
      }
      //set barcode length automatically from known barcodes..
      int bl=(known_barcodes.begin())->length();
      search_pattern.temp_barcode=string(bl,'?');
      Rcpp::Rcout << "Setting barcode length automatically to " << bl << "\n";
      params+=2;
      break;     
    }
    case 'i':{
      remove_barcodes=get_bool_opt_arg(optarg);
      Rcpp::Rcout << "Setting read IDs to be replaced: "<< remove_barcodes << "\n";
      params+=2;
      break;
    }
    case 'e':{
      edit_distance=atoi(optarg);
      Rcpp::Rcout << "Setting max barcode edit distance to "<< edit_distance << "\n";
      params+=2;
      break;
    }
    case 'f':{
      flank_edit_distance=atoi(optarg);
      Rcpp::Rcout << "Setting max flanking sequence edit distance to "<< flank_edit_distance << "\n";
      params+=2;
      break;
    }
    case 'l':{
      search_pattern.primer=optarg;
      Rcpp::Rcout << "Setting primer to search for: " << search_pattern.primer << "\n";
      params+=2;
      break;
    }
    case 'r':{
      search_pattern.polyA=optarg;
      Rcpp::Rcout << "Setting polyT to search for: " << search_pattern.polyA << "\n";
      params+=2;
      break;
    }
    case 'u':{
      int ul=atoi(optarg);
      search_pattern.umi_seq=string(ul,'?');
      Rcpp::Rcout << "Setting UMI length to " << ul << "\n";
      params+=2;
      break;
    }
    case 'b':{
      int bl=atoi(optarg);
      search_pattern.temp_barcode=string(bl,'?');
      Rcpp::Rcout << "Setting barcode length to " << bl << "\n";
      params+=2;
      break;
    }
    case 'h':{
      print_usage();
      exit(1);
    }
    case 'n':{
      out_filename_prefix=optarg;
      Rcpp::Rcout << "Setting output filename prefix to: " << out_filename_prefix << "\n";
      params+=2;
      break;
    }
    case 's':{
      split_file_by_barcode=get_bool_opt_arg(optarg);
      Rcpp::Rcout << "Split read output into separate files by barcode: " << split_file_by_barcode << "\n";
      int max_split_bc=50;
      if(known_barcodes.size()>max_split_bc){
	Rcpp::Rcout << "Too many barcodes to split into separate files: "<< known_barcodes.size()
	     << "> "<< max_split_bc<< "\n";
	split_file_by_barcode=false;
      }
      params+=2;
      break;
    }
    case 'p':{
      n_threads=atoi(optarg);
      Rcpp::Rcout << "Setting number of threads to "<< n_threads << "\n";
      params+=2;
      break;
    }
    case '?': //other unknown options
      Rcpp::Rcout << "Unknown option.. stopping" << "\n";
      print_usage();
      exit(1);
    }
  }
  
  Rcpp::Rcout << "For usage information type: flexiplex -h" << "\n";
  
  istream * in;
  FILE * ifile;
    
  //check that a read file is given
  if(params>=argc){
    Rcpp::Rcout << "No filename given... getting reads from stdin..." << "\n";
    in=&std::cin;
  } else {
    // check that the reads fileis okay
    string reads_file=argv[params];
    file.open(reads_file);
    if(!(file.good())){
      Rcpp::Rcout << "Unable to open file " << reads_file << "\n";
      print_usage();
      exit(1);
    }
    in=&file;
  }
  
  /********* FIND BARCODE IN READS ********/
  string sequence;
  int bc_count=0;
  int r_count=0;
  int multi_bc_count=0;
  
  ofstream out_stat_file;
  out_stat_filename=out_filename_prefix+"_"+out_stat_filename;
  out_bc_filename=out_filename_prefix+"_"+out_bc_filename;
  params+=2;

  if(known_barcodes.size()>0){
    out_stat_file.open(out_stat_filename);
    out_stat_file << "Read\tCellBarcode\tFlankEditDist\tBarcodeEditDist\tUMI"<<"\n";
  }
  Rcpp::Rcout << "Searching for barcodes..." << "\n";
  bool is_fastq=true;
  unordered_map< string, int > barcode_counts; 
  string read_id_line;
  if(getline (*in,read_id_line)){ //check the first line for file type
    if(read_id_line[0]=='>'){ is_fastq=false;
    } else if (read_id_line[0] == '@'){ //fasta
    } else {
      Rcpp::Rcout << "Unknown read format... exiting" << "\n"; exit(1);
    }
  }
  
  while ( getline (*in,line) ){
    const int buffer_size = 2000; //number of reads to pass to one thread.
    vector<vector<SearchResult>> sr_v(n_threads);
    for(int i=0; i<n_threads; i++)
      sr_v[i]=vector<SearchResult>(buffer_size); 
    vector<thread> threads(n_threads);
    for(int t=0; t < n_threads; t++){ //get n_threads*buffer number or reads..
      for(int b=0 ; b < buffer_size ; b++){ 
	SearchResult & sr = sr_v[t][b];
	sr.line=line;
	string read_id;
	//sr.read_id= read_id_line.substr(1,read_id_line.find_first_not_of(" \t")-1);      
	istringstream line_stream(read_id_line);
	line_stream >> sr.read_id;
	sr.read_id.erase(0,1);

	    
	//      string qual_scores="";
	if(!is_fastq){ //fastq (account for multi-lines per read)
	  string buffer_string; 
	  while(getline(*in,buffer_string) && buffer_string[0]!='>')
	    sr.line+=buffer_string;
	  read_id_line=buffer_string;
	} else { //fastq (get quality scores)
	  for(int s=0; s<2; s++) getline(*in,sr.qual_scores);
	  getline(*in,read_id_line);
	}
      
	r_count++; //progress counter
	if(r_count % 100000 == 0)
	  Rcpp::Rcout << r_count/((double) 1000000 ) << " million reads processed.." << "\n";

	//this is quite ugly, must be a better way to do this..
	if(b==buffer_size-1 && t==n_threads-1){
	  break; //if it's the last in the chunk don't getline as this happens in the while statement
	} else if( !getline(*in,line)){ //case we are at the end of the reads.
	  sr_v[t].resize(b+1);
	  threads[t]=std::thread(search_read,ref(sr_v[t]),ref(known_barcodes),flank_edit_distance,edit_distance);
	  for(int t2=t+1; t2 < n_threads ; t2++) sr_v[t2].resize(0);
	  goto print_result; //advance the line
	}
      }
      // send reads to the thread
      threads[t]=std::thread(search_read,ref(sr_v[t]),ref(known_barcodes),flank_edit_distance,edit_distance);
    }
  print_result:
    
    for(int t=0; t < sr_v.size(); t++){ //loop over the threads and print out ther results
      if(sr_v[t].size()>0) threads[t].join(); // wait for the threads to finish before printing
      
      for(int r=0; r< sr_v[t].size(); r++){ // loop over the reads
	
	for(int b=0; b<sr_v[t][r].vec_bc_for.size(); b++)
	  barcode_counts[sr_v[t][r].vec_bc_for.at(b).barcode]++;
	for(int b=0; b<sr_v[t][r].vec_bc_rev.size(); b++)
	  barcode_counts[sr_v[t][r].vec_bc_rev.at(b).barcode]++;
	
	if((sr_v[t][r].vec_bc_for.size()+sr_v[t][r].vec_bc_rev.size())>0)
	  bc_count++;
	if((sr_v[t][r].vec_bc_for.size()+sr_v[t][r].vec_bc_rev.size())>1 ){
	  multi_bc_count++;
	}
	
	if(known_barcodes.size()!=0){ // if we are just looking for all possible barcodes don't output reads etc.
	  
	  print_stats(sr_v[t][r].read_id, sr_v[t][r].vec_bc_for, out_stat_file);
	  print_stats(sr_v[t][r].read_id, sr_v[t][r].vec_bc_rev, out_stat_file);
	  
	  print_read(sr_v[t][r].read_id+"_+",sr_v[t][r].line,sr_v[t][r].qual_scores,sr_v[t][r].vec_bc_for,
		     out_filename_prefix,split_file_by_barcode,found_barcodes,remove_barcodes);
	  reverse(sr_v[t][r].qual_scores.begin(),sr_v[t][r].qual_scores.end());
	  if(remove_barcodes || sr_v[t][r].vec_bc_for.size()==0) //case we just want to print read once if multiple bc found.
	    print_read(sr_v[t][r].read_id+"_-",sr_v[t][r].rev_line,sr_v[t][r].qual_scores,sr_v[t][r].vec_bc_rev,
		       out_filename_prefix,split_file_by_barcode,found_barcodes,remove_barcodes);
	}
      }
    }
  }
  file.close();
  
  Rcpp::Rcout << "Number of reads processed: " << r_count << "\n";
  Rcpp::Rcout << "Number of reads where a barcode was found: " << bc_count << "\n";
  Rcpp::Rcout << "Number of reads where more than one barcode was found: " << multi_bc_count << "\n";
  Rcpp::Rcout << "All done!" << "\n";

  if(known_barcodes.size()>0){
    out_stat_file.close();
    return(0);
  }
  
  if(barcode_counts.size()==0)
    return(0);
  
  typedef std::pair<std::string, int> pair;
  vector<pair> bc_vec;
  copy(barcode_counts.begin(),barcode_counts.end(), back_inserter<vector<pair>>(bc_vec));
  sort(bc_vec.begin(), bc_vec.end(),[](const pair &l, const pair &r){
    if (l.second != r.second)
      return l.second > r.second;
    return l.first < r.first;
  });
  vector<int> hist(bc_vec[0].second);
  ofstream out_bc_file;
  out_bc_file.open(out_bc_filename);
  for (auto const &bc_pair: bc_vec){
    out_bc_file << bc_pair.first << "\t" << bc_pair.second << "\n";
    hist[bc_pair.second-1]++;
  }
  out_bc_file.close();

  cout << "Reads\tBarcodes" << "\n";
  for(int i=hist.size()-1; i>=0; i--)
    cout << i+1 << "\t" << hist[i] << "\n";
    
}
