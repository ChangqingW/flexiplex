# flexiplex
The Flexible Demultiplexer

This is the Rcpp version of flexiplex.  
Please vistion the main repository at [DavidsonGroup/flexiplex](https://github.com/DavidsonGroup/flexiplex) for the standard flexiplex release.

Usage:  
```
flexiplex::flexiplex(
  reads_in="in.fa", # input
  barcodes_file="bc.txt",  # barcodes list
  reads_out="r-fq.out.fasta",  # \
  stats_out="r-stats.txt",     # outputs
  bc_out="r-bc.txt",           # /
  bc_as_readid=T, 
  max_flank_editdistance=8, 
  max_bc_editdistance=2, 
  pattern=list(),  # wip
  n_threads=2)
```
