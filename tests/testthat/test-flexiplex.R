# functions to generate random bases
sample_N <- function(n, bases = c("A", "T", "C", "G")) {
  if (length(n) > 1) {
    n <- sample(n, 1)
  }
  paste0(
    sample(bases, n, replace = T),
    collapse = ""
  )
}

sample_Y <- function(n) {
  sample_N(n = n, bases = c("C", "T"))
}

sample_R <- function(n) {
  sample_N(n = n, bases = c("A", "G"))
}

edit_seq <- function(sequence = "CTACACGACGCTCTTCCGATCT", sites = 1, p_indel = 0.2) {
  edited_sequence <- strsplit(sequence, "")[[1]]
  while (sites > 0) {
    if (runif(1) < p_indel) { # indel
      if (runif(1) < 0.5) { # del
        edited_sequence <- edited_sequence[-sample(seq_along(edited_sequence), 1)]
      } else { # insertion
        edited_sequence <- append(edited_sequence, sample_N(1), sample(seq_along(edited_sequence), 1))
      }
    } else {
      edited_sequence[sample(seq_along(edited_sequence), 1)] <- sample_N(1)
    }
    sites <- sites - 1
  }
  return(paste0(edited_sequence, collapse = ""))
}

sample_BC <- function(len, n) {
  sapply(rep(len, n), sample_N)
}

# -x CTACACGACGCTCTTCCGATCT -b NNNNNNNNNNNNNNNN -x YR -u NNNNNNNNNNNN -x TTTTTTTTT
sample_BC_example_1 <- function(len = 16, n = 100){
  return(sample_BC(len, n))
}
example_1 <- function(BCs, ...) {
  BC <- sample(BCs, 1)
  UMI <- sample_N(12)
  read <- paste0(
    c(
      sample_N(5:10),
      edit_seq(sequence = "CTACACGACGCTCTTCCGATCT", sites = 2, p_indel = 0.2),
      edit_seq(sequence = BC, sites = 1, p_indel = 0.2), sample_Y(1), sample_R(1), UMI,
      edit_seq(sequence = "TTTTTTTTTTTTT", sites = 2, p_indel = 0.2),
      sample_N(5:10)
    ),
    collapse = ""
  )
  return(c(read, BC, UMI))
}

write_example <- function(example_fn, bc_fn, filename_reads, filename_bc) {
  BCs <- bc_fn()
  df <- sapply(1:1000, function(x){example_fn(BCs = BCs)}) |>
    t() |>
    as.data.frame() |>
    setNames(c("read", "BC", "UMI"))

  fp <- file(filename_reads, open = "wt")
  for (i in seq_len(nrow(df))) {
    writeLines(
      c(
        paste0(c(">", df$BC[i], "_", df$UMI[i]), collapse = ""),
        df$read[i]
      ),
      fp
    )
  }
  close(fp)

  # also write a BC allow list for flexiplex
  fp <- file(filename_bc, open = "wt")
  writeLines(BCs, fp)
  close(fp)
}

test_that("flexiplex_performance", {
  outdir <- tempdir()
  fasta <- file.path(outdir, "test.fasta")
  bc <- file.path(outdir, "bc.txt")
  write_example(example_1, sample_BC_example_1, fasta, bc)

  if (file.exists(file.path(outdir, "stats.tsv"))) {
    file.remove(file.path(outdir, "stats.tsv"))
  }
  if (file.exists("flexiplex_reads_barcodes.txt")) {
    file.remove("flexiplex_reads_barcodes.txt")
  }

  find_barcode(
    max_bc_editdistance = 2, max_flank_editdistance = 6,
    fastq = fasta,
    barcodes_file = bc,
    reads_out = file.path(outdir, "out.fq"),
    stats_out = file.path(outdir, "stats.tsv"),
    threads = 1, pattern = c(
      primer = "CTACACGACGCTCTTCCGATCT",
      BC = paste0(rep("N", 16), collapse = ""),
      yr = "YR",
      UMI = paste0(rep("N", 12), collapse = ""),
      polyT = paste0(rep("T", 9), collapse = "")
    )
  )

  res <- read.delim(file.path(outdir, "stats.tsv"))
  file.remove(file.path(outdir, "stats.tsv"))
  cat("UMI performance:")
  print(table(res$UMI == gsub("^.*_", "", res$Read)))
  cat("BC performance:")
  print(table(res$CellBarcode == gsub("_.*$", "", res$Read)))

  expect_true(all(res$CellBarcode == gsub("_.*$", "", res$Read)))

  # binary tests
  expect_true(file.exists(test_path("flexiplex-bin")))
  system2(
    command = file.path(".", test_path("flexiplex-bin")),
    args = c(
      "-k", bc,
      "-x CTACACGACGCTCTTCCGATCT -b '????????????????' -x YR -u '????????????' -x TTTTTTTTT",
      "-f 6 -e 2",
      fasta
    ),
    stdout = FALSE # discard reads
  )
  res_bin <- read.delim("flexiplex_reads_barcodes.txt")
  file.remove("flexiplex_reads_barcodes.txt")
  cat("binary file UMI performance:")
  print(table(res_bin$UMI == gsub("^.*_", "", res_bin$Read)))
  cat("binary file BC performance:")
  print(table(res_bin$CellBarcode == gsub("_.*$", "", res_bin$Read)))
  expect_equal(
    sum(res_bin$UMI == gsub("^.*_", "", res_bin$Read)),
    sum(res$UMI == gsub("^.*_", "", res$Read))
  )
  expect_true(all(res_bin$CellBarcode == gsub("_.*$", "", res_bin$Read)))
})

# " >ACAACCGGATCGGGGC_TATCTTCTAGAT                                                    "
# "                                               | UMI      |                        "
# " CCTTTCTACACGACGCTCTTCCGATCGCACAACCGGATCGGGGCGTTATCTTCTAGATTTTTTTTTTAAGAAGTCC      "
# "      |                    /                  | reported |                         "
# "      CTACACGACGCTCTTCCGATCT  T edited to GC,                                      "
# "                              but aligment probably did not add a insertion here   "
# "                              hence the rest of the aligment shifted 1 nt          "
