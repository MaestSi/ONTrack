#
# Copyright 2019 Simone Maestri. All rights reserved.
# Simone Maestri <simone.maestri@univr.it>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

args = commandArgs(trailingOnly=TRUE)

if (args[1] == "-h" | args[1] == "--help") {
  cat("", sep = "\n")
  cat(paste0("Usage: Rscript ONTrack.R <home_dir> <fast5_dir> <sequencing_summary.txt>"), sep = "\n")
  cat(paste0("Note that config_ONTrack.R must be in the same directory of ONTrack.R"), sep = "\n")
  cat(paste0("<home_dir>: directory containing fastq and fasta files for each sample"), sep = "\n")
  cat(paste0("<fast5_dir>: directory containing raw fast5 files for nanopolish polishing, optional"), sep = "\n")
  cat(paste0("<sequencing_summary.txt>: sequencing summary file generated during base-calling, used to speed-up polishing, optional"), sep = "\n")
  stop(simpleError(sprintf("\r%s\r", paste(rep(" ", getOption("width")-1L), collapse=" "))))
}

if (length(args) == 1) {
  home_dir <- args[1]
  if (!dir.exists(home_dir)) {
    stop(paste0(home_dir, " directory does not exist!"))
  }
} else if (length(args) == 2) {
  home_dir <- args[1]
  fast5_dir <- args[2]
  if (!dir.exists(home_dir)) {
    stop(paste0(home_dir, " directory does not exist!"))
  } else if (!dir.exists(fast5_dir)) {
    stop(paste0(fast5_dir, " directory does not exist!"))
  }
} else if (length(args) == 3) {
  home_dir <- args[1]
  fast5_dir <- args[2]
  sequencing_summary <- args[3]
  if (!dir.exists(home_dir)) {
    stop(paste0(home_dir, " directory does not exist!"))
  } else if (!dir.exists(fast5_dir)) {
    stop(paste0(fast5_dir, " directory does not exist!"))
  }
  
} else {
  stop("At least one input argument must be provided")
}

PIPELINE_DIR <- dirname(strsplit(commandArgs(trailingOnly = FALSE)[4],"=")[[1]][2])

CONFIG_FILE <- paste0(PIPELINE_DIR, "/config_MinION_mobile_lab.R")
source(CONFIG_FILE)

if (!exists("sequencing_summary")) {
  seq_sum_flag <- 0
} else {
  seq_sum_flag <- 1
}

if (!exists("do_blast_flag")) {
  do_blast_flag <- 0
}

if (!exists("num_threads")) {
  num_threads <- 8
}

if (!exists("num_iterations")) {
  num_iterations <- 1
}

if (!exists("primers_length")) {
  primers_length <- 25
}

if (!exists("pair_strands_flag")) {
  pair_strands_flag <- 0
}

if (!exists("do_clustering_flag")) {
  do_clustering_flag <- 1
}

if (!exists("majority_rule_full_seq_flag")) {
  majority_rule_full_seq_flag <- 1
}

logfile <- paste0(home_dir, "/logfile.txt")

fasta_files <- list.files(path = home_dir, pattern = "BC\\d+\\.fasta", full.names = TRUE)
fastq_files <- paste0(home_dir, "/", gsub(pattern = "\\.fasta$", replacement = "\\.fastq", x = basename(fasta_files)))

if (length(fasta_files) > 0) {
  cat(text = paste0("Processing fasta files ", paste0(basename(fasta_files), collapse = ", ")), sep = "\n")
} else {
  stop(paste0("No fasta files in directory ", home_dir))
}

for (i in 1:length(fasta_files)) {
  sample_dir <- gsub(pattern = "\\.fasta", replacement = "", x = fasta_files[i])
  sample_name <- basename(sample_dir)
  dir.create(sample_dir)
  decont_fasta <- paste0(sample_dir, "/", sample_name, "_decont.fasta")
  decont_fastq <- paste0(sample_dir, "/", sample_name, "_decont.fastq")
  if (do_clustering_flag == 1) {
    system(paste0(DECONT, " ", fasta_files[i], " ", VSEARCH, " ", SEQTK))
    system(paste0("mv ", home_dir, "/decontam_tmp_", sample_name, " ", sample_dir))
    system(paste0("mv ", home_dir, "/", sample_name, "_decont.fasta ", sample_dir))
    system(paste0("mv ", home_dir, "/", sample_name, "_decont.fastq ", sample_dir))
  } else {
    system(paste0("cp ", home_dir, "/", sample_name, ".fasta ", sample_dir, "/", sample_name, "_decont.fasta"))
    system(paste0("cp ", home_dir, "/", sample_name, ".fastq ", sample_dir, "/", sample_name, "_decont.fastq"))
  }  
  num_reads_mac <- as.double(system(paste0("cat ", decont_fasta, " | grep \"^>\" | wc -l"), intern=TRUE))
  target_reads_contig <- 200
  target_reads_polishing <- 200
    
  if (num_reads_mac < target_reads_contig) {
    target_reads_contig <- num_reads_mac
    target_reads_polishing <- num_reads_mac
    if (do_clustering_flag == 1) {
      cat(text = paste0("WARNING: Only ", num_reads_mac, " reads available for sample ", sample_name, " after contaminants removal"), sep = "\n")
      cat(text = paste0("WARNING: Only ", num_reads_mac, " reads available for sample ", sample_name, " after contaminants removal"),  file = logfile, sep = "\n", append = TRUE)
    } else {
      cat(text = paste0("WARNING: Only ", num_reads_mac, " reads available for sample ", sample_name), sep = "\n")
      cat(text = paste0("WARNING: Only ", num_reads_mac, " reads available for sample ", sample_name),  file = logfile, sep = "\n", append = TRUE)
    }
  }

  plurality_value <- 0.15*target_reads_contig

  for (j in 1:num_iterations) {
    draft_contig_curr_tmp1 <- paste0(sample_dir, "/", sample_name, "_non_polished.contig_", j, "_tmp1.fasta")
    draft_contig_curr_tmp2 <- paste0(sample_dir, "/", sample_name, "_non_polished.contig_", j, "_tmp2.fasta")
    draft_contig_curr <- paste0(sample_dir, "/", sample_name, "_non_polished.contig_", j, ".fasta")
    draft_contig_curr_trimmed <- paste0(sample_dir, "/", sample_name, "_non_polished.contig_", j, "_trimmed.fasta")
    polished_contig_curr <- paste0(sample_dir, "/", sample_name, ".contigs_", j, ".fasta")
    polished_contig_curr_trimmed <- paste0(sample_dir, "/", sample_name, ".contigs_", j, "_trimmed.fasta")
    blast_results <- paste0(home_dir, "/", sample_name, ".blastn.txt")
    sequences <- readDNAStringSet(fasta_files[i], "fasta")
    ws <- width(sequences)
    amplicon_length <- ceiling(mean(ws))
    draft_reads_fq_curr <- paste0(sample_dir, "/", sample_name, "_draft_", target_reads_contig, "_reads_", j, ".fastq")
    draft_reads_fa_curr <- paste0(sample_dir, "/", sample_name, "_draft_", target_reads_contig, "_reads_", j, ".fasta")
    #target_bases_contig <- target_reads_contig*amplicon_length
    #system(paste0(FILTLONG, " --target_bases ", target_bases_contig, " --mean_q_weight 10 ", decont_fastq, " > ", draft_reads_fq_curr))
    seed_curr <- j*2-1
    system(paste0(SEQTK, " sample -s ", seed_curr , " ", decont_fastq, " ",  target_reads_contig, " > ", draft_reads_fq_curr))
    system(paste0(SEQTK, " seq -A ", draft_reads_fq_curr, " > ", draft_reads_fa_curr))
    mfa_file_curr <- gsub(pattern = "\\.fasta$", replacement = ".mfa", x = draft_reads_fa_curr)
    #system(paste0(MAFFT ,"-linsi --thread ", num_threads, " --threadit 0 --adjustdirectionaccurately ", draft_reads_fa_curr, " > ", mfa_file_curr)) #threadit 0 if you need reproducible results
    system(paste0(MAFFT ,"-linsi --thread ", num_threads, " --adjustdirectionaccurately ", draft_reads_fa_curr, " > ", mfa_file_curr))
    system(paste0(CONS, " -sequence ", mfa_file_curr, " -plurality ", plurality_value, " -outseq ", draft_contig_curr_tmp1))
    system(paste0("sed 's/[nN]//g' ", draft_contig_curr_tmp1, " > ", draft_contig_curr_tmp2))
    DNAStringSet_obj <- readDNAStringSet(draft_contig_curr_tmp2, "fasta")
    DNAStringSet_obj_renamed <- DNAStringSet_obj
    original_headers <- names(DNAStringSet_obj)
    sequences <- seq(DNAStringSet_obj)
    new_headers <- c()
    prefix <- "contig"
    new_headers <- paste0(prefix, "_", j)
    names(DNAStringSet_obj_renamed) <- new_headers
    writeXStringSet(x = DNAStringSet_obj_renamed, filepath = draft_contig_curr, format = "fasta", width = 20000)
    system(paste0(SEQTK, " trimfq ", draft_contig_curr, " -b ", primers_length, " -e ", primers_length, " > ", draft_contig_curr_trimmed))
    if (exists("fast5_dir") && pair_strands_flag != 1) {
      qual_filter <- 0
      reads_polishing_fq_curr <- paste0(sample_dir, "/", sample_name, "_polishing_", target_reads_polishing, "_reads_", j, ".fastq")
      reads_polishing_fa_curr <- paste0(sample_dir, "/", sample_name, "_polishing_", target_reads_polishing, "_reads_", j, ".fasta")  
      #target_bases_polishing <- target_reads_polishing*amplicon_length
      #system(paste0(FILTLONG, " --target_bases ", target_bases_polishing, " --mean_q_weight 10 ", decont_fastq, " > ", reads_polishing_fq_curr))
      seed_polishing_curr <- j*2 
      system(paste0(SEQTK, " sample -s ", seed_polishing_curr, " ", decont_fastq, " ", target_reads_polishing, " > ", reads_polishing_fq_curr))
      system(paste0(SEQTK, " seq -A ", reads_polishing_fq_curr, " > ", reads_polishing_fa_curr))
      bam_file_curr <- paste0(sample_dir, "/", sample_name, "_", j, ".bam")
      system(paste0(MINIMAP2, " -ax map-ont ", draft_contig_curr, " ", reads_polishing_fa_curr, " | ", SAMTOOLS, " view -h -q 55 -F 2048 | " , SAMTOOLS, " sort -o ", bam_file_curr, " -T reads.tmp"))
      system(paste0(SAMTOOLS," index ", bam_file_curr))
      cat(text = paste0("Running nanopolish for sample ", sample_name, " - iteration num. ", j, " out of ", num_iterations), sep = "\n")
      cat(text = "Indexing reads", sep = "\n")
      if (seq_sum_flag == 1) {
        system(paste0(NANOPOLISH, " index -d ", fast5_dir, " -s ", sequencing_summary, " ", reads_polishing_fa_curr))
      } else {
        system(paste0(NANOPOLISH, " index -d ", fast5_dir, " ", reads_polishing_fa_curr))
      }
      cat(text = paste0("Running nanopolish for consensus polishing of sample ", sample_name), sep = "\n")
      output_vcf_curr <- paste0(sample_dir, "/", "nanopolish_output_", j, ".vcf")
      output_vcf_curr_filtered <- paste0(sample_dir, "/", "nanopolish_output_", j, "_filtered.vcf")
      system(paste0(NANOPOLISH, " variants --consensus --reads ", reads_polishing_fa_curr, " --bam ", bam_file_curr, " --genome ", draft_contig_curr, " -p 1 --threads ", num_threads, " -o ", output_vcf_curr))
      system(paste0("cat ", output_vcf_curr, " | grep \"^#\" > ", output_vcf_curr_filtered))
      system(paste0("cat ", output_vcf_curr, " | grep -v \"^#\" | awk '$6 > ", qual_filter, " {print}' >> ", output_vcf_curr_filtered))
      system(paste0(NANOPOLISH, " vcf2fasta -g ", draft_contig_curr, " ", output_vcf_curr_filtered, " > ", polished_contig_curr))
      system(paste0(SEQTK, " trimfq ", polished_contig_curr, " -b ", primers_length, " -e ", primers_length, " > ", polished_contig_curr_trimmed))
    } else {
      system(paste0("cp ", draft_contig_curr, " ", polished_contig_curr))
    }
  }
  final_contig <- paste0(home_dir, "/", sample_name, ".contigs.fasta")
  if (num_iterations > 1) {
    aggregate_contigs <- paste0(sample_dir, "/", sample_name, ".contigs_", num_iterations, "iterations.fasta")
    polished_contig_agg <- list.files(path = sample_dir, pattern = paste0(sample_name, "\\.contigs_\\d*\\.fasta"), full.names = TRUE)
    system(paste0("cat ", paste0(polished_contig_agg, collapse = " "), " > ", aggregate_contigs))
    mfa_file_agg <- gsub(pattern = "\\.fasta$", replacement = ".mfa", x = aggregate_contigs)
    system(paste0(MAFFT, " --op 0 --thread ", num_threads, " --adjustdirectionaccurately ", aggregate_contigs, " > ", mfa_file_agg))
    if (majority_rule_full_seq_flag == 1) {
      derep_file_agg <- gsub(pattern = "\\.fasta$", replacement = "_derep.fasta", x = aggregate_contigs)
      most_freq_seq <- paste0(sample_dir, "/", sample_name, ".contigs_tmp1.fasta")
      aggregate_contigs_onestrand <- gsub(pattern = "\\.fasta$", replacement = "_onestrand.fasta", x = aggregate_contigs)
      system(paste0("cat ", mfa_file_agg, " | sed 's/-//g' > ", aggregate_contigs_onestrand))
      system(paste0(VSEARCH, " --derep_prefix ", aggregate_contigs_onestrand , " --output ", derep_file_agg, " --sizeout --fasta_width 0"))
      num_supp_seq <- as.double(system(paste0("head -n2 ", derep_file_agg, " | grep \"^>\" | sed 's/.*size=//'"), intern=TRUE))
      if (num_supp_seq == 1) {
        system(paste0("head -n2 ", aggregate_contigs, " > ", most_freq_seq))
        cat(text = paste0("WARNING: No common sequence identified for sample ", sample_name, "; the sequence from the 1st iteration was selected"), sep = "\n")
        cat(text = paste0("WARNING: No common sequence identified for sample ", sample_name, "; the sequence from the 1st iteration was selected"),  file = logfile, sep = "\n", append = TRUE)
      } else {
        system(paste0("head -n2 ", derep_file_agg, " > ", most_freq_seq))
      }
      system(paste0(SEQTK, " trimfq ", most_freq_seq, " -b ", primers_length, " -e ", primers_length, " > ", final_contig))
    } else {
      plurality_value_agg <- 0.5*num_iterations
      final_contig_tmp1 <- paste0(sample_dir, "/", sample_name, ".contigs_tmp1.fasta")
      final_contig_tmp2 <- paste0(sample_dir, "/", sample_name, ".contigs_tmp2.fasta")
      final_contig_tmp3 <- paste0(sample_dir, "/", sample_name, ".contigs_tmp3.fasta")
      system(paste0(CONS, " -sequence ", mfa_file_agg, " -plurality ", plurality_value_agg, " -outseq ", final_contig_tmp1))
      system(paste0("sed 's/[nN]//g' ", final_contig_tmp1, " > ", final_contig_tmp2))
      DNAStringSet_obj <- readDNAStringSet(final_contig_tmp2, "fasta")
      DNAStringSet_obj_renamed <- DNAStringSet_obj
      original_headers <- names(DNAStringSet_obj)
      sequences <- seq(DNAStringSet_obj)
      new_header <- paste0("contig_", num_iterations, "_iterations")
      names(DNAStringSet_obj_renamed) <- new_header
      writeXStringSet(x = DNAStringSet_obj_renamed, filepath = final_contig_tmp3, format = "fasta", width = 20000)
      system(paste0(SEQTK, " trimfq ", final_contig_tmp3, " -b ", primers_length, " -e ", primers_length, " > ", final_contig))
    }
  } else {
    system(paste0(SEQTK, " trimfq ", polished_contig_curr, " -b ", primers_length, " -e ", primers_length, " > ", final_contig))
  }
  
  ONtoBAR_compatibility <- 0
  if (do_blast_flag ==1 ) {
    cat(text = paste0("BLASTing consensus for sample ", sample_name, " and saving results to ", basename(blast_results)), sep = "\n")
    if (ONtoBAR_compatibility != 1) {
      system(paste0(BLASTN, " -num_threads ", num_threads, " -db ", NTDB, " -query ", final_contig, " > ", blast_results))
    } else {
      blast_results_alt_sort <- paste0(sample_dir, "/", sample_name, "_sorted_by_perc_id.blastn.txt")
      blast_results_raw <- paste0(sample_dir, "/", sample_name, ".blastn_output.txt")
      system(paste0(BLASTN, " -db ", NTDB, " -query ", final_contig, " > ", blast_results_raw))
      cat(text = paste0("BLASTing consensus for sample ", sample_name, " and saving results to ", basename(blast_results)), sep = "\n")
      system(paste0(BLASTN, " -num_threads ", num_threads, " -db ", NTDB, " -query ", final_contig, " > ", blast_results_raw))
      input_file <- file(blast_results_raw, open = "r")
      input <- readLines(input_file)
      close(input_file)
      gb_vec <- vector(mode = "character")
      name_vec <- vector(mode = "character")
      alignment_id_vec <- vector(mode = "character")
      perc_id_vec <- vector(mode = "integer")
      score_vec <- vector(mode = "integer")
      counter <- 0
      for (j in 1:length(input)) {
        if (length(grep(pattern = "^>", x = input[j])) == 1) {
          counter <- counter + 1
          match_id <- input[j]
          match_id_split <- strsplit(x = match_id, split = "\\|")[[1]]
          gb <- match_id_split[2]
          gb_vec[counter] <- gb
          full_name <- match_id_split[3]
          name <- substr(x = full_name, start = 1, stop = 50)
          name_vec[counter] <- name
        } else if (length(grep(pattern = "^\\s*Score\\s*=\\s*\\d+\\s+bits.+$", x = input[j])) == 1) {
          score <- sub(pattern = "^\\s*Score\\s*=\\s*(\\d+)\\s+bits.+$", "\\1", input[j])
          score_vec[counter] <- score
        } else if (length(grep(pattern = "^\\s*Identities\\s*=\\s+(.+),\\s+.+$", x = input[j])) == 1) {
          alignment_id <-  sub(pattern = "^\\s*Identities\\s*=\\s+(.+),\\s+.+$", "\\1", x = input[j])
          alignment_id_vec[counter] <- alignment_id
          perc_id <- sub(pattern = "^\\s*Identities\\s*=\\s+\\d+\\/\\d+\\s*\\((\\d+)%\\),\\s+.+$", "\\1", input[j])
          perc_id_vec[counter] <- perc_id
        }
      }
      ind_sort <- sort.int(score_vec, decreasing = TRUE, index.return = TRUE)$ix
      ind_sort_alt <- sort.int(perc_id_vec, decreasing = TRUE, index.return = TRUE)$ix
      output_data <- data.frame(gb_vec[ind_sort], name_vec[ind_sort], alignment_id_vec[ind_sort], score_vec[ind_sort], perc_id_vec[ind_sort])
      write.table(x = output_data, file = blast_results, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)  
      output_data_alt_sort <- data.frame(gb_vec[ind_sort_alt], name_vec[ind_sort_alt], alignment_id_vec[ind_sort_alt], score_vec[ind_sort_alt], perc_id_vec[ind_sort_alt])
      write.table(x = output_data_alt_sort, file = blast_results_alt_sort, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
    }
  }
}
