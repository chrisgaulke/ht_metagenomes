#Title: combine_data.R
#Author(s): Christopher A Gaulke
#Date: 2020-03-23
#Project: HT Metagenomes

#This script will combine plexwell replicates into one single replicate

# Environment: Packages and options ----------------------------------------
options("stringsAsFactors"=F)

# # Environment: Functions ---------------------------------------------------------------

make_tres <- function(vec) {
  #make sure all numbers are three digits
  vec <- vec
  for( i in 1:length(vec)){
    if(nchar(vec[i]) < 4){
      vec[i] <- gsub(pattern = "s", x = vec[i], replacement = "s0")
    }
  }
  return(vec)
}

# Data: Import Kraken2 -------------------------------------------------------------

original_kraken <- read.table("/Users/gaulkec/Chris/dev/R_projects/ht_metagenomes/data/corrected_analysis/kraken2_merged.txt",
                              sep = "\t",
                              header = T,
                              row.names = NULL,
                              comment.char = "",
                              quote = ""
                              )

# Data: Import Metaphlan -------------------------------------------------------------

coral.metaphlan <- read.table(
  "/Users/gaulkec/Chris/dev/R_projects/ht_metagenomes/data/corrected_analysis/combined_tables/coral_combined_metaphlan_families_cpm.tab",
  sep = "\t",
  header = T,
  comment.char = ""
)

soil.metaphlan <- read.table(
  "/Users/gaulkec/Chris/dev/R_projects/ht_metagenomes/data/corrected_analysis/combined_tables/soil_combined_metaphlan_families_cpm.tab",
  sep = "\t",
  header = T,
  comment.char = ""
)

fecal.metaphlan <- read.table(
  "/Users/gaulkec/Chris/dev/R_projects/ht_metagenomes/data/corrected_analysis/combined_tables/feces_combined_metaphlan_families_cpm.tab",
  sep = "\t",
  header = T,
  comment.char = ""
)

mock.metaphlan <- read.table(
  "/Users/gaulkec/Chris/dev/R_projects/ht_metagenomes/data/corrected_analysis/combined_tables/mock_combined_metaphlan_families_cpm.tab",
  sep = "\t",
  header = T,
  comment.char = ""
)

# Data: Import Gene Data -------------------------------------------------------------

coral.gene <- read.table(
  "/Users/gaulkec/Chris/dev/R_projects/ht_metagenomes/data/corrected_analysis/combined_tables/coral_combined_gene_families_cpm.tab",
  sep = "\t",
  header = T,
  comment.char = ""
)

soil.gene <- read.table(
  "/Users/gaulkec/Chris/dev/R_projects/ht_metagenomes/data/corrected_analysis/combined_tables/soil_combined_gene_families_cpm.tab",
  sep = "\t",
  header = T,
  comment.char = ""
)

fecal.gene <- read.table(
  "/Users/gaulkec/Chris/dev/R_projects/ht_metagenomes/data/corrected_analysis/combined_tables/feces_combined_gene_families_cpm.tab",
  sep = "\t",
  header = T,
  comment.char = ""
)

mock.gene <- read.table(
  "/Users/gaulkec/Chris/dev/R_projects/ht_metagenomes/data/corrected_analysis/combined_tables/mock_combined_gene_families_cpm.tab",
  sep = "\t",
  header = T,
  comment.char = ""
)

# Data: Import Metadata -----------------------------------------------------

original_metadata <- read.table("/Users/gaulkec/Chris/dev/R_projects/ht_metagenomes/metadata/2020_01_16_metadata_all.txt",
                                   sep = "\t",
                                   header = T
)

# Data: Import Sequence Statistics ------------------------------------------

htmeta_seq_stats.df <- read.table("/Users/gaulkec/Chris/dev/R_projects/ht_metagenomes/data/corrected_analysis/clean_multiqc/output/R1/multiqc_data/multiqc_general_stats.txt",
                                  sep = "\t",
                                  header = T,
                                  row.names = 1,
                                  stringsAsFactors = F)

# Data: Import Multiqc ------------------------------------------

multiqc_raw.df <- read.table("/Users/gaulkec/Chris/dev/R_projects/ht_metagenomes/data/corrected_analysis/raw_multiqc/R1/multiqc_data/multiqc_general_stats.txt",
                                  sep = "\t",
                                  header = T,
                                  row.names = 1,
                                  stringsAsFactors = F)

multiqc_clean.df <- read.table("/Users/gaulkec/Chris/dev/R_projects/ht_metagenomes/data/corrected_analysis/clean_multiqc/output/R1/multiqc_data/multiqc_general_stats.txt",
                                  sep = "\t",
                                  header = T,
                                  row.names = 1,
                                  stringsAsFactors = F)

# Modify: Metadata --------------------------------------------------------

original.names <- paste0("s", original_metadata[,"Sample.ID"])
original.names <- make_tres(original.names)
original_metadata$Sample.ID <- original.names
rownames(original_metadata) <- original_metadata$Sample.ID

plexwell.metadata <-
  original_metadata[which(original_metadata$Prep.Kit == "PlexWell96"),]

#gather information for filtering
plexwell.names <- plexwell.metadata$Sample.ID
sample_types <- unique(original_metadata$Sample.Type)

# Modify: Kraken2 --------------------------------------------------------

#remove the plexwell samples
combined_kraken <- original_kraken[,-which(colnames(original_kraken) %in% plexwell.names)]

#now add back averaged plexwells one at a time
for(i in 1:length(sample_types)){
  #subset the metadata table to allow for easy sample name collection
  t.metadata <- plexwell.metadata[which(plexwell.metadata$Sample.Type == sample_types[i]),]
  for(j in 1:5) {
    #get sample names and average each taxa
    i.vec <- seq(from = j, to = 30, by = 5)
    s.vec <- t.metadata[i.vec, "Sample.ID"]
    i.name <- paste0(sample_types[i],"_plexWell96_", j )
    combined_kraken[,i.name] <- rowMeans(original_kraken[,s.vec])
  }
}

# Modify: Make Combined Metadata ---------------------------------------------

#combined data frame is complete now update the metadata

t_metadata <- original_metadata[-which(rownames(original_metadata) %in% plexwell.names),]

combined_metadata <- data.frame(Sample.ID  = c(t_metadata$Sample.ID,
                                            colnames(combined_kraken)[82:101]),
                                Prep.Kit   = c(t_metadata$Prep.Kit,
                                  rep("plexWell96", times = 20)
                                  ),
                                Input..ng. = c(t_metadata$Input..ng.,
                                  rep(c(1.0,1.0,1.0,0.5,5.0), times = 4)
                                  ),
                                Sample.Type = c(t_metadata$Sample.Type,
                                  rep(c("Coral", "Soil", "Feces", "Mock"),
                                      times = c(5,5,5,5))
                                  )
                                )

# Modify: Multiqc ------------------------------------------

combined_multiqc <- multiqc_raw.df[,"FastQC_mqc.generalstats.fastqc.total_sequences", drop = F]
colnames(combined_multiqc) <- "raw_seqs"
combined_multiqc$clean_seqs <- multiqc_clean.df$FastQC_mqc.generalstats.fastqc.total_sequences
combined_multiqc$filtered <- combined_multiqc$raw_seqs - combined_multiqc$clean_seqs
combined_multiqc$per_retained <- (combined_multiqc$clean_seqs / combined_multiqc$raw_seqs) * 100
combined_multiqc$per_filtered <- (1 - (combined_multiqc$clean_seqs / combined_multiqc$raw_seqs)) * 100

#mod rownames to be consistent with all other data

multiqc.names <- sapply(rownames(combined_multiqc),
                        FUN = function(x){
                          strsplit(x = x, split = "\\-")[[1]][2]
                        }
)

names(multiqc.names) <- NULL

#adjust rownames
rownames(combined_multiqc) <- multiqc.names

combined_multiqc <- t(combined_multiqc)

#pool seqwell
tcombined_multiqc <- combined_multiqc[,-which(colnames(combined_multiqc) %in% plexwell.names)]
tcombined_multiqc <- as.data.frame(combined_multiqc)

#now add back averaged plexwells one at a time
for(i in 1:length(sample_types)){
  #subset the metadata table to allow for easy sample name collection
  t.metadata <- plexwell.metadata[which(plexwell.metadata$Sample.Type == sample_types[i]),]

  for(j in 1:5) {
    #get sample names and average each taxa
    i.vec <- seq(from = j, to = 30, by = 5)
    s.vec <- t.metadata[i.vec, "Sample.ID"]
    i.name <- paste0(sample_types[i],"_plexWell96_", j )
    tcombined_multiqc[,i.name] <- rowMeans(combined_multiqc[,s.vec])
  }
}

combined_multiqc <- as.data.frame(t(tcombined_multiqc))

#now add some metadata
combined_multiqc$sample_id <- rownames(combined_multiqc)
combined_multiqc <- combined_multiqc[combined_metadata$Sample.ID,]
combined_multiqc$type <- combined_metadata$Sample.Type
combined_multiqc$prep <- combined_metadata$Prep.Kit
combined_multiqc$conc <- combined_metadata$Input..ng.

# Modify: Sequence Statistics --------------------------------------------

#clean up ridic colnames
colnames(htmeta_seq_stats.df) <- c("Dups",
                                   "Length",
                                   "M.Seqs",
                                   "Failed",
                                   "GC")

#make total seqs into total seqs (millions)
htmeta_seq_stats.df$M.Seqs <- htmeta_seq_stats.df$M.Seqs/1E6

#mod rownames to be consistent with all other data

seqstat.names <- sapply(rownames(htmeta_seq_stats.df),
                        FUN = function(x){
                          strsplit(x = x, split = "\\-")[[1]][2]
                        }
)

names(seqstat.names) <- NULL

#adjust rownames
rownames(htmeta_seq_stats.df) <- seqstat.names

htmeta_seq_stats.df <- t(htmeta_seq_stats.df)

#pool seqwell
combined.seq_stats <- htmeta_seq_stats.df[,-which(colnames(htmeta_seq_stats.df) %in% plexwell.names)]
combined.seq_stats <- as.data.frame(combined.seq_stats)

#now add back averaged plexwells one at a time
for(i in 1:length(sample_types)){
  #subset the metadata table to allow for easy sample name collection
  t.metadata <- plexwell.metadata[which(plexwell.metadata$Sample.Type == sample_types[i]),]

  for(j in 1:5) {
    #get sample names and average each taxa
    i.vec <- seq(from = j, to = 30, by = 5)
    s.vec <- t.metadata[i.vec, "Sample.ID"]
    i.name <- paste0(sample_types[i],"_plexWell96_", j )
    combined.seq_stats[,i.name] <- rowMeans(htmeta_seq_stats.df[,s.vec])
  }
}

combined.seq_stats <- as.data.frame(t(combined.seq_stats))

#now add some metadata
combined.seq_stats$sample_id <- rownames(combined.seq_stats)
combined.seq_stats <- combined.seq_stats[combined_metadata$Sample.ID,]
combined.seq_stats$type <- combined_metadata$Sample.Type
combined.seq_stats$prep <- combined_metadata$Prep.Kit
combined.seq_stats$conc <- combined_metadata$Input..ng.
combined.seq_stats$per_retained <- combined_multiqc$per_retained
combined.seq_stats$per_filtered <- combined_multiqc$per_filtered


# Modify: Metaphlan ---------------------------------------------------------------

#merge
merge.metaphlan <- merge(mock.metaphlan, coral.metaphlan, by = 1, all = T)
merge.metaphlan <- merge(merge.metaphlan, fecal.metaphlan, by = 1, all = T)
merge.metaphlan <- merge(merge.metaphlan, soil.metaphlan, by = 1, all = T)

merge_metaphlan.names <- sapply(colnames(merge.metaphlan),
                                FUN = function(x){
                                  strsplit(x = x, split = "\\.")[[1]][2]
                                  }
                                )

names(merge_metaphlan.names) <- NULL

colnames(merge.metaphlan) <- merge_metaphlan.names
merge.metaphlan[is.na(merge.metaphlan)] <- 0

#remove the plexwell samples
combined.metaphlan <- merge.metaphlan[,-which(colnames(merge.metaphlan) %in% plexwell.names)]

#now add back averaged plexwells one at a time
for(i in 1:length(sample_types)){
  #subset the metadata table to allow for easy sample name collection
  t.metadata <- plexwell.metadata[which(plexwell.metadata$Sample.Type == sample_types[i]),]
  for(j in 1:5) {
    #get sample names and average each taxa
    i.vec <- seq(from = j, to = 30, by = 5)
    s.vec <- t.metadata[i.vec, "Sample.ID"]
    i.name <- paste0(sample_types[i],"_plexWell96_", j )
    combined.metaphlan[,i.name] <- rowMeans(merge.metaphlan[,s.vec])
  }
}

# Modify: Gene ---------------------------------------------------------------

merge.gene <- merge(mock.gene, coral.gene, by = 1, all = T)
merge.gene <- merge(merge.gene, fecal.gene, by = 1, all = T)
merge.gene <- merge(merge.gene, soil.gene, by = 1, all = T)

colnames(merge.gene)[1] <- "X.GeneID"
merge_gene.names <- sapply(colnames(merge.gene),
                           FUN = function(x){
                             strsplit(x = x, split = "\\.")[[1]][2]
                           }
)

names(merge_gene.names) <- NULL

colnames(merge.gene) <- merge_gene.names
merge.gene[is.na(merge.gene)] <- 0

#remove the plexwell samples
combined.gene <- merge.gene[,-which(colnames(merge.gene) %in% plexwell.names)]

#now add back averaged plexwells one at a time
for(i in 1:length(sample_types)){
  #subset the metadata table to allow for easy sample name collection
  t.metadata <- plexwell.metadata[which(plexwell.metadata$Sample.Type == sample_types[i]),]
  for(j in 1:5) {
    #get sample names and average each taxa
    i.vec <- seq(from = j, to = 30, by = 5)
    s.vec <- t.metadata[i.vec, "Sample.ID"]
    i.name <- paste0(sample_types[i],"_plexWell96_", j )
    combined.gene[,i.name] <- rowMeans(merge.gene[,s.vec])
  }
}

# Write: Combined Files ------------------------------------------------

write.table(file = "/Users/gaulkec/Chris/dev/R_projects/ht_metagenomes/analysis/flat_files/combined_metadata.txt",
            x = combined_metadata,
            sep = "\t",
            row.names = F,
            col.names = T,
            quote = F
)

write.table(file = "/Users/gaulkec/Chris/dev/R_projects/ht_metagenomes/analysis/flat_files/combined_kraken2.txt",
            x = combined_kraken,
            sep = "\t",
            row.names = F,
            col.names = T,
            quote = F
)

write.table(file = "/Users/gaulkec/Chris/dev/R_projects/ht_metagenomes/analysis/flat_files/combined_metaphlan.txt",
            x = combined.metaphlan,
            sep = "\t",
            row.names = F,
            col.names = T,
            quote = F
)

write.table(file = "/Users/gaulkec/Chris/dev/R_projects/ht_metagenomes/analysis/flat_files/combined_gene.txt",
            x = combined.gene,
            sep = "\t",
            row.names = F,
            col.names = T,
            quote = F
)


write.table(file = "/Users/gaulkec/Chris/dev/R_projects/ht_metagenomes/analysis/flat_files/cleaned_seq_stats.txt",
            x = combined.seq_stats,
            sep = "\t",
            row.names = T,
            col.names = T,
            quote = F
)

write.table(file = "/Users/gaulkec/Chris/dev/R_projects/ht_metagenomes/analysis/flat_files/multiqc_stats.txt",
            x = combined_multiqc,
            sep = "\t",
            row.names = T,
            col.names = T,
            quote = F
)
