#Title: ht_metagenomes_mock_analysis.R
#Author(s): Christopher A Gaulke
#Date: 2020-02-24
#Project: HT Metagenomes

# Environment: Packages and options ----------------------------------------
library(ggplot2)
library(vegan)
library(reshape2)
options("stringsAsFactors"=F)

# Environment: Function ---------------------------------------------------

metaphlan_parser <- function(df) {

  #Description:

  # function to parse kraken report into a form that will be easy to work with in R
  # (long overdue).

  #Usage:
  #kraken_parser(df)

  #Arguments
  #df: kraken 1.0 mpa report output

  #Details:
  # The input report should have only 1 counts column. If more than one sample
  # is being interogated then additional columns will need to be added to the
  # dataframe. This feature is easy to write and I will get to it soon.

  #Value:
  #a dataframe with the following columns:
  #full: full tax string
  #k: kingdom
  #p: phylum
  #c: class
  #o: order
  #f: family
  #g: genus
  #s: species
  #t: strain
  #level: lowest tax level classified (e.g., species)
  #count: the number of classified reads that were binned in this taxa label


  rows <- nrow(df)
  ndf <- as.data.frame(matrix(NA, ncol = ncol(df) + 9, nrow = rows))
  colnames(ndf) <- c("full",
                     "k",
                     "p",
                     "c",
                     "o",
                     "f",
                     "g",
                     "s",
                     "t",
                     "level",
                     colnames(df)[2:ncol(df)])


  for( i in 1:nrow(df)){
    #print(i) #for debug'n
    x <- unlist(strsplit(x = df[i,1], split = "\\||__"))
    #print(x) #for debug'n
    if(x[1] == "unclassified"){
      ndf[i,1:10] <- "unclassified"
      ndf[i,11:ncol(ndf)]  <- df[i,2:ncol(df)] #don't want the sample id col
      next
    }else{
      y <- x[seq(from=2, to=length(x), by = 2)] #these are tax labels
      yy <- x[seq(from=1, to=length(x), by = 2)] #these are colnames
      #print(y)
      ndf[i,yy] <- y
      ndf[i,"full"]  <- df[i,1]
      ndf[i,"level"] <- x[length(x)-1]
      ndf[i,11:ncol(ndf)]  <- df[i,2:ncol(df)] #don't want the sample id col
    }
  }

  return(ndf)

}


# Data: Import data -------------------------------------------------------
mock.metaphlan <- read.table(
  "/Users/gaulkec/Chris/dev/R_projects/ht_metagenomes/analysis/flat_files/combined_metaphlan.txt",
  sep = "\t",
  header = T,
  comment.char = ""
)

mock_metadata.df <- read.table(
    "/Users/gaulkec/Chris/dev/R_projects/ht_metagenomes/analysis/flat_files/combined_metadata.txt",
  sep = "\t",
  header = T,
  comment.char = ""
)
rownames(mock_metadata.df) <- mock_metadata.df$Sample.ID

#get sample names for the mock community
mock.names <-  rownames(mock_metadata.df[which(mock_metadata.df$Sample.Type == "Mock"),])

mock_metadata.df <- mock_metadata.df[which(rownames(mock_metadata.df) %in% mock.names),]
mock_metadata.df["truth",] <- c("truth","truth","truth","truth")

# Data: Parse -------------------------------------------------------

mock.parsed <- metaphlan_parser(mock.metaphlan)

mock.species <- mock.parsed[which(mock.parsed$level == "s"),]
rownames(mock.species) <- mock.species$s
mock.species <- mock.species[,11:ncol(mock.species)]
mock.species <- mock.species[, which(colnames(mock.species) %in% mock.names)]
mock.species <- mock.species[which(rowSums(mock.species) > 0),]

#get a list of species in the mock community
keeps <- c("Pseudomonas_aeruginosa",
           "Escherichia_coli",
           "Salmonella_enterica",
           "Lactobacillus_fermentum",
           "Enterococcus_faecalis",
           "Staphylococcus_aureus",
           "Listeria_monocytogenes",
           "Bacillus_subtilis",
           "Saccharomyces_cerevisiae",
           "Cryptococcus_neoformans")

mock.species <- mock.species[which(rownames(mock.species) %in% keeps),]
mock.truth <- c(103000,
           139000,
           152000,
           146000,
           216000,
           85000,
           87000,
           61000,
           5700,
           3700)

mock.species$truth <- mock.truth

# Analysis: PCA -----------------------------------------------------------

mock_species.pca <- prcomp(t(mock.species), center = T, scale = T)
mock_species_pca.df <- mock_species.pca$x
mock_species_pca.df <- as.data.frame(mock_species_pca.df)
mock_species_pca.df$group <- mock_metadata.df[rownames(mock_species_pca.df), "Prep.Kit"]


pdf("~/Chris/dev/R_projects/ht_metagenomes/analysis/mock_metaphlan_taxa_vs_truth_prcomp.pdf",width = 10)

mock_species_pca.plot <- ggplot(mock_species_pca.df, aes(x = PC1,
                                                         y = PC2,
                                                         color = group))

mock_species_pca.plot +
  geom_point(size = 4, alpha = .8)+
  theme(text = element_text(size=18, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        strip.background = element_blank(),
        legend.key       = element_blank() ,
        legend.text      = element_text(size  = 16),
        legend.title     = element_blank()
  ) +
  scale_color_brewer(palette = "Dark2",direction = 1)
dev.off()

# Analysis: Correlation ---------------------------------------------------

library_truth_name <- mock_metadata.df$Prep.Kit
library_truth_conc <- mock_metadata.df$Input..ng.
library_truth_est  <- NULL
library_truth_pval <- NULL


for( i in 1:ncol(mock.species)){
  x <- cor.test(mock.species[,i], mock.species[,"truth"])
  library_truth_est <- c(library_truth_est,x$estimate)
  library_truth_pval <- c(library_truth_pval, x$p.value)
}

lib_truth_cor.df <- data.frame(type = library_truth_name,
                               conc = library_truth_conc,
                               estimate = library_truth_est,
                               pval = library_truth_pval)

lib_truth_cor.df <- lib_truth_cor.df[order(lib_truth_cor.df$type), ]
#prep data frame for correlation plots
mock.species_truth.melt <- as.data.frame(t(mock.species))
mock.species_truth.melt$id <- c(mock_metadata.df$Sample.ID)
mock.species_truth.melt$type <- c(mock_metadata.df$Prep.Kit)
mock.species_truth.melt$conc <- c(mock_metadata.df$Input..ng.)
mock.species_truth.melt <- melt(mock.species_truth.melt, id.vars = c("id", "type", "conc"))

#pull out truth counts and rep so we have a y axis to plot against
mock.species_truth.melt$gcn <-
  rep(mock.species_truth.melt[
    which(mock.species_truth.melt$type == "truth"),"value"],
    times = rep(26, times = 10)
    )
#remove the auto corrs
mock.species_truth.melt <-
  mock.species_truth.melt[-which(mock.species_truth.melt$type == "truth"),,drop=F]

mock.species_truth.melt$type <- factor(mock.species_truth.melt$type)

mock.species_truth.melt <- mock.species_truth.melt[
  order(mock.species_truth.melt$variable, mock.species_truth.melt$type),]

mock.species_truth.melt$id <- factor(mock.species_truth.melt$id,
                                     levels = mock.species_truth.melt$id[1:25]
                                    )

#my_est <- library_truth_est[1:(length(library_truth_est) -1)]
my_est <- lib_truth_cor.df$estimate[1:(length(lib_truth_cor.df$estimate) -1)]

mock_corr_plots.plot <- ggplot(mock.species_truth.melt,
                               aes(x = value,
                                   y = gcn,
                                   color = type,
                                   shape = conc)
)

mock_corr_plots.plot <-
  mock_corr_plots.plot +
  geom_point(size = 3) +
  scale_color_brewer(palette = "Dark2")+
  geom_abline(slope = 1,intercept = 0,alpha = .7)+
  facet_wrap(.~id, ncol =10)+
  annotate("text", label = "test", size = 4, x = 25000, y = 200000)+
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.text = element_text(colour = "black"),
    #axis.text.x = element_text(angle = 45, hjust =1 ),
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    legend.key = element_rect(fill = NA),
    legend.text = element_text(size =16),
    legend.title = element_text(size = 18,face = "bold" )
  ) +
  #scale_color_hue(h = c(0,360),l =45, c =100 )+
  ylab("Count (x 10,000)")+
  xlab("Count (x 10,000)")+
  labs(shape = "Input (ng)", color= "Method") +
  scale_y_continuous(breaks = c(0,50000,100000,150000,200000),
                     labels = c("0", "5", "10", "15", "20"))+
  scale_x_continuous(breaks = c(0,50000,100000,150000,200000),
                     labels = c("0", "5", "10", "15", "20"))


mock_corr_plots.build <- ggplot_build(mock_corr_plots.plot)
mock_corr_plots.build$data[[3]]$label <- round(my_est[1:(length(my_est))],2)
x <- ggplot_gtable(mock_corr_plots.build)

pdf("~/Chris/dev/R_projects/ht_metagenomes/analysis/mock_metaphlan_taxa_vs_truth.pdf", width = 16)
plot(x)
dev.off()



# Data parse all -------------------------------------------------

mock.parsed <- metaphlan_parser(mock.metaphlan)

mock.species <- mock.parsed[which(mock.parsed$level == "s"),]
rownames(mock.species) <- mock.species$s
mock.species <- mock.species[,11:ncol(mock.species)]
mock.species <- mock.species[, which(colnames(mock.species) %in% mock.names)]

mock.species <- mock.species[which(rowSums(mock.species) > 0),]

#get a list of species in the mock community
# keeps <- c("Pseudomonas_aeruginosa",
#            "Escherichia_coli",
#            "Salmonella_enterica",
#            "Lactobacillus_fermentum",
#            "Enterococcus_faecalis",
#            "Staphylococcus_aureus",
#            "Listeria_monocytogenes",
#            "Bacillus_subtilis",
#            "Saccharomyces_cerevisiae",
#            "Cryptococcus_neoformans")
#
# mock.species <- mock.species[which(rownames(mock.species) %in% keeps),]
mock.truthx <- c(103000,
                139000,
                152000,
                146000,
                216000,
                85000,
                0,
                87000,
                0,
                61000,
                0,
                0,
                5700,
                3700)

mock.species$truth <- mock.truthx


# plot  -------------------------------------------------------------------
library_truth_name <- mock_metadata.df$Prep.Kit
library_truth_conc <- mock_metadata.df$Input..ng.
library_truth_est  <- NULL
library_truth_pval <- NULL


for( i in 1:ncol(mock.species)){
  x <- cor.test(mock.species[,i], mock.species[,"truth"])
  library_truth_est <- c(library_truth_est,x$estimate)
  library_truth_pval <- c(library_truth_pval, x$p.value)
}

lib_truth_cor.df <- data.frame(type = library_truth_name,
                               conc = library_truth_conc,
                               estimate = library_truth_est,
                               pval = library_truth_pval)

lib_truth_cor.df <- lib_truth_cor.df[order(lib_truth_cor.df$type), ]
#prep data frame for correlation plots
mock.species_truth.melt <- as.data.frame(t(mock.species))
mock.species_truth.melt$id <- c(mock_metadata.df$Sample.ID)
mock.species_truth.melt$type <- c(mock_metadata.df$Prep.Kit)
mock.species_truth.melt$conc <- c(mock_metadata.df$Input..ng.)
mock.species_truth.melt <- melt(mock.species_truth.melt, id.vars = c("id", "type", "conc"))

#pull out truth counts and rep so we have a y axis to plot against
mock.species_truth.melt$gcn <-
  rep(mock.species_truth.melt[
    which(mock.species_truth.melt$type == "truth"),"value"],
    times = rep(26, times = 14)
  )
#remove the auto corrs
mock.species_truth.melt <-
  mock.species_truth.melt[-which(mock.species_truth.melt$type == "truth"),,drop=F]

mock.species_truth.melt$type <- factor(mock.species_truth.melt$type)

mock.species_truth.melt <- mock.species_truth.melt[
  order(mock.species_truth.melt$variable, mock.species_truth.melt$type),]

mock.species_truth.melt$id <- factor(mock.species_truth.melt$id,
                                     levels = mock.species_truth.melt$id[1:25]
)

#my_est <- library_truth_est[1:(length(library_truth_est) -1)]
my_est <- lib_truth_cor.df$estimate[1:(length(lib_truth_cor.df$estimate) -1)]

mock_corr_plots.plot <- ggplot(mock.species_truth.melt,
                               aes(x = value,
                                   y = gcn,
                                   color = type,
                                   shape = conc)
)

mock_corr_plots.plot <-
  mock_corr_plots.plot +
  geom_point(size = 3) +
  scale_color_brewer(palette = "Dark2")+
  geom_abline(slope = 1,intercept = 0,alpha = .7)+
  facet_wrap(.~id, ncol =10)+
  annotate("text", label = "test", size = 4, x = 25000, y = 200000)+
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.text = element_text(colour = "black"),
    #axis.text.x = element_text(angle = 45, hjust =1 ),
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    legend.key = element_rect(fill = NA),
    legend.text = element_text(size =16),
    legend.title = element_text(size = 18,face = "bold" )
  ) +
  #scale_color_hue(h = c(0,360),l =45, c =100 )+
  ylab("Count (x 10,000)")+
  xlab("Count (x 10,000)")+
  labs(shape = "Input (ng)", color= "Method") +
  scale_y_continuous(breaks = c(0,50000,100000,150000,200000),
                     labels = c("0", "5", "10", "15", "20"))+
  scale_x_continuous(breaks = c(0,50000,100000,150000,200000),
                     labels = c("0", "5", "10", "15", "20"))


mock_corr_plots.build <- ggplot_build(mock_corr_plots.plot)
mock_corr_plots.build$data[[3]]$label <- round(my_est[1:(length(my_est))],2)
x <- ggplot_gtable(mock_corr_plots.build)

pdf("~/Chris/dev/R_projects/ht_metagenomes/analysis/mock_metaphlan_taxa_vs_truth_all_taxa.pdf", width = 16)
plot(x)
dev.off()
