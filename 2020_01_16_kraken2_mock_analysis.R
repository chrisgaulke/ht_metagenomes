#Title: 2020_01_16_kraken2_mock_analysis.R
#Author(s): Christopher A Gaulke
#Date: 2020-02-24 revised 2020-12-03
#Project: HT Metagenomes

#This script quantifies associations between the kraken2 microbial abundances
# and the mock community theoretical abundances.

# Environment: Packages and options ----------------------------------------

library(ggplot2)
library(vegan)
library(reshape2)
options("stringsAsFactors"=F)


#functions

kraken2_parser <- function(df) {

  #Description:

  # function to parse kraken report into a long form data frame

  #Usage:
  #kraken_parser(df)

  #Arguments
  #df: kraken 2 mpa report output

  #Details:
  # The input report should have the tax string in the first column
  # NOT rownames. The remaining columns should con1 counts column

  #Value:
  #a dataframe with the following columns:
  #full: full tax string
  #d: domain
  #k: kingdom
  #p: phylum
  #c: class
  #o: order
  #f: family
  #g: genus
  #s: species
  #level: lowest tax level classified (e.g., species)
  #count: the number of classified reads that were binned in this taxa label


  rows <- nrow(df)
  #ndf <- as.data.frame(matrix(NA, ncol = 10, nrow = rows))
  ndf <- as.data.frame(matrix(NA, ncol = ncol(df) + 9, nrow = rows))


  colnames(ndf) <- c("full",
                     "d",
                     "k",
                     "p",
                     "c",
                     "o",
                     "f",
                     "g",
                     "s",
                     "level",
                     colnames(df)[2:ncol(df)])


  for( i in 1:nrow(df)){
    #print(i) #for debuging
    x <- unlist(strsplit(x = df[i,1], split = "\\||__"))
    y <- x[seq(from=2, to=length(x), by = 2)] #these are tax labels
    yy <- x[seq(from=1, to=length(x), by = 2)] #these are colnames

    # kraken has some funny tax strings that include multiple "s__" columns
    # below is a hackey trick to remove the final instance of "s__" and keep
    # only species
    if(length(yy) > 1 && yy[length(yy)] == yy[length(yy) - 1 ]){
      yy <- yy[-length(yy)]
      y <- y[-length(y)]
    }
    ndf[i,yy] <- y
    ndf[i,"full"]  <- df[i,1]
    ndf[i,"level"] <- x[length(x)-1]
    ndf[i,11:ncol(ndf)]  <- df[i,2:ncol(df)] #don't want the sample id col

  }

  return(ndf)

}

#make sure all numbers are three digits
make_tres <- function(vec) {
  vec <- vec
  for( i in 1:length(vec)){
    if(nchar(vec[i]) < 4){
      vec[i] <- gsub(pattern = "s", x = vec[i], replacement = "s0")
    }
  }
  return(vec)
}


# IMPORT: Data ------------------------------------------------------------

ht_metagenomes_metadata.df <- read.table("analysis/flat_files/combined_metadata.txt",
                                         sep = "\t",
                                         header = T
                                         )

ht_meta_kraken2.df <- read.table("analysis/flat_files/combined_kraken2.txt",
                                 sep = "\t",
                                 header = T,
                                 row.names = NULL,
                                 comment.char = "",
                                 quote = "")


# DATA: Parse Data --------------------------------------------------------

ht_meta_kraken2.parsed <- kraken2_parser(ht_meta_kraken2.df)

#grab names to subset df
mock.names <-
   ht_metagenomes_metadata.df[which(ht_metagenomes_metadata.df$Sample.Type == "Mock"),
                                                   "Sample.ID"]

#add in the tax levels
mock.names <- c(mock.names, colnames(ht_meta_kraken2.parsed)[1:10])

#subset
ht_meta_kraken2_parsed.mock <-
  ht_meta_kraken2.parsed[,which(colnames(ht_meta_kraken2.parsed) %in% mock.names)]

# when combined we give an average of the plexwell samples which produces a
# decimal. To keep things consistent we will round this decimal.

for(i in 31:35){
  ht_meta_kraken2_parsed.mock[,i] <- round(ht_meta_kraken2_parsed.mock[,i], digits = 0)
}

#get rid of 0 sum taxa
ht_meta_kraken2_parsed.mock <-
  ht_meta_kraken2_parsed.mock[which(rowSums(ht_meta_kraken2_parsed.mock[11:ncol(ht_meta_kraken2_parsed.mock)]) > 0),]

#keep only species level data
ht_meta_kraken2_parsed.mock <-
  ht_meta_kraken2_parsed.mock[which(ht_meta_kraken2_parsed.mock$level == "s"),]

#rarefy (this is a little hacky but works)
set.seed(731)
ht_meta_kraken2_mock.species <-
  rrarefy(t(ht_meta_kraken2_parsed.mock[,11:ncol(ht_meta_kraken2_parsed.mock)]),
          sample = 100000)

#recombine data
ht_meta_kraken2_parsed.mock[,11:ncol(ht_meta_kraken2_parsed.mock)] <-
 t(ht_meta_kraken2_mock.species[1:nrow(ht_meta_kraken2_mock.species),])

#filter 0 sum taxa again because we rarefied
ht_meta_kraken2_parsed.mock <-
  ht_meta_kraken2_parsed.mock[which(rowSums(ht_meta_kraken2_parsed.mock[11:ncol(ht_meta_kraken2_parsed.mock)]) > 0),]


#get rid of the next step to allow for all taxa to be considered
keeps <- c("Pseudomonas aeruginosa",
  "Escherichia coli",
  "Salmonella enterica",
  "Lactobacillus fermentum",
  "Enterococcus faecalis",
  "Staphylococcus aureus",
  "Listeria monocytogenes",
  "Bacillus subtilis",
  "Saccharomyces cerevisiae",
  "Cryptococcus neoformans")

ht_meta_kraken2_parsed.final <- ht_meta_kraken2_parsed.mock[which(ht_meta_kraken2_parsed.mock$s %in% keeps),]
rownames(ht_meta_kraken2_parsed.final) <- ht_meta_kraken2_parsed.final$s
ht_meta_kraken2_parsed.final <- ht_meta_kraken2_parsed.final[keeps,]

mock.actual <- c(6100, 8500, 8700, 21600, 14600, 15200, 13900, 10300, 570, 370)

ht_meta_kraken2_parsed.est <- apply(ht_meta_kraken2_parsed.final[,11:ncol(ht_meta_kraken2_parsed.final)],
      2,
      function(x) { cor.test(x, mock.actual)$estimate})

hist(ht_meta_kraken2_parsed.est)

ht_meta_kraken2_parsed.final <- ht_meta_kraken2_parsed.final[11:ncol(ht_meta_kraken2_parsed.final)]
ht_meta_kraken2_parsed.final$tax <- rownames(ht_meta_kraken2_parsed.final)
ht_meta_kraken2_parsed.melt <- melt(ht_meta_kraken2_parsed.final)


mock_meta.df <- ht_metagenomes_metadata.df[which(ht_metagenomes_metadata.df$Sample.Type == "Mock"),]
mock_meta.df$Sample.ID <- mock.names[1:25]
rownames(mock_meta.df) <- mock_meta.df$Sample.ID

ht_meta_kraken2_parsed.melt$kit <- rep(mock_meta.df[unique(as.character(ht_meta_kraken2_parsed.melt$variable)), 2], times = rep(10,times = 25))
ht_meta_kraken2_parsed.melt$conc <- rep(mock_meta.df[unique(as.character(ht_meta_kraken2_parsed.melt$variable)), 3], times = rep(10,times = 25))
ht_meta_kraken2_parsed.melt$variable <- factor(ht_meta_kraken2_parsed.melt$variable, mock.names[1:25])

ht_meta_kraken2_parsed.melt <- ht_meta_kraken2_parsed.melt[order(ht_meta_kraken2_parsed.melt$variable),]
ht_meta_kraken2_parsed.melt$truth <- rep(mock.actual, times = 25)


#pdf("~/Desktop/test.pdf", width = 11)
mock_corr_kraken2.plot <- ggplot(ht_meta_kraken2_parsed.melt,
                               aes(x = value,
                                   y = truth,
                                   color = kit,
                                   shape = factor(conc))
)

  mock_corr_kraken2.plot +
  geom_point(size = 2) +
  geom_abline(slope = 1,intercept = 0,alpha = .7)+
  facet_wrap(.~variable, ncol =10)+
 # annotate("text", label = "test", size = 4, x = 25000, y = 200000)+
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.text = element_text(colour = "black"),
    axis.text.x = element_text(angle = 45, hjust =1 ),
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    legend.key = element_rect(fill = NA),
    legend.text = element_text(size =16),
    legend.title = element_text(size = 18,face = "bold" )
  ) +
  scale_color_hue(h = c(0,360),l =45, c =100 )+
  ylab("Count")+
  xlab("Count")+
  labs(shape = "Input (ng)", color= "Method")
#dev.off()


# ANALYSIS: Full correlation ----------------------------------------------


ht_meta_kraken2_parsed.mock <-
  ht_meta_kraken2.parsed[,which(colnames(ht_meta_kraken2.parsed) %in% mock.names)]

#round the seqwell samples because they are averages
for(i in 31:35){
  ht_meta_kraken2_parsed.mock[,i] <- round(ht_meta_kraken2_parsed.mock[,i], digits = 0)
}

ht_meta_kraken2_parsed.mock <-
  ht_meta_kraken2_parsed.mock[which(rowSums(ht_meta_kraken2_parsed.mock[11:ncol(ht_meta_kraken2_parsed.mock)]) > 0),]

ht_meta_kraken2_parsed.mock <-
  ht_meta_kraken2_parsed.mock[which(ht_meta_kraken2_parsed.mock$level == "s"),]

set.seed(731)
ht_meta_kraken2_mock.species <-
  rrarefy(t(ht_meta_kraken2_parsed.mock[,11:ncol(ht_meta_kraken2_parsed.mock)]),
          sample = 100000)


ht_meta_kraken2_parsed.mock[,11:ncol(ht_meta_kraken2_parsed.mock)] <-
  t(ht_meta_kraken2_mock.species[1:nrow(ht_meta_kraken2_mock.species),])

ht_meta_kraken2_parsed.mock <-
  ht_meta_kraken2_parsed.mock[which(rowSums(ht_meta_kraken2_parsed.mock[11:ncol(ht_meta_kraken2_parsed.mock)]) > 0),]


#get rid of the next step to allow for all taxa to be considered
keeps <- c("Pseudomonas aeruginosa",
           "Escherichia coli",
           "Salmonella enterica",
           "Lactobacillus fermentum",
           "Enterococcus faecalis",
           "Staphylococcus aureus",
           "Listeria monocytogenes",
           "Bacillus subtilis",
           "Saccharomyces cerevisiae",
           "Cryptococcus neoformans")

keeps.index <- sapply(keeps, function(x) {grep(x,ht_meta_kraken2_parsed.mock$s)})
mock.actual <- c(6100, 8500, 8700, 21600, 14600, 15200, 13900, 10300, 570, 370)

#make a vector and populate it with the actual values from the mock community.
# all other values get 0s
mock.full <- rep(0,nrow(ht_meta_kraken2_parsed.mock))

for(i in 1:length(keeps.index)){
  mock.full[keeps.index[i]] <- mock.actual[i]
}

#ht_meta_kraken2_parsed.final <- ht_meta_kraken2_parsed.mock[which(ht_meta_kraken2_parsed.mock$s %in% keeps),]
ht_meta_kraken2_parsed.final <- ht_meta_kraken2_parsed.mock

rownames(ht_meta_kraken2_parsed.final) <- ht_meta_kraken2_parsed.final$s
#ht_meta_kraken2_parsed.final <- ht_meta_kraken2_parsed.final[keeps,]


ht_meta_kraken2_parsed.est.full <- apply(ht_meta_kraken2_parsed.final[,11:ncol(ht_meta_kraken2_parsed.final)],
                                    2,
                                    function(x) { cor.test(x, mock.full)$estimate})

hist(ht_meta_kraken2_parsed.est)
hist(ht_meta_kraken2_parsed.est.full)

ht_meta_kraken2_parsed.final <- ht_meta_kraken2_parsed.final[11:ncol(ht_meta_kraken2_parsed.final)]
ht_meta_kraken2_parsed.final$tax <- rownames(ht_meta_kraken2_parsed.final)
ht_meta_kraken2_parsed.melt <- melt(ht_meta_kraken2_parsed.final)


mock_meta.df <- ht_metagenomes_metadata.df[which(ht_metagenomes_metadata.df$Sample.Type == "Mock"),]
mock_meta.df$Sample.ID <- mock.names[1:25]
rownames(mock_meta.df) <- mock_meta.df$Sample.ID

ht_meta_kraken2_parsed.melt$kit <- rep(mock_meta.df[unique(as.character(ht_meta_kraken2_parsed.melt$variable)), 2], times = rep(nrow(ht_meta_kraken2_parsed.mock),times = 25))
ht_meta_kraken2_parsed.melt$conc <- rep(mock_meta.df[unique(as.character(ht_meta_kraken2_parsed.melt$variable)), 3], times = rep(nrow(ht_meta_kraken2_parsed.mock),times = 25))
ht_meta_kraken2_parsed.melt$variable <- factor(ht_meta_kraken2_parsed.melt$variable, mock.names[1:25])

ht_meta_kraken2_parsed.melt <- ht_meta_kraken2_parsed.melt[order(ht_meta_kraken2_parsed.melt$variable),]
ht_meta_kraken2_parsed.melt$truth <- rep(mock.full, times = 25)

ht_meta_kraken2_parsed.melt$kit <- factor(ht_meta_kraken2_parsed.melt$kit)

ht_meta_kraken2_parsed.melt <- ht_meta_kraken2_parsed.melt[
  order(ht_meta_kraken2_parsed.melt$tax, ht_meta_kraken2_parsed.melt$kit),]

ht_meta_kraken2_parsed.melt$variable <- factor(ht_meta_kraken2_parsed.melt$variable,
                                               levels = ht_meta_kraken2_parsed.melt$variable[1:25]
                                               )


mock_corr_kraken2.plot <- ggplot(ht_meta_kraken2_parsed.melt,
                                 aes(x = value,
                                     y = truth,
                                     color = kit,
                                     shape = factor(conc))
)

mock_corr_kraken2.plot <-
  mock_corr_kraken2.plot +
  geom_point(size = 3) +
  scale_color_brewer(palette = "Dark2")+
  geom_abline(slope = 1,intercept = 0,alpha = .7)+
  facet_wrap(.~variable, ncol =10)+
  annotate("text", label = "test", size = 4, x = 5000, y = 25000)+
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.text = element_text(colour = "black"),
    axis.text.x = element_text(angle = 45, hjust =1 ),
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    legend.key = element_rect(fill = NA),
    legend.text = element_text(size =16),
    legend.title = element_text(size = 18,face = "bold" )
  ) +
  ylab("Count")+
  xlab("Count")+
  labs(shape = "Input (ng)", color= "Method")


mock_corr_kraken2_plots.build <- ggplot_build(mock_corr_kraken2.plot)
mock_corr_kraken2_plots.build$data[[3]]$label <- round(ht_meta_kraken2_parsed.est.full[1:(length(ht_meta_kraken2_parsed.est.full))],2)
x <- ggplot_gtable(mock_corr_kraken2_plots.build)

pdf("analysis/figs/mock_kraken_taxa_vs_truth_all_taxa.pdf", width = 16)
plot(x)
dev.off()


