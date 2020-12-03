#Title: 2020_03_23_ht_metagenomes_metaphlan_analysis.R
#Author(s): Christopher A Gaulke
#Date: 2020-03-23
#Project: HT Metagenomes

#This script will conduct high level analyses of metaphlan data

# Environment: Packages and options ----------------------------------------
library(ggplot2)
library(vegan)
library(reshape2)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
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

split_parsed <- function(df, metadata, cn="Type", xtype){
  #df: parsed metaphlan data frame
  #metadata: metadata data frame that contains the variables for split
  #cn: colname of the split
  #xtype: Is the type being selected (e.g., Mock)
  #value returned is an object with the desired metadata and data

  temp.obj <- NULL

  select.names <- rownames(metadata[which(metadata[,cn] == xtype),])
  temp.metadf <- combined_metadata.df[which(rownames(combined_metadata.df) %in% select.names),]
  my.df <- combined.parsed[,which(colnames(combined.parsed) %in%
                                    c(colnames(combined.parsed)[1:10], select.names))]

  my.df <- my.df[which(rowSums(my.df[,11:ncol(my.df)]) > 0),]

  temp.obj$names <- select.names
  temp.obj$metadata <- temp.metadf
  temp.obj$df <- my.df
  return(temp.obj)
}

# Data: Import data -------------------------------------------------------
combined_metaphlan.df <- read.table(
  "/Users/gaulkec/Chris/dev/R_projects/ht_metagenomes/analysis/flat_files/combined_metaphlan.txt",
  sep = "\t",
  header = T,
  comment.char = ""
)

combined_metadata.df <- read.table(
  "/Users/gaulkec/Chris/dev/R_projects/ht_metagenomes/analysis/flat_files/combined_metadata.txt",
  sep = "\t",
  header = T,
  comment.char = ""
)
rownames(combined_metadata.df) <- combined_metadata.df$Sample.ID
colnames(combined_metadata.df) <- c("Sample", "Kit", "Input", "Type")

# Analysis: Parse Metaphlan --------------------------------------------

combined.parsed <- metaphlan_parser(combined_metaphlan.df)
combined.parsed <- combined.parsed[,c(colnames(combined.parsed)[1:10],
                                      rownames(combined_metadata.df))]
# Analysis: Split by type  --------------------------------------------

mock.split <- split_parsed(df = combined.parsed,
                           metadata = combined_metadata.df,
                           cn = "Type",
                           xtype = "Mock")

feces.split <- split_parsed(df = combined.parsed,
                           metadata = combined_metadata.df,
                           cn = "Type",
                           xtype = "Feces")

coral.split <- split_parsed(df = combined.parsed,
                           metadata = combined_metadata.df,
                           cn = "Type",
                           xtype = "Coral")

soil.split <- split_parsed(df = combined.parsed,
                           metadata = combined_metadata.df,
                           cn = "Type",
                           xtype = "Soil")

# Analysis: Mock PCA  --------------------------------------------

mock.split$species <- mock.split$df[which(mock.split$df$level == "s"),]
rownames(mock.split$species) <- mock.split$species$s
mock.split$species <- mock.split$species[,11:ncol(mock.split$species)]

#PCA
mock.split$pca <- prcomp(t(mock.split$species),scale. = T, center = T)
mock.split$pca.scores <- data.frame(scores(mock.split$pca)[,1:5],
                                    kit = mock.split$metadata$Kit,
                                    conc = mock.split$metadata$Input)

# PCA plot
mock.pca_plot <- ggplot(mock.split$pca.scores,
                        aes(x = PC1,
                            y = PC2,
                            color = kit,
                            shape = factor(conc)
                            )
                        )
ra <-
  rowAnnotation(df = data.frame(
    kit = mock.split$metadata$Kit,
    input = factor(mock.split$metadata$Input)
  ),
  col = list(
    kit = c(
      "NextFlex RAPID XP" = "#E41A1C",
      "NexteraFlex Full" = "#377EB8",
      "NexteraFlex Reduced" = "#4DAF4A",
      "QIASeqFX" = "#984EA3",
      "plexWell96" = "#FF7F00"
    ),
    input = c("0.5" = "#DEEBF7",
              "1" = "#9ECAE1",
              "5" = "#3182BD")
  ),
  show_legend = FALSE

  )

pdf("~/Chris/dev/R_projects/ht_metagenomes/analysis/mock_metaphlan_pca.pdf")
mock.pca_plot +
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
        legend.text      = element_text(size  = 16)
  ) +
  scale_color_brewer("Kit", palette = "Dark2",direction = 1)+
  scale_shape_discrete("Input (ng)")

dev.off()
# Analysis: Mock Correlation Heatmap  ------------------------------------

mock.split$cor <- cor(mock.split$species)

#make the annotation layer
ta <-
  HeatmapAnnotation(
    df = data.frame(
      kit = mock.split$metadata$Kit,
      input = factor(mock.split$metadata$Input)
    ),
    col = list(
      kit = c(
        "NextFlex RAPID XP" = "#E41A1C",
        "NexteraFlex Full" = "#377EB8",
        "NexteraFlex Reduced" = "#4DAF4A",
        "QIASeqFX" = "#984EA3",
        "plexWell96" = "#FF7F00"
      ),
      input = c("0.5" = "#DEEBF7",
                "1" = "#9ECAE1",
                "5" = "#3182BD")
    ),
    annotation_legend_param = list(kit   = list(title = "Kit",
                                              title_gp = gpar(fontsize = 14)),
                                   input = list(title = "Input",
                                                title_gp = gpar(fontsize = 14))
                                  )
  )

pdf("~/Chris/dev/R_projects/ht_metagenomes/analysis/mock_cor_heatmap.pdf")
Heatmap(mock.split$cor,
        top_annotation = ta,
        col = RColorBrewer::brewer.pal(9, "Purples"),
        name = "r",
        show_row_names = FALSE,
        show_column_names = FALSE,
        heatmap_legend_param = list(title_gp = gpar(fontsize = 14)
        )
) +
  ra
dev.off()
# Analysis: Mock Distance Plots  --------------------------------------------

mock.split$dist <- as.matrix(
  vegdist(t(mock.split$species), method = "bray")
)

#mark bottom tri for removal
mock.split$dist[lower.tri(mock.split$dist,diag = T)] <- NA

#melt and label
mock.split$dist.df <- melt(mock.split$dist)
mock.split$dist.df$v1kit <- rep(mock.split$metadata$Kit,
                                length(rownames(mock.split$dist)))
mock.split$dist.df$v2kit <-
  rep(mock.split$metadata$Kit,
      rep(length(rownames(mock.split$dist)),
          length(rownames(mock.split$dist))))


mock.split$dist.df$v1in <- rep(mock.split$metadata$Input,
                                length(rownames(mock.split$dist)))
mock.split$dist.df$v2in <-
  rep(mock.split$metadata$Input,
      rep(length(rownames(mock.split$dist)),
          length(rownames(mock.split$dist))))

#now remove NAs (lower tri)
mock.split$dist.df <- na.omit(mock.split$dist.df)

#so not exactly intra individual variation, but more like intra-replicate
mock.split$dist.intra <-
  mock.split$dist.df[which(
    mock.split$dist.df$v1kit == mock.split$dist.df$v2kit &
      mock.split$dist.df$v1in == mock.split$dist.df$v2in
  ), ]

#so not exactly interindividual variation, but more like inter-replicate

mock.split$dist.inter <-
  mock.split$dist.df[which(
    mock.split$dist.df$v1kit == mock.split$dist.df$v2kit &
      mock.split$dist.df$v1in != mock.split$dist.df$v2in
  ), ]

#plot technical variation
mock_intra.plot <- ggplot(mock.split$dist.intra,
                          aes(x = v1kit,
                              y = value,
                              fill = v1kit))

pdf("~/Chris/dev/R_projects/ht_metagenomes/analysis/mock_intraind_div.pdf")

mock_intra.plot +
  geom_boxplot() +
  theme(text = element_text(size = 16, color = "black"),
        panel.grid.major = element_line(colour = "grey96"),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        axis.line = element_line(colour = "black"),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text = element_text(colour = "black"),
        legend.key = element_rect(fill = NA),
        legend.title = element_blank()
  ) +
  ylab("Dissimilarity (bray)")+
  xlab("")+
  scale_fill_brewer(palette = "Dark2",direction = 1)
dev.off()

mock_inter.plot <- ggplot(mock.split$dist.inter,
                          aes(x = v1kit,
                              y = value,
                              fill = v1kit))

pdf("~/Chris/dev/R_projects/ht_metagenomes/analysis/mock_interind_div.pdf")

mock_inter.plot +
  geom_boxplot() +
  theme(text = element_text(size = 16, color = "black"),
    panel.grid.major = element_line(colour = "grey96"),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    axis.line = element_line(colour = "black"),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    axis.text = element_text(colour = "black"),
    legend.key = element_rect(fill = NA),
    legend.title = element_blank()
  ) +
  ylab("Dissimilarity (bray)")+
  xlab("")+
  scale_fill_brewer(palette = "Dark2",direction = 1)

dev.off()

# Analysis: feces PCA  --------------------------------------------

feces.split$species <- feces.split$df[which(feces.split$df$level == "s"),]
rownames(feces.split$species) <- feces.split$species$s
feces.split$species <- feces.split$species[,11:ncol(feces.split$species)]

#PCA
feces.split$pca <- prcomp(t(feces.split$species),scale. = T, center = T)
feces.split$pca.scores <- data.frame(scores(feces.split$pca)[,1:5],
                                     kit = feces.split$metadata$Kit,
                                     conc = feces.split$metadata$Input)

# PCA plot
feces.pca_plot <- ggplot(feces.split$pca.scores,
                         aes(x = PC1,
                             y = PC2,
                             color = kit,
                             shape = factor(conc)
                         )
)

pdf("~/Chris/dev/R_projects/ht_metagenomes/analysis/feces_metaphlan_pca.pdf")
feces.pca_plot +
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
        legend.text      = element_text(size  = 16)
  ) +
  scale_color_brewer("Kit", palette = "Dark2",direction = 1)+
  scale_shape_discrete("Input (ng)")
dev.off()

# Analysis: feces Correlation Heatmap  ------------------------------------

feces.split$cor <- cor(feces.split$species)

#make the annotation layer
ta <-
  HeatmapAnnotation(
    df = data.frame(
      kit = feces.split$metadata$Kit,
      input = factor(feces.split$metadata$Input)
    ),
    col = list(
      kit = c(
        "NextFlex RAPID XP" = "#E41A1C",
        "NexteraFlex Full" = "#377EB8",
        "NexteraFlex Reduced" = "#4DAF4A",
        "QIASeqFX" = "#984EA3",
        "plexWell96" = "#FF7F00"
      ),
      input = c("0.5" = "#DEEBF7",
                "1" = "#9ECAE1",
                "5" = "#3182BD")
    ),
    annotation_legend_param = list(kit   = list(title = "Kit",
                                                title_gp = gpar(fontsize = 14)),
                                   input = list(title = "Input",
                                                title_gp = gpar(fontsize = 14))
    )
  )

pdf("~/Chris/dev/R_projects/ht_metagenomes/analysis/feces_cor_heatmap.pdf")
Heatmap(feces.split$cor,
        top_annotation = ta,
        col = RColorBrewer::brewer.pal(9, "Purples"),
        name = "r",
        show_row_names = FALSE,
        show_column_names = FALSE,
        heatmap_legend_param = list(title_gp = gpar(fontsize = 14)
        )
)
dev.off()
# Analysis: feces Distance Plots  --------------------------------------------

feces.split$dist <- as.matrix(
  vegdist(t(feces.split$species), method = "bray")
)

#mark bottom tri for removal
feces.split$dist[lower.tri(feces.split$dist,diag = T)] <- NA

#melt and label
feces.split$dist.df <- melt(feces.split$dist)
feces.split$dist.df$v1kit <- rep(feces.split$metadata$Kit,
                                 length(rownames(feces.split$dist)))
feces.split$dist.df$v2kit <-
  rep(feces.split$metadata$Kit,
      rep(length(rownames(feces.split$dist)),
          length(rownames(feces.split$dist))))


feces.split$dist.df$v1in <- rep(feces.split$metadata$Input,
                                length(rownames(feces.split$dist)))
feces.split$dist.df$v2in <-
  rep(feces.split$metadata$Input,
      rep(length(rownames(feces.split$dist)),
          length(rownames(feces.split$dist))))

#now remove NAs (lower tri)
feces.split$dist.df <- na.omit(feces.split$dist.df)

#so not exactly intra individual variation, but more like intra-replicate
feces.split$dist.intra <-
  feces.split$dist.df[which(
    feces.split$dist.df$v1kit == feces.split$dist.df$v2kit &
      feces.split$dist.df$v1in == feces.split$dist.df$v2in
  ), ]

#so not exactly interindividual variation, but more like inter-replicate

feces.split$dist.inter <-
  feces.split$dist.df[which(
    feces.split$dist.df$v1kit == feces.split$dist.df$v2kit &
      feces.split$dist.df$v1in != feces.split$dist.df$v2in
  ), ]

#plot technical variation
feces_intra.plot <- ggplot(feces.split$dist.intra,
                           aes(x = v1kit,
                               y = value,
                               fill = v1kit))

pdf("~/Chris/dev/R_projects/ht_metagenomes/analysis/feces_intraind_div.pdf")

feces_intra.plot +
  geom_boxplot() +
  theme(text = element_text(size = 16, color = "black"),
        panel.grid.major = element_line(colour = "grey96"),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        axis.line = element_line(colour = "black"),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text = element_text(colour = "black"),
        legend.key = element_rect(fill = NA),
        legend.title = element_blank()
  ) +
  ylab("Dissimilarity (bray)")+
  xlab("")+
  scale_fill_brewer(palette = "Dark2",direction = 1)
dev.off()

feces_inter.plot <- ggplot(feces.split$dist.inter,
                           aes(x = v1kit,
                               y = value,
                               fill = v1kit))

pdf("~/Chris/dev/R_projects/ht_metagenomes/analysis/feces_interind_div.pdf")

feces_inter.plot +
  geom_boxplot() +
  theme(text = element_text(size = 16, color = "black"),
        panel.grid.major = element_line(colour = "grey96"),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        axis.line = element_line(colour = "black"),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text = element_text(colour = "black"),
        legend.key = element_rect(fill = NA),
        legend.title = element_blank()
  ) +
  ylab("Dissimilarity (bray)")+
  xlab("")+
  scale_fill_brewer(palette = "Dark2",direction = 1)

dev.off()

# Analysis: soil PCA  --------------------------------------------

soil.split$species <- soil.split$df[which(soil.split$df$level == "s"),]
rownames(soil.split$species) <- soil.split$species$s
soil.split$species <- soil.split$species[,11:ncol(soil.split$species)]

#PCA
soil.split$pca <- prcomp(t(soil.split$species),scale. = T, center = T)
soil.split$pca.scores <- data.frame(scores(soil.split$pca)[,1:5],
                                    kit = soil.split$metadata$Kit,
                                    conc = soil.split$metadata$Input)

# PCA plot
soil.pca_plot <- ggplot(soil.split$pca.scores,
                        aes(x = PC1,
                            y = PC2,
                            color = kit,
                            shape = factor(conc)
                        )
)

pdf("~/Chris/dev/R_projects/ht_metagenomes/analysis/soil_metaphlan_pca.pdf")
soil.pca_plot +
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
        legend.text      = element_text(size  = 16)
  ) +
  scale_color_brewer("Kit", palette = "Dark2",direction = 1)+
  scale_shape_discrete("Input (ng)")

dev.off()
# Analysis: soil Correlation Heatmap  ------------------------------------

soil.split$cor <- cor(soil.split$species)

#make the annotation layer
ta <-
  HeatmapAnnotation(
    df = data.frame(
      kit = soil.split$metadata$Kit,
      input = factor(soil.split$metadata$Input)
    ),
    col = list(
      kit = c(
        "NextFlex RAPID XP" = "#E41A1C",
        "NexteraFlex Full" = "#377EB8",
        "NexteraFlex Reduced" = "#4DAF4A",
        "QIASeqFX" = "#984EA3",
        "plexWell96" = "#FF7F00"
      ),
      input = c("0.5" = "#DEEBF7",
                "1" = "#9ECAE1",
                "5" = "#3182BD")
    ),
    annotation_legend_param = list(kit   = list(title = "Kit",
                                                title_gp = gpar(fontsize = 14)),
                                   input = list(title = "Input",
                                                title_gp = gpar(fontsize = 14))
    )
  )

ra <-
  rowAnnotation(df = data.frame(
    kit = soil.split$metadata$Kit,
    input = factor(soil.split$metadata$Input)
  ),
  col = list(
   kit = c(
    "NextFlex RAPID XP" = "#E41A1C",
    "NexteraFlex Full" = "#377EB8",
    "NexteraFlex Reduced" = "#4DAF4A",
    "QIASeqFX" = "#984EA3",
    "plexWell96" = "#FF7F00"
   ),
    input = c("0.5" = "#DEEBF7",
            "1" = "#9ECAE1",
            "5" = "#3182BD")
    ),
  show_legend = FALSE

)

pdf("~/Chris/dev/R_projects/ht_metagenomes/analysis/soil_cor_heatmap.pdf")
Heatmap(soil.split$cor,
        top_annotation = ta,
        col = colorRamp2(c(-1, 0, 1), c("#00441B", "white", "#40004B")),
        #col = base::rev(RColorBrewer::brewer.pal(11, "PRGn")),
        name = "r",
        show_row_names = FALSE,
        show_column_names = FALSE,
        heatmap_legend_param = list(title_gp = gpar(fontsize = 14)
        )
) +
  ra
dev.off()
# Analysis: soil Distance Plots  --------------------------------------------

soil.split$dist <- as.matrix(
  vegdist(t(soil.split$species), method = "bray")
)

#mark bottom tri for removal
soil.split$dist[lower.tri(soil.split$dist,diag = T)] <- NA

#melt and label
soil.split$dist.df <- melt(soil.split$dist)
soil.split$dist.df$v1kit <- rep(soil.split$metadata$Kit,
                                length(rownames(soil.split$dist)))
soil.split$dist.df$v2kit <-
  rep(soil.split$metadata$Kit,
      rep(length(rownames(soil.split$dist)),
          length(rownames(soil.split$dist))))


soil.split$dist.df$v1in <- rep(soil.split$metadata$Input,
                               length(rownames(soil.split$dist)))
soil.split$dist.df$v2in <-
  rep(soil.split$metadata$Input,
      rep(length(rownames(soil.split$dist)),
          length(rownames(soil.split$dist))))

#now remove NAs (lower tri)
soil.split$dist.df <- na.omit(soil.split$dist.df)

#so not exactly intra individual variation, but more like intra-replicate
soil.split$dist.intra <-
  soil.split$dist.df[which(
    soil.split$dist.df$v1kit == soil.split$dist.df$v2kit &
      soil.split$dist.df$v1in == soil.split$dist.df$v2in
  ), ]

#so not exactly interindividual variation, but more like inter-replicate

soil.split$dist.inter <-
  soil.split$dist.df[which(
    soil.split$dist.df$v1kit == soil.split$dist.df$v2kit &
      soil.split$dist.df$v1in != soil.split$dist.df$v2in
  ), ]

#plot technical variation
soil_intra.plot <- ggplot(soil.split$dist.intra,
                          aes(x = v1kit,
                              y = value,
                              fill = v1kit))

pdf("~/Chris/dev/R_projects/ht_metagenomes/analysis/soil_intraind_div.pdf")

soil_intra.plot +
  geom_boxplot() +
  theme(text = element_text(size = 16, color = "black"),
        panel.grid.major = element_line(colour = "grey96"),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        axis.line = element_line(colour = "black"),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text = element_text(colour = "black"),
        legend.key = element_rect(fill = NA),
        legend.title = element_blank()
  ) +
  ylab("Dissimilarity (bray)")+
  xlab("")+
  scale_fill_brewer(palette = "Dark2",direction = 1)
dev.off()

soil_inter.plot <- ggplot(soil.split$dist.inter,
                          aes(x = v1kit,
                              y = value,
                              fill = v1kit))

pdf("~/Chris/dev/R_projects/ht_metagenomes/analysis/soil_interind_div.pdf")

soil_inter.plot +
  geom_boxplot() +
  theme(text = element_text(size = 16, color = "black"),
        panel.grid.major = element_line(colour = "grey96"),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        axis.line = element_line(colour = "black"),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text = element_text(colour = "black"),
        legend.key = element_rect(fill = NA),
        legend.title = element_blank()
  ) +
  ylab("Dissimilarity (bray)")+
  xlab("")+
  scale_fill_brewer(palette = "Dark2",direction = 1)

dev.off()

# Analysis: coral PCA  --------------------------------------------

coral.split$species <- coral.split$df[which(coral.split$df$level == "s"),]
rownames(coral.split$species) <- coral.split$species$s
coral.split$species <- coral.split$species[,11:ncol(coral.split$species)]

#PCA
coral.split$pca <- prcomp(t(coral.split$species),scale. = T, center = T)
coral.split$pca.scores <- data.frame(scores(coral.split$pca)[,1:5],
                                     kit = coral.split$metadata$Kit,
                                     conc = coral.split$metadata$Input)

# PCA plot
coral.pca_plot <- ggplot(coral.split$pca.scores,
                         aes(x = PC1,
                             y = PC2,
                             color = kit,
                             shape = factor(conc)
                         )
)

pdf("~/Chris/dev/R_projects/ht_metagenomes/analysis/coral_metaphlan_pca.pdf")
coral.pca_plot +
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
        legend.text      = element_text(size  = 16)
  ) +
  scale_color_brewer("Kit", palette = "Dark2",direction = 1)+
  scale_shape_discrete("Input (ng)")

dev.off()
# Analysis: coral Correlation Heatmap  ------------------------------------

coral.split$cor <- cor(coral.split$species)

#make the annotation layer
ta <-
  HeatmapAnnotation(
    df = data.frame(
      kit = coral.split$metadata$Kit,
      input = factor(coral.split$metadata$Input)
    ),
    col = list(
      kit = c(
        "NextFlex RAPID XP" = "#E41A1C",
        "NexteraFlex Full" = "#377EB8",
        "NexteraFlex Reduced" = "#4DAF4A",
        "QIASeqFX" = "#984EA3",
        "plexWell96" = "#FF7F00"
      ),
      input = c("0.5" = "#DEEBF7",
                "1" = "#9ECAE1",
                "5" = "#3182BD")
    ),
    annotation_legend_param = list(kit   = list(title = "Kit",
                                                title_gp = gpar(fontsize = 14)),
                                   input = list(title = "Input",
                                                title_gp = gpar(fontsize = 14))
    )
  )

ra <-
  rowAnnotation(df = data.frame(
    kit = coral.split$metadata$Kit,
    input = factor(coral.split$metadata$Input)
  ),
  col = list(
    kit = c(
      "NextFlex RAPID XP" = "#E41A1C",
      "NexteraFlex Full" = "#377EB8",
      "NexteraFlex Reduced" = "#4DAF4A",
      "QIASeqFX" = "#984EA3",
      "plexWell96" = "#FF7F00"
    ),
    input = c("0.5" = "#DEEBF7",
              "1" = "#9ECAE1",
              "5" = "#3182BD")
  ),
  show_legend = FALSE

  )

pdf("~/Chris/dev/R_projects/ht_metagenomes/analysis/coral_cor_heatmap.pdf")
Heatmap(coral.split$cor,
        top_annotation = ta,
        col = colorRamp2(c(-1, 0, 1), c("#00441B", "white", "#40004B")),
        #col = base::rev(RColorBrewer::brewer.pal(11, "PRGn")),
        name = "r",
        show_row_names = FALSE,
        show_column_names = FALSE,
        heatmap_legend_param = list(title_gp = gpar(fontsize = 14)
        )
) +
  ra
dev.off()
# Analysis: coral Distance Plots  --------------------------------------------

coral.split$dist <- as.matrix(
  vegdist(t(coral.split$species), method = "bray")
)

#mark bottom tri for removal
coral.split$dist[lower.tri(coral.split$dist,diag = T)] <- NA

#melt and label
coral.split$dist.df <- melt(coral.split$dist)
coral.split$dist.df$v1kit <- rep(coral.split$metadata$Kit,
                                 length(rownames(coral.split$dist)))
coral.split$dist.df$v2kit <-
  rep(coral.split$metadata$Kit,
      rep(length(rownames(coral.split$dist)),
          length(rownames(coral.split$dist))))


coral.split$dist.df$v1in <- rep(coral.split$metadata$Input,
                                length(rownames(coral.split$dist)))
coral.split$dist.df$v2in <-
  rep(coral.split$metadata$Input,
      rep(length(rownames(coral.split$dist)),
          length(rownames(coral.split$dist))))

#now remove NAs (lower tri)
coral.split$dist.df <- na.omit(coral.split$dist.df)

#so not exactly intra individual variation, but more like intra-replicate
coral.split$dist.intra <-
  coral.split$dist.df[which(
    coral.split$dist.df$v1kit == coral.split$dist.df$v2kit &
      coral.split$dist.df$v1in == coral.split$dist.df$v2in
  ), ]

#so not exactly interindividual variation, but more like inter-replicate

coral.split$dist.inter <-
  coral.split$dist.df[which(
    coral.split$dist.df$v1kit == coral.split$dist.df$v2kit &
      coral.split$dist.df$v1in != coral.split$dist.df$v2in
  ), ]

#plot technical variation
coral_intra.plot <- ggplot(coral.split$dist.intra,
                           aes(x = v1kit,
                               y = value,
                               fill = v1kit))

pdf("~/Chris/dev/R_projects/ht_metagenomes/analysis/coral_intraind_div.pdf")

coral_intra.plot +
  geom_boxplot() +
  theme(text = element_text(size = 16, color = "black"),
        panel.grid.major = element_line(colour = "grey96"),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        axis.line = element_line(colour = "black"),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text = element_text(colour = "black"),
        legend.key = element_rect(fill = NA),
        legend.title = element_blank()
  ) +
  ylab("Dissimilarity (bray)")+
  xlab("")+
  scale_fill_brewer(palette = "Dark2",direction = 1)
dev.off()

coral_inter.plot <- ggplot(coral.split$dist.inter,
                           aes(x = v1kit,
                               y = value,
                               fill = v1kit))

pdf("~/Chris/dev/R_projects/ht_metagenomes/analysis/coral_interind_div.pdf")

coral_inter.plot +
  geom_boxplot() +
  theme(text = element_text(size = 16, color = "black"),
        panel.grid.major = element_line(colour = "grey96"),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        axis.line = element_line(colour = "black"),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text = element_text(colour = "black"),
        legend.key = element_rect(fill = NA),
        legend.title = element_blank()
  ) +
  ylab("Dissimilarity (bray)")+
  xlab("")+
  scale_fill_brewer(palette = "Dark2",direction = 1)

dev.off()

