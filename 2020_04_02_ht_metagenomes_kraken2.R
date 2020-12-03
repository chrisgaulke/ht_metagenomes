#Title: 2020_04_02_ht_metagenomes_kraken2.R
#Author(s): Christopher A Gaulke
#Date: 2020-04-02 revised
#Project: HT Metagenomes

#This script will conduct analyses of kraken2 data

# Environment: Packages and options ----------------------------------------
library(ggplot2)
library(vegan)
library(reshape2)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
options("stringsAsFactors"=F)

# Environment: Function ---------------------------------------------------

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

split_parsed <- function(df, metadata, cn="Type", xtype){
  #df: parsed metaphlan data frame
  #metadata: metadata data frame that contains the variables for split
  #cn: colname of the split
  #xtype: Is the type being selected (e.g., Mock)
  #value returned is an object with the desired metadata and data

  temp.obj <- NULL

  select.names <- rownames(metadata[which(metadata[,cn] == xtype),])
  temp.metadf <- metadata[which(rownames(metadata) %in% select.names),]
  my.df <- df[,which(colnames(df) %in%
                                    c(colnames(df)[1:10], select.names))]

  my.df <- my.df[which(rowSums(my.df[,11:ncol(my.df)]) > 0),]

  temp.obj$names <- select.names
  temp.obj$metadata <- temp.metadf
  temp.obj$df <- my.df
  return(temp.obj)
}


# Data: Import data -------------------------------------------------------
combined_kraken2.df <- read.table(
  "analysis/flat_files/combined_kraken2.txt",
  sep = "\t",
  header = T,
  row.names = NULL,
  comment.char = "",
  quote = ""
  )

combined_metadata.df <- read.table(
  "analysis/flat_files/combined_metadata.txt",
  sep = "\t",
  header = T,
  comment.char = ""
)
rownames(combined_metadata.df) <- combined_metadata.df$Sample.ID
colnames(combined_metadata.df) <- c("Sample", "Kit", "Input", "Type")

# Analysis: Parse Kraken2 --------------------------------------------

combined_kraken2.parsed <- kraken2_parser(combined_kraken2.df)

# Analysis: Split by type  --------------------------------------------

mock.split <- split_parsed(df = combined_kraken2.parsed,
                           metadata = combined_metadata.df,
                           cn = "Type",
                           xtype = "Mock")

feces.split <- split_parsed(df = combined_kraken2.parsed,
                            metadata = combined_metadata.df,
                            cn = "Type",
                            xtype = "Feces")

coral.split <- split_parsed(df = combined_kraken2.parsed,
                            metadata = combined_metadata.df,
                            cn = "Type",
                            xtype = "Coral")

soil.split <- split_parsed(df = combined_kraken2.parsed,
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

pdf("analysis/figs/mock_kraken2_pca.pdf",
    height = 4)
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
        "NexteraFlex Full" = "#1B9E77",
        "NexteraFlex Reduced" = "#D95F02",
        "NextFlex RAPID XP" = "#7570B3",
        "plexWell96" = "#E7298A",
        "QIASeqFX" = "#66A61E"
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
    input = factor(mock.split$metadata$Input),
    kit = mock.split$metadata$Kit
  ),
  col = list(
    kit = c(
      "NexteraFlex Full" = "#1B9E77",
      "NexteraFlex Reduced" = "#D95F02",
      "NextFlex RAPID XP" = "#7570B3",
      "plexWell96" = "#E7298A",
      "QIASeqFX" = "#66A61E"
    ),
    input = c("0.5" = "#DEEBF7",
              "1" = "#9ECAE1",
              "5" = "#3182BD")
  ),
  show_legend = FALSE

  )

pdf("analysis/figs/mock_cor_heatmap_kraken2.pdf",
    height = 5)
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

pdf("analysis/figs/mock_intraind_div_kraken2.pdf",
    height = 4)

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

pdf("analysis/figs/mock_interind_div_kraken2.pdf",
    height = 4)

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

pdf("analysis/figs/feces_kraken2_pca.pdf",
    height = 4)
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
        "NexteraFlex Full" = "#1B9E77",
        "NexteraFlex Reduced" = "#D95F02",
        "NextFlex RAPID XP" = "#7570B3",
        "plexWell96" = "#E7298A",
        "QIASeqFX" = "#66A61E"
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
    input = factor(feces.split$metadata$Input),
    kit = feces.split$metadata$Kit
  ),
  col = list(
    kit = c(
      "NexteraFlex Full" = "#1B9E77",
      "NexteraFlex Reduced" = "#D95F02",
      "NextFlex RAPID XP" = "#7570B3",
      "plexWell96" = "#E7298A",
      "QIASeqFX" = "#66A61E"
    ),
    input = c("0.5" = "#DEEBF7",
              "1" = "#9ECAE1",
              "5" = "#3182BD")
  ),
  show_legend = FALSE

  )


pdf("analysis/figs/feces_cor_heatmap_kraken2.pdf",
    height = 5)
Heatmap(feces.split$cor,
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

pdf("analysis/figs/feces_intraind_div_kraken2.pdf",
    height = 4)

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

pdf("analysis/figs/feces_interind_div_kraken2.pdf",
    height = 4)

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

pdf("analysis/figs/soil_kraken2_pca.pdf",
    height = 4)
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
        "NexteraFlex Full" = "#1B9E77",
        "NexteraFlex Reduced" = "#D95F02",
        "NextFlex RAPID XP" = "#7570B3",
        "plexWell96" = "#E7298A",
        "QIASeqFX" = "#66A61E"
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
    input = factor(soil.split$metadata$Input),
    kit = soil.split$metadata$Kit
  ),
  col = list(
    kit = c(
      "NexteraFlex Full" = "#1B9E77",
      "NexteraFlex Reduced" = "#D95F02",
      "NextFlex RAPID XP" = "#7570B3",
      "plexWell96" = "#E7298A",
      "QIASeqFX" = "#66A61E"
    ),
    input = c("0.5" = "#DEEBF7",
              "1" = "#9ECAE1",
              "5" = "#3182BD")
  ),
  show_legend = FALSE

  )

pdf("analysis/figs/soil_cor_heatmap_kraken2.pdf",
    height = 5)
Heatmap(soil.split$cor,
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

pdf("analysis/figs/soil_intraind_div_kraken2.pdf",
    height = 4)

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

pdf("analysis/figs/soil_interind_div_kraken2.pdf",
    height = 4)

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

pdf("analysis/figs/coral_kraken2_pca.pdf",
    height = 4)
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
        "NexteraFlex Full" = "#1B9E77",
        "NexteraFlex Reduced" = "#D95F02",
        "NextFlex RAPID XP" = "#7570B3",
        "plexWell96" = "#E7298A",
        "QIASeqFX" = "#66A61E"
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
    input = factor(coral.split$metadata$Input),
    kit = coral.split$metadata$Kit
  ),
  col = list(
    kit = c(
      "NexteraFlex Full" = "#1B9E77",
      "NexteraFlex Reduced" = "#D95F02",
      "NextFlex RAPID XP" = "#7570B3",
      "plexWell96" = "#E7298A",
      "QIASeqFX" = "#66A61E"
    ),
    input = c("0.5" = "#DEEBF7",
              "1" = "#9ECAE1",
              "5" = "#3182BD")
  ),
  show_legend = FALSE

  )

pdf("analysis/figs/coral_cor_heatmap_kraken2.pdf",
    height = 5)
Heatmap(coral.split$cor,
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

pdf("analysis/figs/coral_intraind_div_kraken2.pdf",
    height = 4)

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

pdf("analysis/figs/coral_interind_div_kraken2.pdf",
    height = 4)

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

# Analysis: Adonis  -------------------------------------------------------

set.seed(731)
mock_kraken2.adonis  <- adonis(t(mock.split$species) ~  mock.split$metadata$Kit  * mock.split$metadata$Input,permutations = 5000)
feces_kraken2.adonis <- adonis(t(feces.split$species) ~ feces.split$metadata$Kit * feces.split$metadata$Input,permutations = 5000)
soil_kraken2.adonis  <- adonis(t(soil.split$species) ~  soil.split$metadata$Kit  * soil.split$metadata$Input,permutations = 5000)
coral_kraken2.adonis <- adonis(t(coral.split$species) ~ coral.split$metadata$Kit * coral.split$metadata$Input,permutations = 5000)


# Analysis: intra-individual variation ------------------------------------

mock_kraken2_intra.kruskal <- kruskal.test(mock.split$dist.intra$value, g= factor(mock.split$dist.intra$v1kit))
feces_kraken2_intra.kruskal <- kruskal.test(feces.split$dist.intra$value, g= factor(feces.split$dist.intra$v1kit))
soil_kraken2_intra.kruskal <- kruskal.test(soil.split$dist.intra$value, g= factor(soil.split$dist.intra$v1kit))
coral_kraken2_intra.kruskal <- kruskal.test(coral.split$dist.intra$value, g= factor(coral.split$dist.intra$v1kit))

pairwise.t.test(feces.split$dist.inter$value, g= factor(feces.split$dist.inter$v1kit))


# Analysis: inter-individual variation ------------------------------------

mock_kraken2_inter.kruskal <- kruskal.test(mock.split$dist.inter$value, g= factor(mock.split$dist.inter$v1kit))
feces_kraken2_inter.kruskal <- kruskal.test(feces.split$dist.inter$value, g= factor(feces.split$dist.inter$v1kit))
soil_kraken2_inter.kruskal <- kruskal.test(soil.split$dist.inter$value, g= factor(soil.split$dist.inter$v1kit))
coral_kraken2_inter.kruskal <- kruskal.test(coral.split$dist.inter$value, g= factor(coral.split$dist.inter$v1kit))

#coral had some significance
pairwise.t.test(coral.split$dist.inter$value, g= factor(coral.split$dist.inter$v1kit))
