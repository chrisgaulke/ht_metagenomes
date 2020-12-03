#Title: 2020_03_25_ht_metagenomes_genes.R
#Author(s): Christopher A Gaulke
#Date: 2020-03-25
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

# Data: Import data -------------------------------------------------------
combined_gene.df <- read.table(
  "/Users/gaulkec/Chris/dev/R_projects/ht_metagenomes/analysis/flat_files/combined_gene.txt",
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


# Analysis: Grab top level ---------------------------------------------------

#get rid of species level gene classification
combined_gene_top_level.df <- combined_gene.df[grep(pattern = "\\|", x = combined_gene.df$GeneID, invert = T),]

# Analysis: Split by type  --------------------------------------------

#Mock genes
mock_gene.names <- rownames(combined_metadata.df[which(combined_metadata.df$Type == "Mock"),])
mock_split.gene <- combined_gene_top_level.df[,which(colnames(combined_gene_top_level.df) %in% c("GeneID",mock_gene.names))]

rownames(mock_split.gene) <- mock_split.gene$GeneID
mock_split.gene$GeneID <- NULL

mock_split.gene <- mock_split.gene[which(rowSums(mock_split.gene) > 0),]

#Feces genes
feces_gene.names <- rownames(combined_metadata.df[which(combined_metadata.df$Type == "Feces"),])
feces_split.gene <- combined_gene_top_level.df[,which(colnames(combined_gene_top_level.df) %in% c("GeneID",feces_gene.names))]

rownames(feces_split.gene) <- feces_split.gene$GeneID
feces_split.gene$GeneID <- NULL

feces_split.gene <- feces_split.gene[which(rowSums(feces_split.gene) > 0),]

#Soil genes
soil_gene.names <- rownames(combined_metadata.df[which(combined_metadata.df$Type == "Soil"),])
soil_split.gene <- combined_gene_top_level.df[,which(colnames(combined_gene_top_level.df) %in% c("GeneID",soil_gene.names))]

rownames(soil_split.gene) <- soil_split.gene$GeneID
soil_split.gene$GeneID <- NULL

soil_split.gene <- soil_split.gene[which(rowSums(soil_split.gene) > 0),]

#Coral genes

coral_gene.names <- rownames(combined_metadata.df[which(combined_metadata.df$Type == "Coral"),])
coral_split.gene <- combined_gene_top_level.df[,which(colnames(combined_gene_top_level.df) %in% c("GeneID",coral_gene.names))]

rownames(coral_split.gene) <- coral_split.gene$GeneID
coral_split.gene$GeneID <- NULL

coral_split.gene <- coral_split.gene[which(rowSums(coral_split.gene) > 0),]


# Analysis: Make Objects  --------------------------------------------

mock_gene.obj <-NULL
mock_gene.obj$genes <- mock_split.gene
mock_gene.obj$metadata <- combined_metadata.df[which(combined_metadata.df$Type == "Mock"),]

feces_gene.obj <-NULL
feces_gene.obj$genes <- feces_split.gene
feces_gene.obj$metadata <- combined_metadata.df[which(combined_metadata.df$Type == "Feces"),]

soil_gene.obj <-NULL
soil_gene.obj$genes <- soil_split.gene
soil_gene.obj$metadata <- combined_metadata.df[which(combined_metadata.df$Type == "Soil"),]

coral_gene.obj <-NULL
coral_gene.obj$genes <- coral_split.gene
coral_gene.obj$metadata <- combined_metadata.df[which(combined_metadata.df$Type == "Coral"),]

# Analysis: Mock Gene PCA  --------------------------------------------

mock_gene.obj$pca <- prcomp(t(mock_gene.obj$genes),scale. = T, center = T)
mock_gene.obj$pca.scores <- data.frame(scores(mock_gene.obj$pca)[,1:5],
                                    kit = mock_gene.obj$metadata$Kit,
                                    conc = mock_gene.obj$metadata$Input)

# PCA plot
mock_gene.pca_plot <- ggplot(mock_gene.obj$pca.scores,
                        aes(x = PC1,
                            y = PC2,
                            color = kit,
                            shape = factor(conc)
                        )
)

pdf("~/Chris/dev/R_projects/ht_metagenomes/analysis/mock_gene_pca.pdf",
    height = 4)
mock_gene.pca_plot +
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

mock_gene.obj$cor <- cor(mock_gene.obj$genes)

#make the annotation layer
ta <-
  HeatmapAnnotation(
    df = data.frame(
      kit = mock_gene.obj$metadata$Kit,
      input = factor(mock_gene.obj$metadata$Input)
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
    input = factor(mock_gene.obj$metadata$Input),
    kit = mock_gene.obj$metadata$Kit
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


pdf("~/Chris/dev/R_projects/ht_metagenomes/analysis/mock_gene_cor_heatmap.pdf",
    height = 5)
Heatmap(mock_gene.obj$cor,
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

mock_gene.obj$dist <- as.matrix(
  vegdist(t(mock_gene.obj$genes), method = "bray")
)

#mark bottom tri for removal
mock_gene.obj$dist[lower.tri(mock_gene.obj$dist,diag = T)] <- NA

#melt and label
mock_gene.obj$dist.df <- melt(mock_gene.obj$dist)
mock_gene.obj$dist.df$v1kit <- rep(mock_gene.obj$metadata$Kit,
                                length(rownames(mock_gene.obj$dist)))
mock_gene.obj$dist.df$v2kit <-
  rep(mock_gene.obj$metadata$Kit,
      rep(length(rownames(mock_gene.obj$dist)),
          length(rownames(mock_gene.obj$dist))))


mock_gene.obj$dist.df$v1in <- rep(mock_gene.obj$metadata$Input,
                               length(rownames(mock_gene.obj$dist)))
mock_gene.obj$dist.df$v2in <-
  rep(mock_gene.obj$metadata$Input,
      rep(length(rownames(mock_gene.obj$dist)),
          length(rownames(mock_gene.obj$dist))))

#now remove NAs (lower tri)
mock_gene.obj$dist.df <- na.omit(mock_gene.obj$dist.df)

#so not exactly intra individual variation, but more like intra-replicate
mock_gene.obj$dist.intra <-
  mock_gene.obj$dist.df[which(
    mock_gene.obj$dist.df$v1kit == mock_gene.obj$dist.df$v2kit &
      mock_gene.obj$dist.df$v1in == mock_gene.obj$dist.df$v2in
  ), ]

#so not exactly interindividual variation, but more like inter-replicate

mock_gene.obj$dist.inter <-
  mock_gene.obj$dist.df[which(
    mock_gene.obj$dist.df$v1kit == mock_gene.obj$dist.df$v2kit &
      mock_gene.obj$dist.df$v1in != mock_gene.obj$dist.df$v2in
  ), ]

#plot technical variation
mock_gene_intra.plot <- ggplot(mock_gene.obj$dist.intra,
                          aes(x = v1kit,
                              y = value,
                              fill = v1kit))

pdf("~/Chris/dev/R_projects/ht_metagenomes/analysis/mock_gene_intraind_div.pdf",
    height = 4)

mock_gene_intra.plot +
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

mock_gene_inter.plot <- ggplot(mock_gene.obj$dist.inter,
                          aes(x = v1kit,
                              y = value,
                              fill = v1kit))

pdf("~/Chris/dev/R_projects/ht_metagenomes/analysis/mock_gene_interind_div.pdf",
    height = 4)

mock_gene_inter.plot +
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

#Analysis: Mock Alpha Diversity ---------------------------------------------------------

mock_gene.obj$richness <- specnumber(mock_gene.obj$genes,
                                     MARGIN = 2)

mock_gene.obj$richness <- as.data.frame(mock_gene.obj$richness)
mock_gene.obj$richness$kit <- mock_gene.obj$metadata[rownames(mock_gene.obj$richness),"Kit"]
mock_gene.obj$richness$input <- mock_gene.obj$metadata[rownames(mock_gene.obj$richness),"Input"]
colnames(mock_gene.obj$richness) <- c("richness", "kit", "input")
mock_gene.obj$richness.lm <- #lm(mock_gene.obj$richness$richness ~ mock_gene.obj$richness$kit + mock_gene.obj$richness$input )
  lm(mock_gene.obj$richness$richness ~ mock_gene.obj$richness$kit * mock_gene.obj$richness$input )

anova(mock_gene.obj$richness.lm)

gene_richness.plot <- ggplot(mock_gene.obj$richness,
                             aes(x = kit,
                                 y = richness,
                                 color = kit,
                                 shape = factor(input)))

pdf("~/Chris/dev/R_projects/ht_metagenomes/analysis/mock_gene_richness.pdf",
    height = 4)

gene_richness.plot +
  geom_point(size = 6, alpha = .8)+
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
        legend.title = element_text(face = "bold")
  ) +
  ylab("Richness")+
  xlab("")+
  scale_color_brewer(palette = "Dark2",direction = 1)+
  labs(shape = "Input (ng)", color= "Method")
dev.off()


mock_gene.obj$shannon <- diversity(mock_gene.obj$genes,
                                   MARGIN = 2)
mock_gene.obj$shannon <- as.data.frame(mock_gene.obj$shannon)
mock_gene.obj$shannon$kit <- mock_gene.obj$metadata[rownames(mock_gene.obj$shannon),"Kit"]
mock_gene.obj$shannon$input <- mock_gene.obj$metadata[rownames(mock_gene.obj$shannon),"Input"]
colnames(mock_gene.obj$shannon) <- c("shannon", "kit", "input")

mock_gene.obj$shannon.lm <- #lm(mock_gene.obj$shannon$shannon ~ mock_gene.obj$shannon$kit + mock_gene.obj$shannon$input )
  lm(mock_gene.obj$shannon$shannon ~ mock_gene.obj$shannon$kit * mock_gene.obj$shannon$input )

gene_shannon.plot <- ggplot(mock_gene.obj$shannon,
                            aes(x = kit,
                                y = shannon,
                                color = kit,
                                shape = factor(input)))


pdf("~/Chris/dev/R_projects/ht_metagenomes/analysis/mock_gene_shannon.pdf",
    height = 4)

gene_shannon.plot +
  geom_point(size = 6, alpha = .8)+
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
        legend.title = element_text(face = "bold")
  ) +
  ylab("Shannon Entropy")+
  xlab("")+
  scale_color_brewer(palette = "Dark2",direction = 1)+
  labs(shape = "Input (ng)", color= "Method")
dev.off()

# Analysis: feces Gene PCA  --------------------------------------------

feces_gene.obj$pca <- prcomp(t(feces_gene.obj$genes),scale. = T, center = T)
feces_gene.obj$pca.scores <- data.frame(scores(feces_gene.obj$pca)[,1:5],
                                        kit = feces_gene.obj$metadata$Kit,
                                        conc = feces_gene.obj$metadata$Input)

# PCA plot
feces_gene.pca_plot <- ggplot(feces_gene.obj$pca.scores,
                              aes(x = PC1,
                                  y = PC2,
                                  color = kit,
                                  shape = factor(conc)
                              )
)

pdf("~/Chris/dev/R_projects/ht_metagenomes/analysis/feces_gene_pca.pdf",
    height = 4)
feces_gene.pca_plot +
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

feces_gene.obj$cor <- cor(feces_gene.obj$genes)

#make the annotation layer
ta <-
  HeatmapAnnotation(
    df = data.frame(
      kit = feces_gene.obj$metadata$Kit,
      input = factor(feces_gene.obj$metadata$Input)
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
    input = factor(feces_gene.obj$metadata$Input),
    kit = feces_gene.obj$metadata$Kit

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


pdf("~/Chris/dev/R_projects/ht_metagenomes/analysis/feces_gene_cor_heatmap.pdf",
    height = 5)
Heatmap(feces_gene.obj$cor,
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

feces_gene.obj$dist <- as.matrix(
  vegdist(t(feces_gene.obj$genes), method = "bray")
)

#mark bottom tri for removal
feces_gene.obj$dist[lower.tri(feces_gene.obj$dist,diag = T)] <- NA

#melt and label
feces_gene.obj$dist.df <- melt(feces_gene.obj$dist)
feces_gene.obj$dist.df$v1kit <- rep(feces_gene.obj$metadata$Kit,
                                    length(rownames(feces_gene.obj$dist)))
feces_gene.obj$dist.df$v2kit <-
  rep(feces_gene.obj$metadata$Kit,
      rep(length(rownames(feces_gene.obj$dist)),
          length(rownames(feces_gene.obj$dist))))


feces_gene.obj$dist.df$v1in <- rep(feces_gene.obj$metadata$Input,
                                   length(rownames(feces_gene.obj$dist)))
feces_gene.obj$dist.df$v2in <-
  rep(feces_gene.obj$metadata$Input,
      rep(length(rownames(feces_gene.obj$dist)),
          length(rownames(feces_gene.obj$dist))))

#now remove NAs (lower tri)
feces_gene.obj$dist.df <- na.omit(feces_gene.obj$dist.df)

#so not exactly intra individual variation, but more like intra-replicate
feces_gene.obj$dist.intra <-
  feces_gene.obj$dist.df[which(
    feces_gene.obj$dist.df$v1kit == feces_gene.obj$dist.df$v2kit &
      feces_gene.obj$dist.df$v1in == feces_gene.obj$dist.df$v2in
  ), ]

#so not exactly interindividual variation, but more like inter-replicate

feces_gene.obj$dist.inter <-
  feces_gene.obj$dist.df[which(
    feces_gene.obj$dist.df$v1kit == feces_gene.obj$dist.df$v2kit &
      feces_gene.obj$dist.df$v1in != feces_gene.obj$dist.df$v2in
  ), ]

#plot technical variation
feces_gene_intra.plot <- ggplot(feces_gene.obj$dist.intra,
                                aes(x = v1kit,
                                    y = value,
                                    fill = v1kit))

pdf("~/Chris/dev/R_projects/ht_metagenomes/analysis/feces_gene_intraind_div.pdf",
    height = 4)

feces_gene_intra.plot +
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

feces_gene_inter.plot <- ggplot(feces_gene.obj$dist.inter,
                                aes(x = v1kit,
                                    y = value,
                                    fill = v1kit))

pdf("~/Chris/dev/R_projects/ht_metagenomes/analysis/feces_gene_interind_div.pdf",
    height = 4)

feces_gene_inter.plot +
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
#Analysis: feces Alpha Diversity ---------------------------------------------------------

feces_gene.obj$richness <- specnumber(feces_gene.obj$genes,
                                      MARGIN = 2)

feces_gene.obj$richness <- as.data.frame(feces_gene.obj$richness)
feces_gene.obj$richness$kit <- feces_gene.obj$metadata[rownames(feces_gene.obj$richness),"Kit"]
feces_gene.obj$richness$input <- feces_gene.obj$metadata[rownames(feces_gene.obj$richness),"Input"]
colnames(feces_gene.obj$richness) <- c("richness", "kit", "input")
feces_gene.obj$richness.lm <- #lm(feces_gene.obj$richness$richness ~ feces_gene.obj$richness$kit + feces_gene.obj$richness$input )
  lm(feces_gene.obj$richness$richness ~ feces_gene.obj$richness$kit * feces_gene.obj$richness$input )

anova(feces_gene.obj$richness.lm)

feces_gene_richness.plot <- ggplot(feces_gene.obj$richness,
                                   aes(x = kit,
                                       y = richness,
                                       color = kit,
                                       shape = factor(input)))

pdf("~/Chris/dev/R_projects/ht_metagenomes/analysis/feces_gene_richness.pdf",
    height = 4)

feces_gene_richness.plot +
  geom_point(size = 6, alpha = .8)+
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
        legend.title = element_text(face = "bold")
  ) +
  ylab("Richness")+
  xlab("")+
  scale_color_brewer(palette = "Dark2",direction = 1)+
  labs(shape = "Input (ng)", color= "Method")
dev.off()


feces_gene.obj$shannon <- diversity(feces_gene.obj$genes,
                                    MARGIN = 2)
feces_gene.obj$shannon <- as.data.frame(feces_gene.obj$shannon)
feces_gene.obj$shannon$kit <- feces_gene.obj$metadata[rownames(feces_gene.obj$shannon),"Kit"]
feces_gene.obj$shannon$input <- feces_gene.obj$metadata[rownames(feces_gene.obj$shannon),"Input"]
colnames(feces_gene.obj$shannon) <- c("shannon", "kit", "input")

feces_gene.obj$shannon.lm <- #lm(feces_gene.obj$shannon$shannon ~ feces_gene.obj$shannon$kit + feces_gene.obj$shannon$input )
  lm(feces_gene.obj$shannon$shannon ~ feces_gene.obj$shannon$kit * feces_gene.obj$shannon$input )
feces_gene_shannon.plot <- ggplot(feces_gene.obj$shannon,
                                  aes(x = kit,
                                      y = shannon,
                                      color = kit,
                                      shape = factor(input)))


pdf("~/Chris/dev/R_projects/ht_metagenomes/analysis/feces_gene_shannon.pdf",
    height = 4)

feces_gene_shannon.plot +
  geom_point(size = 6, alpha = .8)+
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
        legend.title = element_text(face = "bold")
  ) +
  ylab("Shannon Entropy")+
  xlab("")+
  scale_color_brewer(palette = "Dark2",direction = 1)+
  labs(shape = "Input (ng)", color= "Method")
dev.off()
# Analysis: soil Gene PCA  --------------------------------------------

soil_gene.obj$pca <- prcomp(t(soil_gene.obj$genes),scale. = T, center = T)
soil_gene.obj$pca.scores <- data.frame(scores(soil_gene.obj$pca)[,1:5],
                                       kit = soil_gene.obj$metadata$Kit,
                                       conc = soil_gene.obj$metadata$Input)

# PCA plot
soil_gene.pca_plot <- ggplot(soil_gene.obj$pca.scores,
                             aes(x = PC1,
                                 y = PC2,
                                 color = kit,
                                 shape = factor(conc)
                             )
)

pdf("~/Chris/dev/R_projects/ht_metagenomes/analysis/soil_gene_pca.pdf",
    height = 4)
soil_gene.pca_plot +
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

soil_gene.obj$cor <- cor(soil_gene.obj$genes)

#make the annotation layer
ta <-
  HeatmapAnnotation(
    df = data.frame(
      kit = soil_gene.obj$metadata$Kit,
      input = factor(soil_gene.obj$metadata$Input)
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
    input = factor(soil_gene.obj$metadata$Input),
    kit = soil_gene.obj$metadata$Kit
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


pdf("~/Chris/dev/R_projects/ht_metagenomes/analysis/soil_gene_cor_heatmap.pdf",
    height = 5)
Heatmap(soil_gene.obj$cor,
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

soil_gene.obj$dist <- as.matrix(
  vegdist(t(soil_gene.obj$genes), method = "bray")
)

#mark bottom tri for removal
soil_gene.obj$dist[lower.tri(soil_gene.obj$dist,diag = T)] <- NA

#melt and label
soil_gene.obj$dist.df <- melt(soil_gene.obj$dist)
soil_gene.obj$dist.df$v1kit <- rep(soil_gene.obj$metadata$Kit,
                                   length(rownames(soil_gene.obj$dist)))
soil_gene.obj$dist.df$v2kit <-
  rep(soil_gene.obj$metadata$Kit,
      rep(length(rownames(soil_gene.obj$dist)),
          length(rownames(soil_gene.obj$dist))))


soil_gene.obj$dist.df$v1in <- rep(soil_gene.obj$metadata$Input,
                                  length(rownames(soil_gene.obj$dist)))
soil_gene.obj$dist.df$v2in <-
  rep(soil_gene.obj$metadata$Input,
      rep(length(rownames(soil_gene.obj$dist)),
          length(rownames(soil_gene.obj$dist))))

#now remove NAs (lower tri)
soil_gene.obj$dist.df <- na.omit(soil_gene.obj$dist.df)

#so not exactly intra individual variation, but more like intra-replicate
soil_gene.obj$dist.intra <-
  soil_gene.obj$dist.df[which(
    soil_gene.obj$dist.df$v1kit == soil_gene.obj$dist.df$v2kit &
      soil_gene.obj$dist.df$v1in == soil_gene.obj$dist.df$v2in
  ), ]

#so not exactly interindividual variation, but more like inter-replicate

soil_gene.obj$dist.inter <-
  soil_gene.obj$dist.df[which(
    soil_gene.obj$dist.df$v1kit == soil_gene.obj$dist.df$v2kit &
      soil_gene.obj$dist.df$v1in != soil_gene.obj$dist.df$v2in
  ), ]

#plot technical variation
soil_gene_intra.plot <- ggplot(soil_gene.obj$dist.intra,
                               aes(x = v1kit,
                                   y = value,
                                   fill = v1kit))

pdf("~/Chris/dev/R_projects/ht_metagenomes/analysis/soil_gene_intraind_div.pdf",
    height = 4)

soil_gene_intra.plot +
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

soil_gene_inter.plot <- ggplot(soil_gene.obj$dist.inter,
                               aes(x = v1kit,
                                   y = value,
                                   fill = v1kit))

pdf("~/Chris/dev/R_projects/ht_metagenomes/analysis/soil_gene_interind_div.pdf",
    height = 4)

soil_gene_inter.plot +
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
#Analysis: soil Alpha Diversity ---------------------------------------------------------

soil_gene.obj$richness <- specnumber(soil_gene.obj$genes,
                                     MARGIN = 2)

soil_gene.obj$richness <- as.data.frame(soil_gene.obj$richness)
soil_gene.obj$richness$kit <- soil_gene.obj$metadata[rownames(soil_gene.obj$richness),"Kit"]
soil_gene.obj$richness$input <- soil_gene.obj$metadata[rownames(soil_gene.obj$richness),"Input"]
colnames(soil_gene.obj$richness) <- c("richness", "kit", "input")
soil_gene.obj$richness.lm <- #lm(soil_gene.obj$richness$richness ~ soil_gene.obj$richness$kit + soil_gene.obj$richness$input )
  lm(soil_gene.obj$richness$richness ~ soil_gene.obj$richness$kit * soil_gene.obj$richness$input )

anova(soil_gene.obj$richness.lm)

soil_gene_richness.plot <- ggplot(soil_gene.obj$richness,
                                  aes(x = kit,
                                      y = richness,
                                      color = kit,
                                      shape = factor(input)))

pdf("~/Chris/dev/R_projects/ht_metagenomes/analysis/soil_gene_richness.pdf",
    height = 4)

soil_gene_richness.plot +
  geom_point(size = 6, alpha = .8)+
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
        legend.title = element_text(face = "bold")
  ) +
  ylab("Richness")+
  xlab("")+
  scale_color_brewer(palette = "Dark2",direction = 1)+
  labs(shape = "Input (ng)", color= "Method")
dev.off()


soil_gene.obj$shannon <- diversity(soil_gene.obj$genes,
                                   MARGIN = 2)
soil_gene.obj$shannon <- as.data.frame(soil_gene.obj$shannon)
soil_gene.obj$shannon$kit <- soil_gene.obj$metadata[rownames(soil_gene.obj$shannon),"Kit"]
soil_gene.obj$shannon$input <- soil_gene.obj$metadata[rownames(soil_gene.obj$shannon),"Input"]
colnames(soil_gene.obj$shannon) <- c("shannon", "kit", "input")

soil_gene.obj$shannon.lm <- #lm(soil_gene.obj$shannon$shannon ~ soil_gene.obj$shannon$kit + soil_gene.obj$shannon$input )
  lm(soil_gene.obj$shannon$shannon ~ soil_gene.obj$shannon$kit * soil_gene.obj$shannon$input )



soil_gene_shannon.plot <- ggplot(soil_gene.obj$shannon,
                                 aes(x = kit,
                                     y = shannon,
                                     color = kit,
                                     shape = factor(input)))


pdf("~/Chris/dev/R_projects/ht_metagenomes/analysis/soil_gene_shannon.pdf",
    height = 4)

soil_gene_shannon.plot +
  geom_point(size = 6, alpha = .8)+
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
        legend.title = element_text(face = "bold")
  ) +
  ylab("Shannon Entropy")+
  xlab("")+
  scale_color_brewer(palette = "Dark2",direction = 1)+
  labs(shape = "Input (ng)", color= "Method")
dev.off()


# Analysis: coral Gene PCA  --------------------------------------------

coral_gene.obj$pca <- prcomp(t(coral_gene.obj$genes),scale. = T, center = T)
coral_gene.obj$pca.scores <- data.frame(scores(coral_gene.obj$pca)[,1:5],
                                        kit = coral_gene.obj$metadata$Kit,
                                        conc = coral_gene.obj$metadata$Input)

# PCA plot
coral_gene.pca_plot <- ggplot(coral_gene.obj$pca.scores,
                              aes(x = PC1,
                                  y = PC2,
                                  color = kit,
                                  shape = factor(conc)
                              )
)

pdf("~/Chris/dev/R_projects/ht_metagenomes/analysis/coral_gene_pca.pdf",
    height = 4)
coral_gene.pca_plot +
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

coral_gene.obj$cor <- cor(coral_gene.obj$genes)

#make the annotation layer
ta <-
  HeatmapAnnotation(
    df = data.frame(
      kit = coral_gene.obj$metadata$Kit,
      input = factor(coral_gene.obj$metadata$Input)
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
    input = factor(coral_gene.obj$metadata$Input),
    kit = coral_gene.obj$metadata$Kit

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


pdf("~/Chris/dev/R_projects/ht_metagenomes/analysis/coral_gene_cor_heatmap.pdf",
    height = 5)
Heatmap(coral_gene.obj$cor,
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

coral_gene.obj$dist <- as.matrix(
  vegdist(t(coral_gene.obj$genes), method = "bray")
)

#mark bottom tri for removal
coral_gene.obj$dist[lower.tri(coral_gene.obj$dist,diag = T)] <- NA

#melt and label
coral_gene.obj$dist.df <- melt(coral_gene.obj$dist)
coral_gene.obj$dist.df$v1kit <- rep(coral_gene.obj$metadata$Kit,
                                    length(rownames(coral_gene.obj$dist)))
coral_gene.obj$dist.df$v2kit <-
  rep(coral_gene.obj$metadata$Kit,
      rep(length(rownames(coral_gene.obj$dist)),
          length(rownames(coral_gene.obj$dist))))


coral_gene.obj$dist.df$v1in <- rep(coral_gene.obj$metadata$Input,
                                   length(rownames(coral_gene.obj$dist)))
coral_gene.obj$dist.df$v2in <-
  rep(coral_gene.obj$metadata$Input,
      rep(length(rownames(coral_gene.obj$dist)),
          length(rownames(coral_gene.obj$dist))))

#now remove NAs (lower tri)
coral_gene.obj$dist.df <- na.omit(coral_gene.obj$dist.df)

#so not exactly intra individual variation, but more like intra-replicate
coral_gene.obj$dist.intra <-
  coral_gene.obj$dist.df[which(
    coral_gene.obj$dist.df$v1kit == coral_gene.obj$dist.df$v2kit &
      coral_gene.obj$dist.df$v1in == coral_gene.obj$dist.df$v2in
  ), ]

#so not exactly interindividual variation, but more like inter-replicate

coral_gene.obj$dist.inter <-
  coral_gene.obj$dist.df[which(
    coral_gene.obj$dist.df$v1kit == coral_gene.obj$dist.df$v2kit &
      coral_gene.obj$dist.df$v1in != coral_gene.obj$dist.df$v2in
  ), ]

#plot technical variation
coral_gene_intra.plot <- ggplot(coral_gene.obj$dist.intra,
                                aes(x = v1kit,
                                    y = value,
                                    fill = v1kit))

pdf("~/Chris/dev/R_projects/ht_metagenomes/analysis/coral_gene_intraind_div.pdf",
    height = 4)

coral_gene_intra.plot +
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

coral_gene_inter.plot <- ggplot(coral_gene.obj$dist.inter,
                                aes(x = v1kit,
                                    y = value,
                                    fill = v1kit))

pdf("~/Chris/dev/R_projects/ht_metagenomes/analysis/coral_gene_interind_div.pdf",
    height = 4)

coral_gene_inter.plot +
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

#Analysis: coral Alpha Diversity ---------------------------------------------------------

coral_gene.obj$richness <- specnumber(coral_gene.obj$genes,
                                      MARGIN = 2)

coral_gene.obj$richness <- as.data.frame(coral_gene.obj$richness)
coral_gene.obj$richness$kit <- coral_gene.obj$metadata[rownames(coral_gene.obj$richness),"Kit"]
coral_gene.obj$richness$input <- coral_gene.obj$metadata[rownames(coral_gene.obj$richness),"Input"]
colnames(coral_gene.obj$richness) <- c("richness", "kit", "input")
coral_gene.obj$richness.lm <- #lm(coral_gene.obj$richness$richness ~ coral_gene.obj$richness$kit + coral_gene.obj$richness$input )
  lm(coral_gene.obj$richness$richness ~ coral_gene.obj$richness$kit * coral_gene.obj$richness$input )

anova(coral_gene.obj$richness.lm)

coral_gene_richness.plot <- ggplot(coral_gene.obj$richness,
                             aes(x = kit,
                                 y = richness,
                                 color = kit,
                                 shape = factor(input)))

pdf("~/Chris/dev/R_projects/ht_metagenomes/analysis/coral_gene_richness.pdf",
    height = 4)

coral_gene_richness.plot +
  geom_point(size = 6, alpha = .8)+
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
        legend.title = element_text(face = "bold")
  ) +
  ylab("Richness")+
  xlab("")+
  scale_color_brewer(palette = "Dark2",direction = 1)+
  labs(shape = "Input (ng)", color= "Method")
dev.off()


coral_gene.obj$shannon <- diversity(coral_gene.obj$genes,
                                    MARGIN = 2)
coral_gene.obj$shannon <- as.data.frame(coral_gene.obj$shannon)
coral_gene.obj$shannon$kit <- coral_gene.obj$metadata[rownames(coral_gene.obj$shannon),"Kit"]
coral_gene.obj$shannon$input <- coral_gene.obj$metadata[rownames(coral_gene.obj$shannon),"Input"]
colnames(coral_gene.obj$shannon) <- c("shannon", "kit", "input")

coral_gene.obj$shannon.lm <- #lm(coral_gene.obj$shannon$shannon ~ coral_gene.obj$shannon$kit + coral_gene.obj$shannon$input )
  lm(coral_gene.obj$shannon$shannon ~ coral_gene.obj$shannon$kit * coral_gene.obj$shannon$input )

anova(coral_gene.obj$shannon.lm )

coral_gene_shannon.plot <- ggplot(coral_gene.obj$shannon,
                            aes(x = kit,
                                y = shannon,
                                color = kit,
                                shape = factor(input)))


pdf("~/Chris/dev/R_projects/ht_metagenomes/analysis/coral_gene_shannon.pdf",
    height = 4)

coral_gene_shannon.plot +
  geom_point(size = 6, alpha = .8)+
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
        legend.title = element_text(face = "bold")
  ) +
  ylab("Shannon Entropy")+
  xlab("")+
  scale_color_brewer(palette = "Dark2",direction = 1)+
  labs(shape = "Input (ng)", color= "Method")
dev.off()

# Analysis: Adonis  -------------------------------------------------------

set.seed(731)
mock_gene.adonis  <- adonis(t(mock_gene.obj$genes) ~  mock_gene.obj$metadata$Kit  * mock_gene.obj$metadata$Input,permutations = 5000)
feces_gene.adonis <- adonis(t(feces_gene.obj$genes) ~ feces_gene.obj$metadata$Kit * feces_gene.obj$metadata$Input,permutations = 5000)
soil_gene.adonis  <- adonis(t(soil_gene.obj$genes) ~  soil_gene.obj$metadata$Kit  * soil_gene.obj$metadata$Input,permutations = 5000)
coral_gene.adonis <- adonis(t(coral_gene.obj$genes) ~ coral_gene.obj$metadata$Kit * coral_gene.obj$metadata$Input,permutations = 5000)


# Analysis: intra-individual variation ------------------------------------

mock_gene_intra.kruskal <- kruskal.test(mock_gene.obj$dist.intra$value, g= factor(mock_gene.obj$dist.intra$v1kit))
feces_gene_intra.kruskal <- kruskal.test(feces_gene.obj$dist.intra$value, g= factor(feces_gene.obj$dist.intra$v1kit))
soil_gene_intra.kruskal <- kruskal.test(soil_gene.obj$dist.intra$value, g= factor(soil_gene.obj$dist.intra$v1kit))
coral_gene_intra.kruskal <- kruskal.test(coral_gene.obj$dist.intra$value, g= factor(coral_gene.obj$dist.intra$v1kit))

pairwise.wilcox.test(mock_gene.obj$dist.intra$value, g= factor(mock_gene.obj$dist.intra$v1kit))
pairwise.wilcox.test(feces_gene.obj$dist.intra$value, g= factor(feces_gene.obj$dist.intra$v1kit))
pairwise.wilcox.test(soil_gene.obj$dist.intra$value, g= factor(soil_gene.obj$dist.intra$v1kit))
pairwise.wilcox.test(coral_gene.obj$dist.intra$value, g= factor(coral_gene.obj$dist.intra$v1kit))

# Analysis: inter-individual variation ------------------------------------

mock_gene_inter.kruskal <- kruskal.test(mock_gene.obj$dist.inter$value, g= factor(mock_gene.obj$dist.inter$v1kit))
feces_gene_inter.kruskal <- kruskal.test(feces_gene.obj$dist.inter$value, g= factor(feces_gene.obj$dist.inter$v1kit))
soil_gene_inter.kruskal <- kruskal.test(soil_gene.obj$dist.inter$value, g= factor(soil_gene.obj$dist.inter$v1kit))
coral_gene_inter.kruskal <- kruskal.test(coral_gene.obj$dist.inter$value, g= factor(coral_gene.obj$dist.inter$v1kit))

pairwise.wilcox.test(mock_gene.obj$dist.inter$value, g= factor(mock_gene.obj$dist.inter$v1kit))
pairwise.wilcox.test(feces_gene.obj$dist.inter$value, g= factor(feces_gene.obj$dist.inter$v1kit))
pairwise.wilcox.test(soil_gene.obj$dist.inter$value, g= factor(soil_gene.obj$dist.inter$v1kit))
pairwise.wilcox.test(coral_gene.obj$dist.inter$value, g= factor(coral_gene.obj$dist.inter$v1kit))


# Analysis: Alpha Diversity Linear Models ---------------------------------

anova(coral_gene.obj$richness.lm)
anova(coral_gene.obj$shannon.lm)

anova(feces_gene.obj$richness.lm)
anova(feces_gene.obj$shannon.lm)

anova(soil_gene.obj$richness.lm)
anova(soil_gene.obj$shannon.lm)

anova(mock_gene.obj$richness.lm)
anova(mock_gene.obj$shannon.lm)

#summary
summary(coral_gene.obj$richness.lm)
summary(coral_gene.obj$shannon.lm)

summary(feces_gene.obj$richness.lm)
summary(feces_gene.obj$shannon.lm)

summary(soil_gene.obj$richness.lm)
summary(soil_gene.obj$shannon.lm)

summary(mock_gene.obj$richness.lm)
summary(mock_gene.obj$shannon.lm)
# Analysis: Absolute Variance Mock --------------------------

#need to make a data frame that averages the 3 replicates from the 1ng
mock_gene.obj$genes.lm <- mock_gene.obj$genes


for(name in unique(mock_gene.obj$metadata$Kit)){
  select.df <- mock_gene.obj$metadata[which(mock_gene.obj$metadata$Kit == name),]
  select.names <- rownames(select.df[which(select.df$Input == 1.0),])
  select.means <- rowMeans(mock_gene.obj$genes.lm[,select.names])
  mock_gene.obj$genes.lm <- mock_gene.obj$genes.lm[,-which(colnames(mock_gene.obj$genes.lm) %in% select.names)]
  mock_gene.obj$genes.lm[,select.names[1]] <- select.means
}

#now make a new metadate table
mock_gene.obj$metadata.lm <- mock_gene.obj$metadata[which(rownames(mock_gene.obj$metadata) %in% colnames(mock_gene.obj$genes.lm)),]

lm_scatter.df <- t(mock_gene.obj$genes.lm)
lm_scatter.df <- as.data.frame(lm_scatter.df)
lm_scatter.df <- lm_scatter.df[rownames(mock_gene.obj$metadata.lm),]

lm_scatter.df$kit <- mock_gene.obj$metadata.lm[rownames(lm_scatter.df),"Kit"]
lm_scatter.df$input <- mock_gene.obj$metadata.lm[rownames(lm_scatter.df),"Input"]

lm_scatter.melt <- melt(lm_scatter.df, id.vars = c("kit", "input"))
lm_scatter.melt <- lm_scatter.melt[order(lm_scatter.melt$input),]

temp_nextflex   <- lm_scatter.melt[which(lm_scatter.melt$kit == "NexteraFlex Full"),]
lm_scatter.melt <- lm_scatter.melt[-which(lm_scatter.melt$kit == "NexteraFlex Full"),]
lm_scatter.melt$y <- rep(temp_nextflex$value, times = rep(4, times = nrow(temp_nextflex)))

cols <- c("NexteraFlex Reduced" = "#D95F02",
          "NextFlex RAPID XP" = "#7570B3",
          "plexWell96" = "#E7298A",
          "QIASeqFX" = "#66A61E")

png("~/Chris/dev/R_projects/ht_metagenomes/analysis/mock_comparison_to_netera_flex_full.png",
    res = 300, height = 5, width = 7, units = "in")
scatter_genes.plot <- ggplot(lm_scatter.melt, aes(x = log(value),
                                                  y = log(y),
                                                  shape = factor(input),
                                                  color = kit))

scatter_genes.plot +
  geom_point(alpha = .3) +
  geom_abline(slope = 1,intercept = 0,alpha = .7)+
  facet_wrap(.~kit*input, ncol =3)+
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.text = element_text(colour = "black"),
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    legend.key = element_rect(fill = NA),
    legend.text = element_text(size =16),
    legend.title = element_text(size = 18,face = "bold" )
  ) +
  scale_color_manual(values = cols)+
  ylab("Count (log)")+
  xlab("Count (log)")+
  labs(shape = "Input (ng)", color= "Method")
dev.off()



# now get the absolute value of the difference here. While I am not sure that
# it matters (give abs()) which is subtracted I will use y-x.

y05 <- rownames(mock_gene.obj$metadata.lm[which(mock_gene.obj$metadata.lm$Kit == "NexteraFlex Full" & mock_gene.obj$metadata.lm$Input == 0.5),])
y1  <- rownames(mock_gene.obj$metadata.lm[which(mock_gene.obj$metadata.lm$Kit == "NexteraFlex Full" & mock_gene.obj$metadata.lm$Input == 1.0),])
y5  <- rownames(mock_gene.obj$metadata.lm[which(mock_gene.obj$metadata.lm$Kit == "NexteraFlex Full" & mock_gene.obj$metadata.lm$Input == 5.0),])


#add small value to enable the next part
mock_gene.obj$genes.lmplus <- mock_gene.obj$genes.lm[] + 0.001

#make null df

mock_gene_abs.df <- NULL

#now populate df

for(name in c("NextFlex RAPID XP",
              "NexteraFlex Reduced",
              "QIASeqFX",
              "plexWell96")){

  x05 <- rownames(mock_gene.obj$metadata.lm[which(mock_gene.obj$metadata.lm$Kit == name & mock_gene.obj$metadata.lm$Input == 0.5),])
  x1  <- rownames(mock_gene.obj$metadata.lm[which(mock_gene.obj$metadata.lm$Kit == name & mock_gene.obj$metadata.lm$Input == 1.0),])
  x5  <- rownames(mock_gene.obj$metadata.lm[which(mock_gene.obj$metadata.lm$Kit == name & mock_gene.obj$metadata.lm$Input == 5.0),])

  name.05 <- paste0(name, "_0.5")
  name.1 <- paste0(name, "_1.0")
  name.5 <- paste0(name, "_5.0")

  #mock_gene_abs.df[[name.05]]  <- (mock_gene.obj$genes.lmplus[,y05] - mock_gene.obj$genes.lmplus[,x05]) / mock_gene.obj$genes.lmplus[,y05]
  #mock_gene_abs.df[[name.1]]   <- (mock_gene.obj$genes.lmplus[,y1]  - mock_gene.obj$genes.lmplus[,x1]) / mock_gene.obj$genes.lmplus[,y1]
  #mock_gene_abs.df[[name.5]]   <- (mock_gene.obj$genes.lmplus[,y5]  - mock_gene.obj$genes.lmplus[,x5]) / mock_gene.obj$genes.lmplus[,y5]
  mock_gene_abs.df[[name.05]]  <-  mock_gene.obj$genes.lmplus[,x05] / mock_gene.obj$genes.lmplus[,y05]
  mock_gene_abs.df[[name.1]]   <-  mock_gene.obj$genes.lmplus[,x1] / mock_gene.obj$genes.lmplus[,y1]
  mock_gene_abs.df[[name.5]]   <-  mock_gene.obj$genes.lmplus[,x5] / mock_gene.obj$genes.lmplus[,y5]
}


#now make a data frame
mock_gene_abs.df <- as.data.frame(mock_gene_abs.df)
rownames(mock_gene_abs.df) <- rownames(mock_gene.obj$genes.lm)
mock_gene_abs.df <- t(mock_gene_abs.df)
mock_gene_abs.df <- as.data.frame(mock_gene_abs.df)

#add metadata
mock_gene_abs.df$kit <- rep(c("NextFlex RAPID XP",
      "NexteraFlex Reduced",
      "QIASeqFX",
      "plexWell96"), times = rep (3,4 ))

mock_gene_abs.df$input <- rep(c("0.5", "1.0", "5.0"), times = 4)

#melt and reorder
mock_gene_abs.melt <- melt(mock_gene_abs.df)
mock_gene_abs.melt <- mock_gene_abs.melt[order(mock_gene_abs.melt$input),]

#now add in the "truth for each

mock_gene_abs.melt$y_abundance <- c(rep(mock_gene.obj$genes.lmplus[,y05], times = rep(4, length(rownames(mock_gene.obj$genes.lmplus)))),
rep(mock_gene.obj$genes.lmplus[,y1], times = rep(4, length(rownames(mock_gene.obj$genes.lmplus)))),
rep(mock_gene.obj$genes.lmplus[,y5], times = rep(4, length(rownames(mock_gene.obj$genes.lmplus))))
)

#get names of genes with that a 2 sd beyond the mean for a group
bottom_limit <- mean(mock_gene_abs.melt$value) - sd(mock_gene_abs.melt$value)*2
top_limit    <- mean(mock_gene_abs.melt$value) + sd(mock_gene_abs.melt$value)*2

select.families <- unique(mock_gene_abs.melt[which(mock_gene_abs.melt$value < bottom_limit | mock_gene_abs.melt$value > top_limit),"variable"])


mock_gene_abs.melt <- mock_gene_abs.melt[which(mock_gene_abs.melt$variable %in% select.families),]


mock_gene_abs.plot <- ggplot(mock_gene_abs.melt,
                             aes(x = log10(value),
                                 y = log10(y_abundance),
                                 color = kit,
                                 shape = factor(input)))


mock_gene_abs.plot +
  geom_point(alpha = 0.1) +
  facet_wrap(.~input)+
  geom_vline(xintercept = -(log10(abs(bottom_limit))), alpha = .7) +
  geom_vline(xintercept = log10(top_limit), alpha = .7)




mock_gene_abs.plot <- ggplot(mock_gene_abs.melt,
                             aes(x = value,
                                 y = y_abundance,
                                 color = kit,
                                 shape = factor(input)))


mock_gene_abs.plot +
  geom_point(alpha = 0.1) +
  facet_wrap(.~input)

geom_vline(xintercept = -(log10(abs(bottom_limit))), alpha = .7) +
  geom_vline(xintercept = log10(top_limit), alpha = .7)

