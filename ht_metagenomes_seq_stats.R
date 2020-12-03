#Title: ht_metagenomes_seq_stats.R
#Author(s): Christopher A Gaulke
#Date: 2020-02-24 Revised 2020-12-03
#Project: HT Metagenomes

#Purpose: To summarize the library sequencing statistics and
# Environment: Packages and options ----------------------------------------

options("stringsAsFactors"=F)
library(ggplot2)
library(dplyr)
# Data: Import Data  -----------------------------------------------------------------

# Import multiqc output stats for post-shotcleaner R1 reads
# Similar results were visually confirmed with R2 reads

htmeta_seq_stats.df <- read.table("analysis/flat_files/cleaned_seq_stats.txt",
                                  sep = "\t",
                                  header = T,
                                  row.names = 1,
                                  stringsAsFactors = F)

#Data from Mark

lib_qc.df <- read.table("data/cgrb_libs_qc_stats_modified.txt",
                        sep = "\t",
                        header=T,
                        stringsAsFactors = F
                        )


# Analysis: Duplication Rates ----------------------------------------------

htmeta_dups.plot <- ggplot(htmeta_seq_stats.df,
                           aes(x = prep,
                               y = Dups,
                               fill = prep
                               )
                           )


pdf("analysis/figs/ht_meta_duplicates_plot.pdf", width = 8)

htmeta_dups.plot +
  geom_boxplot() +
  geom_point()+
  facet_wrap(.~type) +
  theme(text = element_text(size=14, colour = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line      = element_line(colour = "black"),
      legend.key       = element_blank() ,
      legend.position  = "top",
      legend.text      = element_text(size  = 16),
      legend.title     = element_blank(),
      axis.text        = element_text(color = "black"),
      axis.text.x      = element_blank(),
      axis.title.x     = element_blank(),
      axis.ticks.x     = element_blank(),
      strip.background = element_blank(),
      strip.text.x     = element_text(size = 18, face = "bold")
  ) +
  scale_y_continuous(limits = c(0,40))+
  ylab("% Duplicates")+
  xlab("") +
  guides(fill=guide_legend(nrow=2))+
  scale_fill_brewer(palette = "Dark2")

dev.off()

# Analysis: %GC ----------------------------------------------------------

htmeta_gc.plot <- ggplot(htmeta_seq_stats.df,
                         aes(x = prep,
                             y = GC,
                             fill = prep))

pdf("analysis/figs/ht_meta_gc_plot.pdf", width = 8)

htmeta_gc.plot +
  geom_boxplot() +
  geom_jitter() +
  facet_wrap(.~type) +

  theme(text = element_text(size=14, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line      = element_line(colour = "black"),
        legend.key       = element_blank() ,
        legend.position  = "top",
        legend.text      = element_text(size  = 16),
        legend.title     = element_blank(),
        axis.text        = element_text(color = "black"),
        axis.text.x      = element_blank(),
        axis.title.x     = element_blank(),
        axis.ticks.x     = element_blank(),
        strip.background = element_blank(),
        strip.text.x     = element_text(size = 18, face = "bold")
  ) +
  scale_y_continuous(limits = c(35,65))+
  ylab("% GC")+
  xlab("") +
  guides(fill=guide_legend(nrow=2))+
  scale_fill_brewer(palette = "Dark2",direction = 1)

dev.off()

# Analysis: Millions of Sequences ------------------------------------------

htmeta_nseqs.plot <- ggplot(htmeta_seq_stats.df,
                            aes(x = prep,
                                y = M.Seqs,
                                fill = prep))


pdf("analysis/figs/ht_meta_nseqs_plot.pdf",width = 8)

htmeta_nseqs.plot +
  geom_boxplot() +
  geom_jitter() +
  facet_wrap(.~type) +

  theme(text = element_text(size=14, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line      = element_line(colour = "black"),
        legend.key       = element_blank() ,
        legend.position  = "top",
        legend.text      = element_text(size  = 16),
        legend.title     = element_blank(),
        axis.text        = element_text(color = "black"),
        axis.text.x      = element_blank(),
        axis.title.x     = element_blank(),
        axis.ticks.x     = element_blank(),
        strip.background = element_blank(),
        strip.text.x     = element_text(size = 18, face = "bold")
  ) +
  scale_y_continuous(limits = c(0,8))+
  ylab("Reads (Millions)")+
  xlab("") +
  guides(fill=guide_legend(nrow=2))+
  scale_fill_brewer(palette = "Dark2",direction = 1)

dev.off()

# Analysis: Insert Size ---------------------------------------------------

insert_size.plot <- ggplot(lib_qc.df, aes(x = lib,
                                          y = Median.Fragment.Size..bp.,
                                          fill = lib))

pdf("analysis/figs/insert_size_by_kit.pdf",
    width = 8)

insert_size.plot +
  geom_boxplot() +
  ylab("Insert Size (bp)")+
  xlab("")+
  facet_wrap( . ~ type, ncol = 2)+
  theme(text = element_text(size = 16),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.text = element_text(colour = "black"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    legend.key = element_rect(fill = NA),
    legend.text = element_text(size =16),
    strip.background =element_blank(),
    strip.text = element_text(colour = 'black', size = 16,face = "bold"),
    legend.title = element_text(size = 16,face = "bold" ),
    legend.position = "top"
  )+
  guides(fill=guide_legend(nrow=2))+
  labs(fill = "")+
  scale_fill_brewer(palette = "Dark2")
dev.off()

# Analysis: Library Concentration --------------------------------------------

conc.plot <- ggplot(lib_qc.df, aes(x = lib,
                                   y = Conc...ng.uL.,
                                   fill = lib))


pdf("analysis/figs/lib_conc_by_kit.pdf",width = 8)

conc.plot +
  geom_boxplot() +
  ylab("Library Concentration")+
  xlab("")+
  facet_wrap( . ~ type, ncol = 2)+
  theme(text = element_text(size = 16),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.text = element_text(colour = "black"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    legend.key = element_rect(fill = NA),
    legend.text = element_text(size =16),
    strip.background =element_blank(),
    strip.text = element_text(colour = 'black', size = 16,face = "bold"),
    legend.title = element_text(size = 16,face = "bold" ),
    legend.position = "top"
  )+
  guides(fill=guide_legend(nrow=2))+
  labs(fill = "")+
  scale_fill_brewer(palette = "Dark2")
dev.off()

# Analysis: Number of Filtered reads --------------------------------------

pdf("analysis/figs/reads_filtered_by_kit.pdf",width = 8)

filtered.plot <- ggplot(htmeta_seq_stats.df, aes(x = prep,
                                                     y = per_filtered,
                                                     fill = prep)
)

filtered.plot +
  geom_boxplot() +
  facet_wrap( . ~ type, ncol = 2)+
  ylab("Percent of Reads Filtered")+
  xlab("")+
  facet_wrap( . ~ type, ncol = 2)+
  theme(text = element_text(size = 16),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(colour = "black"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.key = element_rect(fill = NA),
        legend.text = element_text(size =16),
        strip.background =element_blank(),
        strip.text = element_text(colour = 'black', size = 16,face = "bold"),
        legend.title = element_text(size = 16,face = "bold" ),
        legend.position = "top"
  )+
  guides(fill=guide_legend(nrow=2))+
  labs(fill = "")+
  scale_fill_brewer(palette = "Dark2")
dev.off()

# Analysis: Linear Models -------------------------------------------------

# A function to build four models and use AIC to determine if interaction terms
# are actually appropriate (i.e., reduce information loss). So with each model
# we ask: does kit choice, library type, and dna input concentration matter to
# library characteristics.

# When a models is selected then we need to rebuild the model that was the best


htmeta_seq_stats.df$type <-
  factor(htmeta_seq_stats.df$type,
         levels = c("Feces", "Coral", "Mock", "Soil"))

lib_qc.df$type <-
  factor(lib_qc.df$type,
         levels = c("Feces", "Coral", "Mock", "Soil"))


model_AIC <- function(df, i, x,y,z){
  fit0 <- lm(df[,i] ~ df[,x] + df[,y] + df[,z])
  fit3 <- lm(df[,i] ~ df[,x] * df[,y] * df[,z])
  fit.aov <- AIC(fit0, fit3)
  fit.aov
}


#inclusion of a double interaction term improved model
model_AIC(lib_qc.df, "Conc...ng.uL.", "input", "lib", "type")
model_AIC(lib_qc.df, "Median.Fragment.Size..bp.", "input", "lib", "type")
model_AIC(lib_qc.df, "Molarity..nM.", "input", "lib", "type")
model_AIC(htmeta_seq_stats.df, "Length", "type", "prep", "conc")
model_AIC(htmeta_seq_stats.df, "GC", "type", "prep", "conc")
model_AIC(htmeta_seq_stats.df, "per_retained", "type", "prep", "conc")
model_AIC(htmeta_seq_stats.df, "per_filtered", "type", "prep", "conc")

#inclusion of interactions did not improve model for Dups
model_AIC(htmeta_seq_stats.df, "Dups", "type", "prep", "conc")


# now build the most parsimonious model for each library characteristic

#make an object to put the stats data in
htmeta_seq_stats.stats <- NULL

htmeta_seq_stats.stats$duplicates <-
  lm(htmeta_seq_stats.df$Dups ~ htmeta_seq_stats.df$type +
       htmeta_seq_stats.df$prep +
       htmeta_seq_stats.df$conc)

htmeta_seq_stats.stats$length <-
  lm(htmeta_seq_stats.df$Length ~ htmeta_seq_stats.df$type *
       htmeta_seq_stats.df$prep *
       htmeta_seq_stats.df$conc)

htmeta_seq_stats.stats$m_seqs <-
  lm(htmeta_seq_stats.df$M.Seqs ~ htmeta_seq_stats.df$type *
       htmeta_seq_stats.df$prep *
       htmeta_seq_stats.df$conc)

htmeta_seq_stats.stats$gc <-
  lm(htmeta_seq_stats.df$GC ~ htmeta_seq_stats.df$type *
       htmeta_seq_stats.df$prep *
       htmeta_seq_stats.df$conc)

htmeta_seq_stats.stats$filtered <-
  lm(htmeta_seq_stats.df$per_filtered ~ htmeta_seq_stats.df$type *
       htmeta_seq_stats.df$prep *
       htmeta_seq_stats.df$conc)

htmeta_seq_stats.stats$lib_conc <-
  lm(lib_qc.df$Conc...ng.uL. ~ lib_qc.df$type *
       lib_qc.df$lib *
       lib_qc.df$input)

htmeta_seq_stats.stats$frag_size <-
  lm(lib_qc.df$Median.Fragment.Size..bp. ~ lib_qc.df$type *
       lib_qc.df$lib *
       lib_qc.df$input)

htmeta_seq_stats.stats$molarity <-
  lm(lib_qc.df$Molarity..nM. ~ lib_qc.df$type *
       lib_qc.df$lib *
       lib_qc.df$input)

#now check up on the stats
anova(htmeta_seq_stats.stats$duplicates) # only one that has the additive model selected

anova(htmeta_seq_stats.stats$length)    ##
anova(htmeta_seq_stats.stats$gc)        ##
anova(htmeta_seq_stats.stats$m_seqs)    ##

anova(htmeta_seq_stats.stats$filtered)  ## all interactions significant

anova(htmeta_seq_stats.stats$lib_conc)  ##
anova(htmeta_seq_stats.stats$frag_size) ## all interactions significant
anova(htmeta_seq_stats.stats$molarity)  ## all interactions significant


# Data: Make combined sequence stats df -----------------------------------

#make a dataframe with sumary stats
htmeta_seq_stats.min  <- NULL
htmeta_seq_stats.max  <- NULL
htmeta_seq_stats.mean <- NULL


for(name in c("Coral","Feces", "Mock","Soil")){
  for(parameter in c("Dups","Length","M.Seqs",
                     "GC", "per_filtered")){

    #mean
    htmeta_seq_stats.mean <- c(htmeta_seq_stats.mean,
                               mean(subset.data.frame(htmeta_seq_stats.df,
                                                      type == name,
                                                      select = parameter,
                                                      drop = T)
                               )
    )
    #min
    htmeta_seq_stats.min <- c(htmeta_seq_stats.min,
                              min(subset.data.frame(htmeta_seq_stats.df,
                                                    type == name,
                                                    select = parameter,
                                                    drop = T)
                              )
    )
    #max
    htmeta_seq_stats.max <- c(htmeta_seq_stats.max,
                              max(subset.data.frame(htmeta_seq_stats.df,
                                                    type == name,
                                                    select = parameter,
                                                    drop = T)
                              )
    )

  }
}

htmeta_seq_stats.min  <- round(htmeta_seq_stats.min, digits = 2)
htmeta_seq_stats.max  <- round(htmeta_seq_stats.max, digits = 2)
htmeta_seq_stats.mean <- round(htmeta_seq_stats.mean, digits = 2)

htmeta_seq_stats.paste <- paste0(htmeta_seq_stats.mean,
                                 " (",
                                 htmeta_seq_stats.min,
                                 "-",
                                 htmeta_seq_stats.max,
                                 ")"
)

htmeta_seq_stats_paste.df <- matrix(htmeta_seq_stats.paste,
                                    ncol = 5,
                                    nrow = 4,
                                    byrow = T
)

rownames(htmeta_seq_stats_paste.df) <- c("Coral",
                                         "Soil",
                                         "Feces",
                                         "Mock"
)

colnames(htmeta_seq_stats_paste.df) <- c("Duplicates",
                                         "Read Length",
                                         "Million of Reads",
                                         "%GC",
                                         "%Reads Filtered"
)

#Now the next set
htmeta_seq_stats.min  <- NULL
htmeta_seq_stats.max  <- NULL
htmeta_seq_stats.mean <- NULL


for(name in  c("Coral","Feces", "Mock","Soil")){
  for(parameter in c("Conc...ng.uL.",
                     "Median.Fragment.Size..bp.",
                     "Molarity..nM."
                     )
      )
  {
    #mean
    htmeta_seq_stats.mean <- c(htmeta_seq_stats.mean,
            mean(subset.data.frame(lib_qc.df,
                         type == name,
                         select = parameter,
                         drop = T)
       )
    )
    #min
    htmeta_seq_stats.min <- c(htmeta_seq_stats.min,
            min(subset.data.frame(lib_qc.df,
                          type == name,
                          select = parameter,
                          drop = T)
        )
    )
    #max
    htmeta_seq_stats.max <- c(htmeta_seq_stats.max,
            max(subset.data.frame(lib_qc.df,
                          type == name,
                          select = parameter,
                          drop = T)
        )
    )
  }
}

htmeta_seq_stats.min  <- round(htmeta_seq_stats.min, digits = 2)
htmeta_seq_stats.max  <- round(htmeta_seq_stats.max, digits = 2)
htmeta_seq_stats.mean <- round(htmeta_seq_stats.mean, digits = 2)


htmeta_seq_stats.paste2 <- paste0(htmeta_seq_stats.mean,
                                " (",
                                htmeta_seq_stats.min,
                                "-",
                                htmeta_seq_stats.max,
                                ")"
                                )

htmeta_seq_stats_paste.df2 <- matrix(htmeta_seq_stats.paste2,
                           ncol = 3,
                           nrow = 4,
                           byrow = T
                           )

rownames(htmeta_seq_stats_paste.df2) <- c("Coral",
                                         "Feces",
                                         "Mock",
                                         "Soil"
                                         )

colnames(htmeta_seq_stats_paste.df2) <- c("Library Concentration (ng/µL)",
                                         "Median Fragment Size",
                                         "Library Molarity (nM)"
                                         )


htmeta_seq_stats_paste.merge <- merge.data.frame(htmeta_seq_stats_paste.df,
                                                 htmeta_seq_stats_paste.df2,
                                                 by = "row.names" )


htmeta_seq_stats_paste.merge <- as.data.frame(htmeta_seq_stats_paste.merge)

row.names(htmeta_seq_stats_paste.merge) <-
  htmeta_seq_stats_paste.merge$Row.names

htmeta_seq_stats_paste.merge$Row.names <- NULL

htmeta_seq_stats_paste.merge <- htmeta_seq_stats_paste.merge[c("Coral",
                                                               "Soil",
                                                               "Feces",
                                                               "Mock"),
                                                             ]



htmeta_seq_stats_paste.merge <-
  htmeta_seq_stats_paste.merge[, c(
    "Million of Reads",
    "Median Fragment Size",
    "Library Concentration (ng/µL)",
    "Library Molarity (nM)",
    "Read Length",
    "Duplicates",
    "%GC",
    "%Reads Filtered"
  )]


write.table(x = htmeta_seq_stats_paste.merge,
            file = "analysis/flat_files/merge_seq_stats_summary.txt",
            sep = "\t",
            quote = F,
            row.names = T,
            col.names = T
            )


# Data: Make combined sequence stats df -----------------------------------

stats.df <- NULL
# for the first half in lib_qc.df
for(parameter in c("Median.Fragment.Size..bp.",
                     "Conc...ng.uL.", "Molarity..nM.")){
  x.df <- lib_qc.df %>%
            group_by(type, lib) %>%
              summarise(mean = mean(.data[[parameter]]),
            max = max(.data[[parameter]]),
            min = min(.data[[parameter]]),
            stdev = sd(.data[[parameter]]))

  stats.df[[parameter]] <- paste0(
    round(x.df$mean, digits = 2),
       " (",
    round(x.df$min, digits = 2),
       "-",
    round(x.df$max, digits = 2),
       ")"
  )
}


#for the rest in htmeta_seq_stats.df

for(parameter in c("M.Seqs","per_filtered","Length",
                   "Dups", "GC")){
  x.df <- htmeta_seq_stats.df %>%
    group_by(type, prep) %>%
    summarise(mean = mean(.data[[parameter]]),
              max = max(.data[[parameter]]),
              min = min(.data[[parameter]]),
              stdev = sd(.data[[parameter]]))

  stats.df[[parameter]] <- paste0(
    round(x.df$mean, digits = 2),
    " (",
    round(x.df$min, digits = 2),
    "-",
    round(x.df$max, digits = 2),
    ")"
  )
}

stats.df <- cbind(type = x.df$type,
                  prep = x.df$prep,
                  as.data.frame(stats.df)
                  )

stats.df <- as.data.frame(stats.df)


write.table(x = stats.df,
            file = "analysis/flat_files/merge_seq_stats_summary_by_kit.txt",
            sep = "\t",
            quote = F,
            row.names = T,
            col.names = T
)



