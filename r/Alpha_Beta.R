library(phyloseq)
library(ape)
library(tidyverse)
library(picante)
library(ggsignif)
library(vegan)
library(ggforce)


load("data/depression_phyloseq.Rdata")
load("data/depression_phyloseq_rare.Rdata")

#samp_dat_wdiv$BDI_category_antidepressant_use <- factor(samp_dat_wdiv$BDI_category_antidepressant_use, levels = c("low+no","low+yes","high+no","high+yes"))


#### alpha diversity ####
### richness measures ###
gg_richnessA <- plot_richness(depression_rare, x = "Antidepressant_use", measure = c("shannon", "Observed"))+
    xlab("Antidepressant Use")+
    geom_boxplot()

gg_richnessB <- plot_richness(depression_rare, x = "BDI_category", measure = c("shannon", "Observed"))+
  xlab("BDI Category")+
  geom_boxplot()

gg_richnessC <- plot_richness(depression_rare, x = "BDI_category_antidepressant_use", measure = "Observed") +
  xlab("BDI Category and Antidepressant Use")+
  geom_boxplot()+
  theme_bw()+
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),   
    plot.title = element_blank()  
    )

ggsave(filename = "figures/richness_plot_C.png", 
       gg_richnessC,
       height = 4, width = 6)

### faith's PD ##

phylo_dist <- pd(t(otu_table(depression_rare)), phy_tree(depression_rare),
                 include.root = F)

sample_data(depression_rare)$PD <- phylo_dist$PD

plot.pdA <- ggplot(sample_data(depression_rare), aes(Antidepressant_use, PD)) +
  geom_boxplot()+
  xlab("Antidepressant Use")+
  ylab("Phylogenetic Diversity")

plot.pdB <- ggplot(sample_data(depression_rare), aes(BDI_category, PD)) +
  geom_boxplot()+
  xlab("BDI Category")+
  ylab("Phylogenetic Diversity")

plot.pdC <- ggplot(sample_data(depression_rare), aes(BDI_category_antidepressant_use, PD)) +
  geom_boxplot()+
  xlab("BDI Category and Antidepressant Use")+
  ylab("Phylogenetic Diversity")


ggsave(filename = "figures/PD_plot_C.png", 
       plot.pdC,
       height = 4, width = 6)

### stats of observered and shannon ###

alphadiv <- estimate_richness(depression_rare)
samp_dat <- sample_data(depression_rare)
samp_dat_wdiv <- data.frame(samp_dat, alphadiv)

wilcox.test(Observed ~ Antidepressant_use, data=samp_dat_wdiv, exact = FALSE)
wilcox.test(Shannon ~ Antidepressant_use, data=samp_dat_wdiv, exact = FALSE)

wilcox.test(Observed ~ BDI_category, data=samp_dat_wdiv, exact = FALSE)
wilcox.test(Shannon ~ BDI_category, data=samp_dat_wdiv, exact = FALSE)

kruskal.test(Observed ~ BDI_category_antidepressant_use, data=samp_dat_wdiv)
kruskal.test(Shannon ~ BDI_category_antidepressant_use, data=samp_dat_wdiv)

### Full stats list of observed, shannon, FPD using anova ###

lm_ob_vs_site_log <- lm(log(Observed) ~ `BDI_category_antidepressant_use`, data=samp_dat_wdiv)
anova_ob_vs_site_log <- aov(lm_ob_vs_site_log)
summary(anova_ob_vs_site_log)
TukeyHSD(anova_ob_vs_site_log)

lm_sh_vs_site_log <- lm(log(Shannon) ~ `BDI_category_antidepressant_use`, data=samp_dat_wdiv)
anova_sh_vs_site_log <- aov(lm_sh_vs_site_log)
summary(anova_sh_vs_site_log)
TukeyHSD(anova_sh_vs_site_log)

lm_pd_vs_site_log <- lm(log(PD) ~ `BDI_category_antidepressant_use`, data=samp_dat_wdiv)
anova_pd_vs_site_log <- aov(lm_pd_vs_site_log)
summary(anova_pd_vs_site_log)
TukeyHSD(anova_pd_vs_site_log)

#### Beta Diversity ####
#### stats ####
dm_wunifrac <- UniFrac(depression_rare, weighted=TRUE) # Weighted UniFrac
dm_unifrac <- UniFrac(depression_rare, weighted=FALSE) # unweighted UniFrac
dm_bray <- vegdist(t(otu_table(depression_rare)), method="bray") # Bray-curtis
dm_jaccard <- vegdist(t(otu_table(depression_rare)), method="jaccard") # Jaccard

adonis2(dm_wunifrac ~ BDI_category*Antidepressant_use, data=samp_dat_wdiv)
adonis2(dm_unifrac ~ BDI_category*Antidepressant_use, data=samp_dat_wdiv)
adonis2(dm_bray ~ BDI_category*Antidepressant_use, data=samp_dat_wdiv)
adonis2(dm_jaccard ~ BDI_category*Antidepressant_use, data=samp_dat_wdiv)

bray <- distance(depression_rare, method = "bray")
jaccard <- distance(depression_rare, method = "jaccard")
unifrac <- distance(depression_rare, method = "unifrac")
wunifrac <- distance(depression_rare, method = "wunifrac")

pcoa_b <- ordinate(depression_rare, method = "PCoA", distance = bray)
pcoa_j <- ordinate(depression_rare, method = "PCoA", distance = jaccard)
pcoa_u <- ordinate(depression_rare, method = "PCoA", distance = unifrac)
pcoa_w <- ordinate(depression_rare, method = "PCoA", distance = wunifrac)

gg_pcoab <- plot_ordination(depression_rare, pcoa_b, color = "BDI_category_antidepressant_use")+
  labs(col = "BDI Category and \nAntidepressant Use")+
  theme(legend.title = element_text(size = 8))+
  geom_mark_ellipse(aes(filter = PD != 9.160527))+
  scale_x_continuous(expand = expansion(mult = 0.26)) +  # Adds 10% space on x-axis
  scale_y_continuous(expand = expansion(mult = 0.2))    # Adds 10% space on y-axis
gg_pcoab

gg_pcoaj <- plot_ordination(depression_rare, pcoa_j, color = "BDI_category_antidepressant_use")+
  labs(col = "BDI Category and \nAntidepressant Use")+
  theme(legend.title = element_text(size = 8))+
  stat_ellipse(type = "norm")

gg_pcoau <- plot_ordination(depression_rare, pcoa_u, color = "BDI_category_antidepressant_use")+
  labs(col = "BDI Category and \nAntidepressant Use")+
  theme(legend.title = element_text(size = 8))+
  stat_ellipse(type = "norm")

gg_pcoaw <- plot_ordination(depression_rare, pcoa_w, color = "BDI_category_antidepressant_use")+
  labs(col = "BDI Category and \nAntidepressant Use")+
  theme(legend.title = element_text(size = 8))+
  stat_ellipse(type = "norm", level = 0.8)

ggsave(filename = "figures/gg_pcoab.png", 
       gg_pcoab,
       height = 4, width = 6)

ggsave(filename = "figures/gg_pcoaj.png", 
       gg_pcoaj,
       height = 4, width = 6)

ggsave(filename = "figures/gg_pcoau.png", 
       gg_pcoau,
       height = 4, width = 6)

ggsave(filename = "figures/beta/gg_pcoaw.png", 
       gg_pcoaw,
       height = 4, width = 6)


#### tax bar plot ####
plot_bar(depression_rare, fill="Phylum") 

depression_RA <- transform_sample_counts(depression_rare, function(x) x/sum(x))

mpt_phylum <- tax_glom(depression_RA, taxrank = "Phylum", NArm=FALSE)

gg_taxaA <- plot_bar(mpt_phylum, fill="Phylum") + 
  facet_wrap(.~Antidepressant_use, scales = "free_x")

gg_taxaB <- plot_bar(mpt_phylum, fill="Phylum") + 
  facet_wrap(.~BDI_category, scales = "free_x")

gg_taxaC <- plot_bar(mpt_phylum, fill="Phylum") + 
  facet_wrap(.~BDI_category_antidepressant_use, scales = "free_x")

ggsave("figures/plot_taxonomyA.png"
       , gg_taxaA
       , height=8, width =12)

ggsave("figures/plot_taxonomyB.png"
       , gg_taxaB
       , height=8, width =12)

ggsave("figures/plot_taxonomyC.png"
       , gg_taxaC
       , height=8, width =12)

