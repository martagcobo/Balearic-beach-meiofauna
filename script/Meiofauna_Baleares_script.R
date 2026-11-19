################################################################################
# BALEARIC ISLANDS MEIOFAUNA STUDY - STATISTICAL & BIOINFORMATIC ANALYSIS
# Script Author: Marta García-Cobo
# Last update: November 11th, 2026
# Paper: Local management shapes meiofaunal communities in Balearic sandy beaches 
################################################################################


# 0. SETUP & PACKAGE LOADING---------------------------------------------------------


# Set working directory
setwd("")

# Load required packages
library(dplyr)
library(tidyr)
library(stringr)
library(glmmTMB)
library(performance)
library(car)
library(lme4)
library(vegan)
library(emmeans)
library(sjPlot)
library(pairwiseAdonis)
library(ggplot2)
library(ggforce)

# Set seed for reproducibility
set.seed(123)



# 1. MORPHOLOGICAL DATA PREPARATION----------------------------------------------

# Load ecological covariates and community matrix
ecol <- read.csv2(here("Baleares samples.csv")) 
comm <- read.csv2("community_bal.csv") 

# Ensure identical sample ordering across datasets
ecol <- ecol[order(ecol$sample_ID), ]
comm <- comm[order(comm$sample_ID), ]

# Assign sample_ID as rownames and clean community matrix
rownames(ecol) <- ecol$sample_ID
rownames(comm) <- comm$sample_ID
comm <- comm[, -1] # Remove sample_ID column from count matrix

# Standardize continuous environmental variables (z-scores)
ecol$OM_s       <- scale(ecol$OM)
ecol$hii_s      <- scale(ecol$hii)
ecol$dist_pos_s <- scale(ecol$distance_pos)
ecol$grain_s    <- scale(ecol$grain_D50)

# Factor level formatting
ecol$protect_gradient <- factor(ecol$protect_gradient, levels = c("0", "1", "2", "3"))
ecol$beach_level      <- factor(ecol$beach_level, 
                                levels = c("upper", "middle", "lower"),
                                labels = c("Swash", "Shoaling", "Subtidal"))


# 2. UNIVARIATE STATISTICAL MODELS (GLMMs)-------------------------------------------

# --- 2.1. Taxonomic Alpha Diversity Model (Poisson) --------------------------
model_richness <- glmmTMB(taxa_richness ~ OM_s + grain_s + hii_s + 
                            beach_level + protect_gradient + (1|beach), 
                          data = ecol, 
                          family = poisson)

# Model diagnostics
check_overdispersion(model_richness)
check_collinearity(model_richness)
performance::check_model(model_richness)

# Model summary and Type II ANOVA
summary(model_richness)
car::Anova(model_richness)

# Post-hoc pairwise comparisons
emmeans(model_richness, pairwise ~ beach_level, type = "response")


# --- 2.2. Total Meiofaunal Density Model (Tweedie) ---------------------------
model_density <- glmmTMB(density_10 ~ OM_s + grain_s + hii_s + 
                             beach_level + protect_gradient + (1|beach), 
                           data = ecol, 
                           family = tweedie())

# Model diagnostics
check_overdispersion(model_density)
check_collinearity(model_density)
performance::check_model(model_density)

# Model summary and Type II ANOVA
summary(model_density)
car::Anova(model_density)

# Post-hoc pairwise comparisons
emmeans(model_density, pairwise ~ beach_level, type = "response")


# --- 2.3. Export Diversity & Density Model Summaries ------------------------
tab_model(model_density, model_richness, 
          show.intercept = TRUE, 
          transform = NULL,      
          show.se = TRUE,        
          show.stat = TRUE,      
          show.icc = TRUE,       
          show.r2 = TRUE,        
          string.est = "Estimate",
          p.style = "numeric",
          digits = 3,            
          file = "Models_richness_abundance.doc")


# --- 2.4. Nematoda / Copepoda Ratio Model (Gamma) ----------------------------
model_n_c <- glmmTMB(Nem_Cop ~ OM_s + grain_s + hii_s + 
                       protect_gradient + beach_level + (1|beach), 
                     data = ecol, 
                     family = Gamma(link = "log"))

# Model diagnostics
check_overdispersion(model_n_c)
check_collinearity(model_n_c)
performance::check_model(model_n_c)

# Model summary and Type II ANOVA
summary(model_n_c)
car::Anova(model_n_c)


# --- 2.5. Dominant Taxa Density Models (Tweedie) -----------------------------
area_core <- 36.32 # Core area in cm2

ecol$density_Nematoda   <- (comm$Nematoda / area_core) * 10
ecol$density_Copepoda   <- (comm$Copepoda / area_core) * 10
ecol$density_Annelida   <- (comm$Annelida / area_core) * 10
ecol$density_Turbellaria <- (comm$Turbellaria / area_core) * 10

# Nematoda
model_density_nematoda <- glmmTMB(density_Nematoda ~ OM_s + grain_s + hii_s + 
                                    beach_level + protect_gradient + (1|beach), 
                                  data = ecol, family = tweedie())

# Copepoda
model_density_copepoda <- glmmTMB(density_Copepoda ~ OM_s + grain_s + hii_s + 
                                    beach_level + protect_gradient + (1|beach), 
                                  data = ecol, family = tweedie())

# Annelida
model_density_annelida <- glmmTMB(density_Annelida ~ OM_s + grain_s + hii_s + 
                                    beach_level + protect_gradient + (1|beach), 
                                  data = ecol, family = tweedie())

# Turbellaria
model_density_turbellaria <- glmmTMB(density_Turbellaria ~ OM_s + grain_s + hii_s + 
                                       beach_level + protect_gradient + (1|beach), 
                                     data = ecol, family = tweedie())

# Export dominant taxa models summary table
tab_model(model_density_nematoda, 
          model_density_copepoda, 
          model_density_annelida, 
          model_density_turbellaria,
          show.intercept = TRUE, 
          transform = NULL,   
          show.se = TRUE,     
          show.stat = TRUE, 
          digits = 3,         
          digits.p = 3,
          string.est = "Estimate",
          p.style = "numeric",
          file = "Modelos_density_by_group.doc")


# 3. MULTIVARIATE COMMUNITY ANALYSES---------------------------------------------------

# --- 3.1. PERMANOVA Analysis (Bray-Curtis on Untransformed Matrix) -----------
beta_diversity <- vegdist(comm, method = "bray")

adonis <- adonis2(beta_diversity ~ OM_s + grain_s + hii_s + beach_level + protect_gradient, 
                  data = ecol, 
                  strata = ecol$beach, 
                  by = "margin")
print(adonis)

# --- 3.2. Sensitivity Check (Square-Root Transformed Matrix) -----------------
comm_sqrt <- sqrt(comm)
beta_sqrt <- vegdist(comm_sqrt, method = "bray")
adonis_sqrt <- adonis2(beta_sqrt ~ OM_s + grain_s + hii_s + beach_level + protect_gradient, 
                       data = ecol, 
                       strata = ecol$beach, 
                       by = "margin")
print(adonis_sqrt)

# Pairwise PERMANOVA across beach zones
pw_adonis_level <- pairwise.adonis(beta_diversity, ecol$beach_level, p.adjust.m = "bonferroni")
print(pw_adonis_level)


# --- 3.3. SIMPER Analysis -----------------------------------------------------
comm_density <- (comm / area_core) * 10
sim_global   <- simper(comm_density, ecol$beach_level)

# Utility function for cleaning SIMPER output
get_clean_simper_global <- function(sim_obj) {
  pairs <- names(sim_obj)
  all_pairs_df <- data.frame()
  
  for(p in pairs) {
    df <- as.data.frame(summary(sim_obj)[[p]])
    df$average <- round(df$average, 3)
    df$sd      <- round(df$sd, 3)
    df$ratio   <- round(df$ratio, 3)
    df$ava     <- round(df$ava, 2)
    df$avb     <- round(df$avb, 2)
    df$cumsum  <- round(df$cumsum, 3)
    df$p_formatted <- ifelse(df$p < 0.001, "< 0.001", as.character(round(df$p, 3)))
    df$Contrast <- p
    df$Taxon    <- rownames(df)
    all_pairs_df <- rbind(all_pairs_df, df)
  }
  return(all_pairs_df)
}

final_table_simper <- get_clean_simper_global(sim_global)
colnames(final_table_simper)[colnames(final_table_simper) == "ava"] <- "mean_density_groupA"
colnames(final_table_simper)[colnames(final_table_simper) == "avb"] <- "mean_density_groupB"

final_table_simper <- final_table_simper[, c("Contrast", "Taxon", "average", "sd", "ratio", 
                                             "mean_density_groupA", "mean_density_groupB", 
                                             "cumsum", "p_formatted")]

write.csv2(final_table_simper, "SIMPER_BeachZone_Density.csv", row.names = FALSE)



# 4. GRAPHICAL VISUALIZATIONS & FIGURES----------------------------------------------------------


# Custom taxonomic color palette
taxon_colors <- c(
  "Tardigrada" = "#90b293", "Acarina"    = "#488050", "Ostracoda"   = "#3d5d41",
  "Copepoda"   = "#1B4D3E", "Peracarida" = "#62a094", "Nematoda"    = "#197278",
  "Kinorhyncha"= "#fffef2", "Turbellaria"= "#EFEDA9", "Gastrotricha"= "#F3EB05",
  "Annelida"   = "#9DAD11", "Mollusca"   = "#664F02"
)

# Protection gradient color palette
colores_prot <- c("0" = "#FCF75F", "1" = "#6E6A02", "2" = "#488050", "3" = "#283618")

# Beach zone color palette
colores_niveles <- c("Swash" = "#A8E9ED", "Shoaling" = "#197278", "Subtidal" = "#283d3b")


# --- 4.1. Stacked Density Histogram across Beach Levels ----------------------
f_conv <- 0.275
ecol_with_comm <- cbind(ecol, comm)

ecol_long_global <- ecol_with_comm %>%
  pivot_longer(cols = any_of(names(taxon_colors)), names_to = "Taxon", values_to = "Raw_Count") %>%
  mutate(
    Density_Taxon = Raw_Count * f_conv,
    Taxon = factor(Taxon, levels = names(taxon_colors))
  )

ggplot(ecol_long_global, aes(x = beach_level, y = Density_Taxon, fill = Taxon)) +
  geom_bar(stat = "summary", fun = "mean", position = "stack", color = "black", linewidth = 0.1) + 
  scale_fill_manual(values = taxon_colors) +
  labs(y = expression(paste("Mean Density (ind. / 10 ", cm^2, ")")), 
       x = "Beach Zone", 
       fill = "Taxonomic Group") +
  theme_classic() +
  theme(axis.text.x = element_text(face = "bold", size = 12))


# --- 4.2. nMDS Ordination Plot -----------------------------------------------
nmds <- metaMDS(comm, distance = "bray", k = 2, trymax = 100, autotransform = FALSE)

fit_env  <- envfit(nmds, ecol[, c("OM", "grain_D50", "hii", "protect_gradient")], permutations = 999)
fit_taxa <- envfit(nmds, comm, permutations = 999)

zone_colors       <- c("#70C1B3", "#247BA0", "#1B365D")
zone_colors_trans <- adjustcolor(zone_colors, alpha.f = 0.35)
zone_shapes       <- c(21, 22, 24)

par(mar = c(5, 5, 4, 2))
plot(nmds, display = "sites", type = "n", main = "nMDS - Meiofaunal Community Structure")

ordiellipse(nmds, ecol$beach_level, kind = "se", conf = 0.95, 
            col = zone_colors_trans, border = zone_colors, draw = "polygon", lwd = 1.5)

points(nmds, display = "sites", 
       col = "black", 
       bg = zone_colors[as.numeric(ecol$beach_level)], 
       pch = zone_shapes[as.numeric(ecol$beach_level)], 
       cex = 0.8)

plot(fit_env, p.max = 0.05, col = "black", lwd = 2)
plot(fit_taxa, p.max = 0.05, col = "#4A4A4A", cex = 0.85, add = TRUE)

legend("topright", 
       legend = levels(ecol$beach_level), 
       col = "black", 
       pt.bg = zone_colors, 
       pch = zone_shapes, 
       pt.cex = 1.3, 
       bty = "n", 
       title = "Beach Zone")


# --- 4.3. Protection Gradient Violin & Boxplots ------------------------------
# Taxa Richness vs Protection Level
ggplot(ecol, aes(x = protect_gradient, y = taxa_richness, fill = protect_gradient)) +
  geom_violin(trim = FALSE, alpha = 0.6, color = NA) + 
  geom_boxplot(width = 0.1, color = "black", fill = "white", alpha = 0.3, outlier.shape = NA) +
  geom_jitter(aes(color = protect_gradient), width = 0.2, size = 2, alpha = 0.8) +
  geom_smooth(aes(x = as.numeric(protect_gradient)), method = "lm", color = "black", se = TRUE, linewidth = 0.8, fill = "gray80") +
  scale_fill_manual(values = colores_prot) +
  scale_color_manual(values = colores_prot) +
  labs(x = "Protection Gradient", y = "Taxa Richness", title = "Meiofauna Diversity across Protection Levels") +
  theme_classic() +
  theme(legend.position = "none", axis.text = element_text(size = 11), axis.title = element_text(size = 12, face = "bold"))

# Total Abundance vs Protection Level
ggplot(ecol, aes(x = protect_gradient, y = log(density_10), fill = protect_gradient)) +
  geom_violin(trim = FALSE, alpha = 0.6, color = NA) + 
  geom_boxplot(width = 0.1, color = "black", fill = "white", alpha = 0.3, outlier.shape = NA) +
  geom_jitter(aes(color = protect_gradient), width = 0.2, size = 2, alpha = 0.8) +
  geom_smooth(aes(x = as.numeric(protect_gradient)), method = "lm", color = "black", se = TRUE, linewidth = 0.8, fill = "gray80") +
  scale_fill_manual(values = colores_prot) +
  scale_color_manual(values = colores_prot) +
  labs(x = "Protection Gradient", y = expression(paste("Log Density (ind. / 10 ", cm^2, ")")), title = "Meiofauna Abundance across Protection Levels") +
  theme_classic() +
  theme(legend.position = "none", axis.text = element_text(size = 11), axis.title = element_text(size = 12, face = "bold"))


# --- 4.4. Beach Zone Violin & Sina Plots ------------------------------------
# Taxa Richness across Beach Levels
ggplot(ecol, aes(x = beach_level, y = taxa_richness, fill = beach_level)) +
  geom_violin(trim = FALSE, alpha = 0.6, color = NA) + 
  geom_boxplot(width = 0.1, color = "black", fill = "white", alpha = 0.3, outlier.shape = NA) +
  geom_sina(aes(color = beach_level), size = 1.8, alpha = 1) +
  geom_smooth(aes(x = as.numeric(beach_level)), method = "lm", color = "black", se = TRUE, linewidth = 0.8, fill = "gray80") +
  scale_fill_manual(values = colores_niveles) +
  scale_color_manual(values = colores_niveles) +
  labs(x = "Beach Zone", y = "Taxa Richness") +
  theme_classic() +
  theme(legend.position = "none", axis.title = element_text(face = "bold"))

# Abundance across Beach Levels
ggplot(ecol, aes(x = beach_level, y = log(density_10), fill = beach_level)) +
  geom_violin(trim = FALSE, alpha = 0.6, color = NA) + 
  geom_boxplot(width = 0.1, color = "black", fill = "white", alpha = 0.3, outlier.shape = NA) +
  geom_sina(aes(color = beach_level), size = 1.8, alpha = 1) +
  geom_smooth(aes(x = as.numeric(beach_level)), method = "lm", color = "black", se = TRUE, linewidth = 0.8, fill = "gray80") +
  scale_fill_manual(values = colores_niveles) +
  scale_color_manual(values = colores_niveles) +
  labs(x = "Beach Zone", y = expression(paste("Log Density (ind. / 10 ", cm^2, ")"))) +
  theme_classic() +
  theme(legend.position = "none", axis.title = element_text(face = "bold"))


# --- 4.5. Continuous Trends: Density vs. Sediment Organic Matter -------------
ggplot(ecol, aes(x = OM, y = density_10, color = beach_level, fill = beach_level)) +
  geom_point(alpha = 0.7, size = 2.6) + 
  geom_smooth(method = "lm", alpha = 0.2) +
  facet_wrap(~ beach_level) +
  scale_color_manual(values = colores_niveles) +
  scale_fill_manual(values = colores_niveles) +
  labs(title = "Trends in Density vs Sediment Organic Matter",
       x = "Sediment Organic Matter (%)", 
       y = expression(paste("Density (ind. / 10 ", cm^2, ")"))) +
  theme_bw() +
  theme(legend.position = "none", 
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(face = "bold"))


# --- 4.6. Density and Richness Boxplots by Beach ------------------------------

# Define explicit ordering for beach factor along X axis
orden_playas <- c("pujols", "calo", "comte", "bassa", "bou", "dor", "mitjana", "son_bou", "fornells")
ecol$beach   <- factor(ecol$beach, levels = orden_playas)

# Shape mapping for beach zones
formas_nivel <- c("Swash" = 15, "Shoaling" = 17, "Subtidal" = 16) # 15 = Square, 17 = Triangle, 16 = Circle

# Density Boxplot across individual beaches
p1_box <- ggplot(ecol, aes(x = beach, y = density_10, fill = protect_gradient)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8, color = "black") + 
  geom_jitter(aes(shape = beach_level), width = 0.2, size = 2.5, alpha = 0.8) + 
  scale_fill_manual(values = colores_prot, name = "Protection\nDegree") +
  scale_shape_manual(values = formas_nivel, name = "Beach\nZone") +
  labs(title = "Meiofaunal Density across Sampling Beaches",
       x = "Beach",
       y = expression(paste("Density (ind. / 10 ", cm^2, ")"))) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"))

# Richness Boxplot across individual beaches
p2_box <- ggplot(ecol, aes(x = beach, y = taxa_richness, fill = protect_gradient)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8, color = "black") +
  geom_jitter(aes(shape = beach_level), width = 0.2, size = 2.5, alpha = 0.8) +
  scale_fill_manual(values = colores_prot, name = "Protection\nDegree") +
  scale_shape_manual(values = formas_nivel, name = "Beach\nZone") +
  labs(title = "Taxonomic Richness across Sampling Beaches",
       x = "Beach",
       y = "Number of Taxa") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"))

print(p1_box)
print(p2_box)



# 5. METABARCODING (18S rDNA) ANALYSIS-------------------------------------------------

# --- 5.1. Load & Prepare Molecular Datasets ----------------------------------
ecol_meta <- read.csv2("Baleares_samples_meta.csv") 
species   <- read.csv2("Data18S_Baleares.csv")
comm_meta <- read.csv2("s0_DADA2.comm_baleares.csv")

# Filter taxonomy table for valid and verified meiofaunal ASVs
species <- species[species$njeval.result == "correct", ] 
species <- species[!duplicated(species$ASVid), ]       
species <- species[species$meio.tax == "TRUE", ]

# Set sample IDs as rownames
rownames(comm_meta) <- comm_meta$sample_ID

# Calculate global sequencing depth prior to filtering
sequencing_depth <- rowSums(comm_meta[, -c(1, 2)])

# Subset community matrix to meiofauna-confirmed ASVs
comm_meta <- comm_meta[, colnames(comm_meta) %in% species$ASVid]

# Read counts summary
total_meio_reads      <- sum(comm_meta)
meio_reads_per_sample <- rowSums(comm_meta)
total_reads_global    <- sum(sequencing_depth)

# Convert community matrix to binary presence/absence
comm_meta[comm_meta > 1] <- 1

# Align molecular community matrix with ecological metadata
common_samples <- intersect(rownames(comm_meta), ecol_meta$sample_ID)
comm_meta      <- comm_meta[common_samples, ]
ecol_meta      <- ecol_meta[match(common_samples, ecol_meta$sample_ID), ]
rownames(ecol_meta) <- ecol_meta$sample_ID

species   <- species[species$ASVid %in% colnames(comm_meta), ]
comm_meta <- comm_meta[, colnames(comm_meta) %in% species$ASVid]

# Add molecular depth metrics to metadata
ecol_meta$richness_meta      <- rowSums(comm_meta)
ecol_meta$seq_depth          <- sequencing_depth[match(ecol_meta$sample_ID, names(sequencing_depth))]
ecol_meta$meio_reads_sample  <- meio_reads_per_sample[match(ecol_meta$sample_ID, names(meio_reads_per_sample))]


# --- 5.2. Descriptive Molecular Statistics -----------------------------------
cat("\n================ METABARCODING SUMMARY ================\n")
cat("Total Reads:", total_reads_global, "\n")
cat("Total Meiofaunal Reads:", total_meio_reads, "\n")
cat("Meiofaunal Read Percentage:", round((total_meio_reads / total_reads_global) * 100, 2), "%\n")
cat("Total Meiofaunal ASVs:", ncol(comm_meta), "\n\n")

cat("ASVs per Major Taxonomic Group:\n")
print(table(species$blast.bestgroup))

cat("\nSample Richness Range (ASVs):\n")
cat("Max Richness:", max(ecol_meta$richness_meta, na.rm = TRUE), "\n")
cat("Min Richness:", min(ecol_meta$richness_meta, na.rm = TRUE), "\n")

# Parasite ASVs count
n_parasite_ASVs <- length(unique(species$ASVid[species$sistertrait.parasite == "1"]))
cat("Parasitic ASVs count:", n_parasite_ASVs, "\n")


# --- 5.3. Sequence Identity (% pident) Distribution ---------------------------
species$pident_num <- as.numeric(str_extract(species$blast.bestlabel, "(?<=pident=)[0-9.]+"))

cat("\n--- SEQUENCE IDENTITY (PIDENT) SUMMARY ---\n")
cat("Min Identity:", min(species$pident_num, na.rm = TRUE), "%\n")
cat("Max Identity:", max(species$pident_num, na.rm = TRUE), "%\n")
cat("Mean Identity:", round(mean(species$pident_num, na.rm = TRUE), 2), "%\n\n")

total_asvs <- sum(!is.na(species$pident_num))
p_mayor_97 <- sum(species$pident_num >= 97, na.rm = TRUE) / total_asvs * 100
p_90_a_97  <- sum(species$pident_num >= 90 & species$pident_num < 97, na.rm = TRUE) / total_asvs * 100
p_menor_90 <- sum(species$pident_num < 90, na.rm = TRUE) / total_asvs * 100

cat("--- THRESHOLD DISTRIBUTION ---\n")
cat(">= 97% Identity (Species/Genus level):", round(p_mayor_97, 2), "%\n")
cat("90% - 96.99% Identity (Family level):", round(p_90_a_97, 2), "%\n")
cat("< 90% Identity (Higher ranks / Unassigned):", round(p_menor_90, 2), "%\n")


# --- 5.4. ASV Richness across Beaches ----------------------------------------
cat("\n--- UNIQUE ASVs PER BEACH ---\n")
comm_playa <- aggregate(comm_meta, by = list(beach = ecol_meta$beach), FUN = sum)
rownames(comm_playa) <- comm_playa$beach
comm_playa <- comm_playa[, -1]

asvs_per_beach <- sort(rowSums(comm_playa > 0), decreasing = TRUE)
print(asvs_per_beach)


# --- 5.5. Checklist Generation ------------------------------------------------
comm_long <- as.data.frame(comm_meta) %>%
  mutate(sample_ID = rownames(comm_meta)) %>%
  pivot_longer(cols = -sample_ID, names_to = "ASVid", values_to = "Presence") %>%
  filter(Presence == 1)

checklist <- comm_long %>%
  left_join(species, by = "ASVid") %>%
  select(ASVid, sample_ID, blast.bestgroup, bestblast.family, bestblast.taxon) %>%
  filter(!is.na(bestblast.family) & bestblast.family != "")

write.csv2(checklist, "Checklist.csv", row.names = FALSE)


# --- 5.6. Taxonomic Heatmap Preparation ---------------------------------------
library(stringr)

df_heatmap <- checklist %>%
  left_join(species %>% select(ASVid, pident_num), by = "ASVid") %>%
  mutate(
    Presence = 1,
    clean_taxon = case_when(
      # 1. Species (>= 99%)
      pident_num >= 99 & !is.na(bestblast.taxon) & bestblast.taxon != "" ~ bestblast.taxon,
      
      # 2. Genus (97% - 98.99%)
      pident_num >= 97 & !is.na(bestblast.taxon) & bestblast.taxon != "" ~ paste0(word(bestblast.taxon, 1), " sp."),
      
      # 3. Family (90% - 96.99%)
      pident_num >= 90 & !is.na(bestblast.family) & bestblast.family != "" ~ paste0(bestblast.family, " sp."),
      
      # 4.  (< 90%)
      TRUE ~ paste0(blast.bestgroup, " sp.")
    ),
    asv_label = paste0(ASVid, " - ", clean_taxon),
    asv_num   = as.numeric(gsub("[^0-9]", "", ASVid))
  ) %>%
  complete(sample_ID, nesting(blast.bestgroup, bestblast.family, clean_taxon, asv_num, asv_label), 
           fill = list(Presence = 0)) %>%
  arrange(blast.bestgroup, bestblast.family, clean_taxon, asv_num) %>%
  mutate(asv_label = factor(asv_label, levels = unique(asv_label)))

p_final <- ggplot(df_heatmap, aes(x = sample_ID, y = asv_label, fill = as.factor(Presence))) +
  geom_tile(color = "grey80", linewidth = 0.2) + 
  scale_fill_manual(values = c("0" = "white", "1" = "#197278"), labels = c("Absence", "Presence")) +
  facet_grid(blast.bestgroup ~ ., scales = "free_y", space = "free_y") +
  labs(x = "Sample ID", y = "Amplicon Sequence Variants (ASVs)", fill = "Occurrence") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
    axis.text.y = element_text(size = 7),
    strip.text = element_text(face = "bold")
  )

print(p_final)



# --- 5.7. ASV Occurrence Heatmap ---------------------------------------------
df_heatmap <- checklist %>%
  left_join(species %>% select(ASVid, pident_num), by = "ASVid") %>%
  mutate(
    Presence = 1,
    clean_taxon = case_when(
      pident_num >= 97 & !is.na(bestblast.taxon) & bestblast.taxon != "" ~ bestblast.taxon,
      !is.na(bestblast.family) & bestblast.family != "" ~ paste0(bestblast.family, " sp."),
      TRUE ~ paste0(blast.bestgroup, " sp.")
    ),
    asv_label = paste0(ASVid, " - ", clean_taxon),
    asv_num   = as.numeric(gsub("[^0-9]", "", ASVid))
  ) %>%
  complete(sample_ID, 
           nesting(blast.bestgroup, bestblast.family, clean_taxon, asv_num, asv_label), 
           fill = list(Presence = 0)) %>%
  arrange(blast.bestgroup, bestblast.family, clean_taxon, asv_num) %>%
  mutate(asv_label = factor(asv_label, levels = unique(asv_label)))

p_final <- ggplot(df_heatmap, aes(x = sample_ID, y = asv_label, fill = as.factor(Presence))) +
  geom_tile(color = "grey80", linewidth = 0.2) + 
  scale_fill_manual(values = c("0" = "white", "1" = "#197278"),
                    labels = c("Absence", "Presence")) +
  facet_grid(blast.bestgroup ~ ., scales = "free_y", space = "free_y", switch = "y") +
  labs(fill = "Occurrence", x = "Beaches", y = "ASV Checklist (Conservative Taxonomic Assignment)") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_text(size = 7, face = "plain"),
    strip.text.y.left = element_text(angle = 0, face = "bold", size = 10),
    strip.background = element_rect(fill = "grey94", color = "grey80"),
    strip.placement = "outside",
    panel.spacing = unit(0.3, "lines"),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "grey80", fill = NA)
  )

print(p_final)



# 6. ABIOTIC ENVIRONMENTAL CHARACTERIZATION-----------------------------------------

# --- 6.1. Linear Models for Grain Size (D50) and Organic Matter (SOM) -------

ecol$protect_gradient <- as.numeric(as.character(ecol$protect_gradient))

# Median Grain Size (D50) Model
lm_grain <- lm(grain_s ~ beach_level + protect_gradient + hii_s, data = ecol)
car::Anova(lm_grain, type = "II")

# Sediment Organic Matter (SOM) Model
lm_som <- lm(OM_s ~ beach_level + protect_gradient + hii_s, data = ecol)
car::Anova(lm_som, type = "II")


# --- 6.2. Environmental Summary Statistics Table -----------------------------
env_summary <- ecol %>%
  group_by(beach_level) %>%
  summarise(
    SOM_mean = mean(OM, na.rm = TRUE),
    SOM_sd   = sd(OM, na.rm = TRUE),
    D50_mean = mean(grain_D50, na.rm = TRUE),
    D50_sd   = sd(grain_D50, na.rm = TRUE)
  )

print(env_summary)

coef(lm_grain)
coef(lm_som)


# 7. SESSION INFORMATION-----------------------------------------------------------

# > sessionInfo()
# R version 4.5.1 (2025-06-13 ucrt)
# Platform: x86_64-w64-mingw32/x64
# Running under: Windows 11 x64 (build 26200)
# 
# Matrix products: default
#   LAPACK version 3.12.1
# 
# locale:
# [1] LC_COLLATE=Spanish_Spain.utf8  LC_CTYPE=Spanish_Spain.utf8   
# [3] LC_MONETARY=Spanish_Spain.utf8 LC_NUMERIC=C                  
# [5] LC_TIME=Spanish_Spain.utf8    
# 
# time zone: Europe/Madrid
# tzcode source: internal
# 
# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#  [1] ggforce_0.5.0        stringr_1.5.1        tidyr_1.3.1          ggplot2_4.0.0       
#  [5] pairwiseAdonis_0.4.1 cluster_2.1.8.1      sjPlot_2.9.0         emmeans_1.11.2      
#  [9] vegan_2.7-1          permute_0.9-8        lme4_1.1-37          Matrix_1.7-3        
# [13] car_3.1-3            carData_3.0-5        performance_0.15.2   glmmTMB_1.1.11      
# [17] dplyr_1.1.4         
# 
# loaded via a namespace (and not attached):
#   [1] splines_4.5.1           palmerpenguins_0.1.1    bitops_1.0-9           
#   [4] tibble_3.3.0            R.oo_1.27.1             polyclip_1.10-7        
#   [7] datawizard_1.2.0        hardhat_1.4.1           pROC_1.18.5            
#  [10] rpart_4.1.24            lifecycle_1.0.4         Rdpack_2.6.4           
#  [13] fastcluster_1.3.0       doParallel_1.0.17       rprojroot_2.1.1        
#  [16] globals_0.18.0          lattice_0.22-7          MASS_7.3-65            
#  [19] insight_1.4.2           magrittr_2.0.3          minqa_1.2.8            
#  [22] RColorBrewer_1.1-3      lubridate_1.9.4         maps_3.4.3             
#  [25] abind_1.4-8             DHARMa_0.4.7            expm_1.0-0             
#  [28] R.cache_0.17.0          quadprog_1.5-8          purrr_1.1.0            
#  [31] R.utils_2.13.0          RCurl_1.98-1.17         nnet_7.3-20            
#  [34] pracma_2.4.4            tweenr_2.0.3            ipred_0.9-15           
#  [37] lava_1.8.1              listenv_0.9.1           terra_1.8-60           
#  [40] proto_1.0.0             parallelly_1.45.1       codetools_0.2-20       
#  [43] tidyselect_1.2.1        farver_2.1.2            effectsize_1.0.1       
#  [46] stats4_4.5.1            caret_7.0-1             e1071_1.7-16           
#  [49] ks_1.15.1               Formula_1.2-5           survival_3.8-3         
#  [52] iterators_1.0.14        foreach_1.5.2           TreeTools_1.15.0       
#  [55] tools_4.5.1             progress_1.2.3          Rcpp_1.1.0             
#  [58] glue_1.8.0              mnormt_2.1.1            prodlim_2025.04.28     
#  [61] xfun_0.52               here_1.0.2              mgcv_1.9-3             
#  [64] phytools_2.4-4          withr_3.0.2             numDeriv_2016.8-1.1    
#  [67] combinat_0.0-8          GGally_2.4.0            boot_1.3-31            
#  [70] pdist_1.2.1             digest_0.6.37           timechange_0.3.0       
#  [73] R6_2.6.1                estimability_1.5.1      colorspace_2.1-1       
#  [76] R.methodsS3_1.8.2       see_0.11.0              utf8_1.2.6             
#  [79] generics_0.1.4          Amelia_1.8.3            data.table_1.17.8      
#  [82] recipes_1.3.1           class_7.3-23            clusterGeneration_1.3.8
#  [85] prettyunits_1.2.0       parameters_0.28.2       scatterplot3d_0.3-44   
#  [88] ggstats_0.11.0          ModelMetrics_1.2.2.2    pkgconfig_2.0.3        
#  [91] gtable_0.3.6            timeDate_4041.110       S7_0.2.0               
#  [94] TMB_1.9.17              scales_1.4.0            gower_1.0.2            
#  [97] reformulas_0.4.1        geometry_0.5.2          knitr_1.50             
# [100] rstudioapi_0.17.1       reshape2_1.4.4          coda_0.19-4.1          
# [103] nlme_3.1-168            magic_1.6-1             nloptr_2.2.1           
# [106] proxy_0.4-27            sjlabelled_1.2.0        DEoptim_2.2-8          
# [109] KernSmooth_2.23-26      parallel_4.5.1          foreign_0.8-90         
# [112] pillar_1.11.0           grid_4.5.1              vctrs_0.6.5            
# [115] PlotTools_0.3.1         xtable_1.8-4            evaluate_1.0.4         
# [118] mvtnorm_1.3-3           cli_3.6.5               compiler_4.5.1         
# [121] rlang_1.1.6             crayon_1.5.3            future.apply_1.20.0    
# [124] labeling_0.4.3          mclust_6.1.1            plyr_1.8.9             
# [127] sjmisc_2.8.11           stringi_1.8.7           psych_2.5.6            
# [130] optimParallel_1.0-2     bayestestR_0.17.0       pacman_0.5.1           
# [133] hms_1.1.3               patchwork_1.3.2         bit64_4.6.0-1          
# [136] future_1.58.0           BAT_2.10.0              rbibutils_2.3          
# [139] igraph_2.1.4            hypervolume_3.1.6       phangorn_2.12.1        
# [142] fastmatch_1.1-6         nls2_0.3-4              bit_4.6.0              
# [145] ape_5.8-1