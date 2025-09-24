#Funguild script

#in bash:
conda create -n funguild python=3.8
conda activate funguild
conda install pandas


git clone https://github.com/UMNFuN/FUNGuild
cd FUNGuild/
  
python Guilds_v1.0.py -otu /Users/username/Documents/project/otu_table.txt -db fungi

# Preparation of OTU table in R

nopoolps.dada2.soil <- readRDS("/data/lastexpansion/danieang/data/trimmed/mergedPlates/nopoolps.dada2.soil.rds")
poolps.dada2.soil <- readRDS("/data/lastexpansion/danieang/data/trimmed/mergedPlates/poolps.dada2.soil.rds")
pspoolps.dada2.soil <- readRDS("/data/lastexpansion/danieang/data/trimmed/mergedPlates/pspoolps.dada2.soil.rds")

np2 <- readRDS("/data/lastexpansion/danieang/data/trimmed/mergedPlates/rg2.poolps.rds")

otu_table_data <- as.data.frame(otu_table(np2))
tax_table_data <- as.data.frame(tax_table(np2))

# Transpose the OTU table
otu_table_data <- as.data.frame(t(otu_table(np2)))
dim(otu_table_data)

# Create a concatenated taxonomy string
tax_table_data$taxonomy_string <- apply(tax_table_data, 1, function(x) {
  paste0("k__", x["Kingdom"], ";p__", x["Phylum"], ";c__", x["Class"], ";o__", x["Order"], ";f__", x["Family"], ";g__", x["Genus"], ";s__", x["Species"])
})

# Combine OTU/ASV IDs with their taxonomy
otu_tax_table <- cbind(otu_table_data, tax_table_data$taxonomy_string)
colnames(otu_tax_table)[ncol(otu_tax_table)] <- "taxonomy"

write.table(otu_tax_table, file = "funguild_input_raw.txt", sep = "\t", quote = FALSE, col.names = NA)

#Back in bash

conda activate funguild
cd FUNGuild/
  
python Guilds_v1.1.py -otu /data/lastexpansion/danieang/data/trimmed/mergedPlates/funguild_input_raw.txt -db fungi

#Result saved to '/data/lastexpansion/danieang/data/trimmed/mergedPlates/funguild_input_raw.guilds.txt'
#FunGuild tried to assign function to 18926 OTUs in '/data/lastexpansion/danieang/data/trimmed/mergedPlates/funguild_input_raw.txt'.
#FUNGuild made assignments on 6808 OTUs.


#Back in R
library(dada2)
library(phyloseq)

# Step 1: Read the FUNGuild data
funguild_data <- read.table("/data/lastexpansion/danieang/data/trimmed/mergedPlates/funguild_input_raw.guilds.txt", 
                            header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Step 2: Check the first few rows of the FUNGuild data
head(funguild_data)

# Step 3: Extract the OTU table and taxonomy table from the phyloseq object
otu_table_data <- as.data.frame(otu_table(np2))
tax_table_data <- as.data.frame(tax_table(np2))

# Step 4: Check column names in the FUNGuild data
colnames(funguild_data)

# Step 5: Extract OTU IDs from the FUNGuild output and phyloseq OTU table
otu_ids_funguild <- funguild_data$X  # Assuming "X" is the column that contains the OTU IDs from FUNGuild
otu_ids_phyloseq <- rownames(otu_table_data)

# Step 6: Add OTU IDs as a column in the phyloseq taxonomy table
tax_table_data$OTU_ID <- rownames(tax_table_data)

# Step 7: Merge the FUNGuild data (both 'Guild' and 'Trophic.Mode') with the phyloseq taxonomy table
merged_tax_table <- merge(tax_table_data, 
                          funguild_data[, c("X", "Guild", "Trophic.Mode")],  # Include both 'Guild' and 'Trophic.Mode'
                          by.x = "OTU_ID", by.y = "X", all.x = TRUE)

# Step 8: Ensure row names are OTU IDs in the merged taxonomy table
rownames(merged_tax_table) <- merged_tax_table$OTU_ID
merged_tax_table$OTU_ID <- NULL  # Remove the extra OTU_ID column

# Step 9: Convert the merged taxonomy table back to a phyloseq tax_table object
tax_table(np2) <- tax_table(as.matrix(merged_tax_table))

# Step 10: Check the first few rows of the updated taxonomy table
head(tax_table(np2))

# Step 11: Save the updated phyloseq object with FUNGuild data merged
saveRDS(np2, file = "np2.funguild.rds")

# Step 12: Load the saved phyloseq object
funguild <- readRDS("np2.funguild.rds")

funguild <- phyloseq_validate(funguild, remove_undetected = TRUE)

#Number of OTUs assigned a functional trait = 2908

root_traits <- subset_samples(funguild, Individual != "S")
root_traits <- phyloseq_validate(root_traits, remove_undetected = TRUE)

soil_traits <- subset_samples(funguild, Individual == "S")
soil_traits <- phyloseq_validate(soil_traits, remove_undetected = TRUE)

# Filter out "-" entries in Trophic.Mode
alldatroot2_filtered <- subset_taxa(root_traits, tax_table(root_traits)[, "Trophic.Mode"] != "-")

# Transform the abundance data to relative abundance
alldatroot2_relabund <- transform_sample_counts(alldatroot2_filtered, function(x) x / sum(x))

#  Aggregate the total relative abundance by Trophic.Mode
abundance_data <- psmelt(alldatroot2_relabund) %>%
  group_by(Trophic.Mode) %>%
  summarise(Total_Abundance = sum(Abundance)) %>%
  arrange(desc(Total_Abundance))

# Select the top 15 most abundant Trophic.Modes
top_15_trophic_modes <- abundance_data %>%
  top_n(15, Total_Abundance) %>%
  pull(Trophic.Mode)

#  Filter the phyloseq object to include only the top 15 Trophic.Modes
alldatroot2_top15 <- subset_taxa(alldatroot2_relabund, Trophic.Mode %in% top_15_trophic_modes)

# Convert Trophic.Mode to a factor, ordered by total relative abundance
psmelt_data <- psmelt(alldatroot2_top15)
psmelt_data$Trophic.Mode <- factor(psmelt_data$Trophic.Mode, 
                                   levels = abundance_data$Trophic.Mode[abundance_data$Trophic.Mode %in% top_15_trophic_modes])

# Create the plot with a smaller legend
png("/data/lastexpansion/danieang/plots/funguild_elevation_site.png", width = 12, height = 8, units = "in", res = 500)
p <- ggplot(psmelt_data, aes(x = site, y = Abundance, fill = Trophic.Mode)) +
  geom_bar(stat = "identity", position = "stack") +  # Stacked bar plot
  facet_wrap(~ elevation_adj, scales = "free_x") +  # Facet by elevation_adj
  labs(x = "site", y = "Abundance", fill = "Trophic.Mode") +  # Update x-axis label and fill label to "Guild.x"
  theme_minimal() +  # Clean minimal theme
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  # Rotate x-axis labels
  guides(fill = guide_legend(reverse = FALSE)) +  # Legend guide
  theme(legend.text = element_text(size = 6),  # Smaller text for legend
        legend.title = element_text(size = 8),  # Smaller title
        legend.key.size = unit(0.4, 'cm'))  # Smaller size for legend keys

print(p)

dev.off()

# Step 1: Filter out "-" entries in Guild.x
alldatroot2_filtered <- subset_taxa(root_traits, tax_table(root_traits)[, "Guild.x"] != "-")

# Step 2: Transform the abundance data to relative abundance
alldatroot2_relabund <- transform_sample_counts(alldatroot2_filtered, function(x) x / sum(x))

# Step 3: Aggregate the total relative abundance by Guild.x
abundance_data <- psmelt(alldatroot2_relabund) %>%
  group_by(Guild.x) %>%
  summarise(Total_Abundance = sum(Abundance)) %>%
  arrange(desc(Total_Abundance))

# Step 4: Select the top 15 most abundant Guild.x categories
top_15_guilds <- abundance_data %>%
  top_n(15, Total_Abundance) %>%
  pull(Guild.x)

# Step 5: Filter the phyloseq object to include only the top 15 Guild.x categories
alldatroot2_top15 <- subset_taxa(alldatroot2_relabund, Guild.x %in% top_15_guilds)

# Step 6: Convert Guild.x to a factor, ordered by total relative abundance
psmelt_data <- psmelt(alldatroot2_top15)
psmelt_data$Guild.x <- factor(psmelt_data$Guild.x, 
                              levels = abundance_data$Guild.x[abundance_data$Guild.x %in% top_15_guilds])

#step 7: Create the plot with a smaller legend
png("/data/lastexpansion/danieang/plots/guild_x_elevation_site_small_legend.png", width = 12, height = 8, units = "in", res = 500)
p <- ggplot(psmelt_data, aes(x = elevation_adj, y = Abundance, fill = Guild.x)) +
  geom_bar(stat = "identity", position = "stack") +  # Stacked bar plot
  facet_wrap(~ site, scales = "free_x") +  # Facet by elevation_adj
  labs(x = "site", y = "Abundance", fill = "Guild.x") +  # Update x-axis label and fill label to "Guild.x"
  theme_minimal() +  # Clean minimal theme
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  # Rotate x-axis labels
  guides(fill = guide_legend(reverse = FALSE)) +  # Legend guide
  theme(legend.text = element_text(size = 6),  # Smaller text for legend
        legend.title = element_text(size = 8),  # Smaller title
        legend.key.size = unit(0.4, 'cm'))  # Smaller size for legend keys

print(p)

dev.off()




# Extract the tax_table and otu_table from your phyloseq object
tax_table_data <- as.data.frame(tax_table(nopoolps.soil.funguild))
otu_table_data <- as.data.frame(otu_table(nopoolps.soil.funguild))

# Add OTU_ID as a column for merging
tax_table_data <- rownames_to_column(tax_table_data, var = "OTU_ID")
otu_table_data <- rownames_to_column(otu_table_data, var = "OTU_ID")

tax_table_data <- tax_table_data %>%
  mutate(Guild = str_trim(Guild, side = "both"))

# Filter out entries where 'Guild' is exactly "-"
tax_table_filtered <- tax_table_data %>%
  filter(Guild != "-")
 
# Verify that '-' guilds are removed
print("Unique guilds after filtering '-':")
print(unique(tax_table_filtered$Guild))

# Proceed to join the OTU table with the filtered tax table
guild_abundances <- otu_table_data %>%
  left_join(tax_table_filtered, by = "OTU_ID") %>%
  filter(!is.na(Guild)) %>%  # Ensure that there are no NAs after the join
  group_by(Guild) %>%
  summarise(Total_Abundance = sum(across(where(is.numeric)))) %>%
  arrange(desc(Total_Abundance))

print(guild_abundances)

# Extract the top 12 most abundant guilds
top_12_guilds <- guild_abundances %>%
  slice(1:12) %>%
  pull(Guild)

print(top_12_guilds)

# Filter the phyloseq object to include only OTUs with the top 12 guilds
filtered_ps <- subset_taxa(nopoolps.soil.funguild, Guild %in% top_12_guilds)

# Create the plot
# Create the plot
guild_plot <- plot_bar(filtered_ps, fill = "Guild", x = "site") +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = "Top 12 Functional Guild Abundances Across Habitats",
       x = "Habitat",
       y = "Abundance") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels
        legend.text = element_text(size = 5))  # Adjust the legend text size

# Save the plot
output_file <- "/data/lastexpansion/danieang/plots/top_12_guild_abundance.pdf"

pdf(output_file, width = 10, height = 6)
print(guild_plot)
dev.off()

fungi_data <- read.table("/data/lastexpansion/danieang/data/trimmed/mergedPlates/funguild_input_raw.guilds.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

head(fungi_data)

# Filter out rows where Trophic.Mode or Confidence.Ranking are "-"
filtered_data <- fungi_data %>%
  filter(Trophic.Mode != "-", Confidence.Ranking != "-")

# Quick check of the filtered data
head(filtered_data)

# Summarize the top guilds after filtering
top_guilds <- filtered_data %>%
  group_by(Guild) %>%
  summarise(Count = n()) %>%
  arrange(desc(Count))

# View the top guilds
print(top_guilds)

# Summarize guilds by their confidence ranking after filtering
confidence_summary <- filtered_data %>%
  group_by(Confidence.Ranking) %>%
  summarise(Count = n()) %>%
  arrange(desc(Count))

# View the confidence summary
print(confidence_summary)

# Plotting top guilds
ggplot(top_guilds, aes(x = reorder(Trophic.Mode, -Count), y = Count)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = "Top Fungal Guilds", x = "Guild", y = "Count")

# Plotting guilds by confidence ranking
ggplot(confidence_summary, aes(x = reorder(Confidence.Ranking, -Count), y = Count)) +
  geom_bar(stat = "identity") +
  labs(title = "Guilds by Confidence Ranking", x = "Confidence Ranking", y = "Count")

#Creating object with only root samples

new_metadata <- read.csv("/data/lastexpansion/danieang/data/trimmed/mergedPlates/ITS_metadata_Soil.csv", row.names = 1)

sample_data(funguild) <- sample_data(new_metadata)
head(sample_data(funguild))

funguild_root <- subset_samples(funguild, Individual != "S")
funguild_root <- prune_taxa(taxa_sums(funguild_root) > 0, funguild_root)  # Prune taxa with zero abundance

#Now ploting agains "site"


tax_table(funguild_root)

# Convert Phyloseq object to a data frame with abundance, taxonomy, and metadata
fungi_df <- psmelt(funguild_root)

# Quick check of the data
head(fungi_df)

# Filter out rows where Guild is "-" and check for the "Site" column (adjust the column names if needed)
filtered_fungi_df <- fungi_df %>%
  filter(Guild != "-") %>%
  filter(!is.na(site))

# Summarize guild abundance by Site
guilds_by_site <- filtered_fungi_df %>%
  group_by(site, Guild) %>%
  summarise(Abundance = sum(Abundance)) %>%
  arrange(desc(Abundance))

# Identify the top 12 most abundant guilds across all sites
top_guilds_overall <- guilds_by_site %>%
  group_by(Guild) %>%
  summarise(Total_Abundance = sum(Abundance)) %>%
  arrange(desc(Total_Abundance)) %>%
  top_n(12, Total_Abundance) %>%
  pull(Guild)

# Filter the data to include only the top 12 guilds
top_guilds_by_site <- guilds_by_site %>%
  filter(Guild %in% top_guilds_overall)

# Plot the top 12 guilds by site
ggplot(top_guilds_by_site, aes(x = reorder(Guild, -Abundance), y = Abundance, fill = site)) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip() +
  labs(title = "Top 12 Fungal Guilds by Site", x = "Guild", y = "Abundance") +
  theme_minimal()
