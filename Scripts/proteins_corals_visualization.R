library(ggplot2)
library(dplyr)
library(grid)
library(tidyverse)
library(here)
library(ggVennDiagram)
library(plotly)
library(pheatmap)

# read in biomineralization-associated proteins data set #
proteins <- read_csv(here("Data", "Biomin_Proteins_Corals.csv"))
proteins_short <- read_csv(here("Data", "proteins_shortlist.csv"))

protein_cols <- colnames(proteins_short)[colnames(proteins_short) != "Coral" & colnames(proteins_short) != "Clade"]

venn_list_1 <- proteins_short %>%
  group_by(Clade) %>%
  summarise(across(all_of(protein_cols), max)) %>%
  column_to_rownames("Clade")

str(venn_list_1)

venn_list_1 <- apply(venn_list_1, 1, function(x) {
  names(x)[x == 1]
})

venn_clade <- ggVennDiagram(venn_list_1, label_alpha = 1, label = "count",
                       category.names = c("Complexa", "Robusta"),
                       set_color = c("blue", "darkgreen")) +
  scale_fill_distiller(palette = "Reds", direction = 1, name = "Protein Count") + 
  labs(title = "Comparison of Biomineralization-Associated Proteins \n in Complexa and Robusta Clades") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_x_continuous(expand = expansion(mult = .2))
venn_clade

ggsave(plot = venn_clade, bg = "white", filename = here("Outputs", "venn_clade.png"), height = 8, width = 8)


venn_data <- process_data(Venn(venn_list_1))
region_data <- venn_region(venn_data) |> as.data.frame()
region_data

items_by_region <- split(region_data$item, region_data$name)
items_by_region
# summarize protein data as bar chart counts # 

protein_long <- proteins_short %>%
  pivot_longer(
    cols = -c(Coral, Clade),
    names_to = "Protein", 
    values_to = "Presence"
  ) %>%
  mutate(Protein_num = as.numeric(sub("P", "", Protein))
         ) %>%
  arrange(Protein_num) %>%
  mutate(Protein = factor(Protein, levels = unique(Protein)))
protein_long

protein_counts <- protein_long %>%
  filter(Presence == 1) %>%
  group_by(Coral) %>%
  summarise(
    Count = n (),
    Proteins = paste(Protein, collapse = ",")
  )

heat_map <- ggplot(protein_long,
                   aes(x = Protein, y = Coral, fill = factor(Presence))) + 
  geom_tile(color = "white") +
  scale_fill_manual(
    values = c("0" = "grey", "1" = "coral"), 
    name = "Protein\nPresence", 
    labels = c("Absent", "Present") 
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
    panel.grid = element_blank()
  ) + 
  labs(
    title = "Biomineralization-Associated Proteins Across 4 Coral Species", 
    x = "Protein", 
    y = "Coral Species"
  )
heat_map

### get top most abundant variables ### 

proteins_short2 <- proteins_short %>%
  select(-("Clade"))

coral_labels <- proteins_short2[["Coral"]]

protein_matrix <- proteins_short2 %>%
  select(-all_of("Coral")) %>%
  mutate(across(everything(), as.numeric)) %>%
  as.matrix()

rownames(protein_matrix) <- coral_labels

top_20 <- 20

col_sums <- colSums(protein_matrix, na.rm = TRUE)
top_proteins <- names(sort(col_sums, decreasing = TRUE)[1:top_20])

protein_top <- protein_matrix[, top_proteins]

top_20_heatmap <- pheatmap(
  protein_top,
  color          = c("grey90", "coral"),   # 0, 1
  breaks         = c(-0.5, 0.5, 1.5),       
  cluster_rows   = TRUE,
  cluster_cols   = TRUE,
  show_rownames  = TRUE,                    
  show_colnames  = TRUE,
  fontsize_row   = 8,
  fontsize_col   = 7,
  angle_col      = 45,
  main           = paste0("Top ", top_20, " Most Abundant Proteins"),
  border_color   = "black",
  legend_breaks  = c(0, 1),
  legend_labels  = c("Absent", "Present")
)

ggsave(plot = top_20_heatmap, bg = "white", filename = here("Outputs", "top_20_heatmap.png"), height = 8, width = 8)

