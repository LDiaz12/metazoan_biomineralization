library(ggplot2)
library(dplyr)
library(grid)
library(tidyverse)
library(here)
library(ggVennDiagram)
library(plotly)
here()
# read in biomineralization-associated proteins data set #
proteins <- read_csv("/Users/aryagothoskar/Desktop/Biomineralization Files/mollusk_proteins_9.csv")

mollusk_test <- proteins %>%
  filter(Species %in% c("B_stephanieae", "E_crispata", "M_gigas", "O_bimaculoides"))
mollusk_test

m_venn_list_1 <- apply(mollusk_test, 1, function(x) {
  names(x)[x == 1]
})

venn_mollusk <- ggVennDiagram(m_venn_list_1, label_alpha = 1, label = "count",
                            category.names = c( "B_stephanieae", "E_crispata", "M_gigas", "O_bimaculoides"),
                            set_color = c("darkblue", "black", "darkgreen", "darkred")) +
  scale_fill_distiller(palette = "Greens", direction = 1, name = "Protein Count") + 
  #labs(title = "Comparison of Biomineralization-Associated Proteins \n in Complexa and Robusta Clades") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_x_continuous(expand = expansion(mult = .2))
venn_mollusk

row_number <- 100
numbers <- as.numeric(proteins[row_number, ])
numbers <- numbers[!is.na(numbers)]

top10 <- sort(numbers, decreasing = TRUE)[1:11]

protein_names <- c("Peptidyl-prolyl cis-trans isomerase B [Crassostrea gigas]", 
                   "EGF and laminin G domain-containing protein-like", 
                   "MAM and LDL-receptor class A domain-containing protein 1", 
                   "Tyrosine-protein phosphatase Lar [Crassostrea gigas]", 
                   "PREDICTED: uncharacterized protein LOC105319624 isoform X1 [Crassostrea gigas]",
                   "hypothetical protein CGI_10008375 [Crassostrea gigas]", 
                   "Chitotriosidase-1 [Crassostrea gigas]",
                   "EGF-like domain-containing protein 2",
                   "Wnt inhibitory factor 1 [Crassostrea gigas]",
                   "Chitotriosidase-1  [Crassostrea gigas]",
                   "PREDICTED: nuclear receptor coactivator 6-like [Crassostrea gigas]")

proteins_top10 <- data.frame(
  protein = str_wrap(protein_names, width = 40),  # wraps long names
  value = top10
)

proteins_top10 <- proteins_top10[!duplicated(proteins_top10$protein), ]

ggplot(proteins_top10, aes(x = reorder(protein, value), y = value, fill = value)) +
  geom_bar(stat = "identity") +
  coord_flip() +                                  # horizontal bars = more room for labels
  scale_fill_distiller(palette = "Greens", direction = 1, name = "Species Count") +
  labs(title = "10 Most Common Proteins",
       x = "Protein",
       y = "Amount of Species") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8),
        plot.title = element_text(hjust = 0.5))

# STOP HERE
m_protein_counts <- proteins %>%
  filter(Presence == 1) %>%
  group_by(Species) %>%
  summarise(
    Count = n (),
    Proteins = paste(Protein, collapse = ",")
  )

##################################################################################

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

#ggsave(plot = venn_clade, bg = "white", filename = here("Outputs", "venn_clade.png"))


venn_data <- process_data(Venn(venn_list_1))
region_data <- venn_region(venn_data) |> as.data.frame()
items_by_region <- split(region_data$item, region_data$name)

# summarize protein data as bar chart counts # 

protein_long <- proteins_short %>%
  pivot_longer(
    cols = -Coral,
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
    values = c("0" = "white", "1" = "coral"), 
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


heat_map_int <- plot_ly( 
  protein_long, 
  x = ~Protein, 
  y = ~Coral, 
  z = ~Presence, 
  type = "heatmap", 
  colors = c("white", "coral"), 
  hovertemplate = paste(
    "<b>Species:</b> %{y}<br>",
    "<b>Protein:</b> %{x}<br>",
    "<b>Presence:</b> %{z}<extra></extra>"))
heat_map_int
