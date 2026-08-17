library(ggplot2)
library(dplyr)
library(grid)
library(tidyverse)
library(here)
library(ggVennDiagram)
library(plotly)
here()
# read in biomineralization-associated proteins data set #
proteins <- read_csv("/Users/aryagothoskar/Desktop/Biomineralization Files/brachiopod_proteins4.csv")    

brach_test <- proteins %>%
  filter(Species %in% c("L_anatina", "N_anomala", "T_transversa"))

b_venn_list_1 <- apply(brach_test, 1, function(x) {
  names(x)[x == 1]
})

venn_brach <- ggVennDiagram(b_venn_list_1, label_alpha = 1, label = "count",
                              category.names = c("L_anatina", "N_anomala", "T_transversa"),
                              set_color = c("darkblue", "black", "darkgreen" )) +
  scale_fill_distiller(palette = "Greens", direction = 1, name = "Protein Count") + 
  #labs(title = "Comparison of Biomineralization-Associated Proteins \n in Brachiopods) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_x_continuous(expand = expansion(mult = .2))
venn_brach
