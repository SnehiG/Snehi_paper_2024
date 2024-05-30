### BEFORE RUNNING THIS SCRIPT ###
#' Check if the file "RData_files/GO_annotation.RData" and
#' "RData_files/clusterProfiler_enrichmentResults.RData"

# Install and load packages
options("repos" = c(CRAN = "https://cran.rstudio.com/",
		    BioCsoft = "https://bioconductor.org/packages/3.17/bioc/",
		    BioCann = "https://bioconductor.org/packages/3.17/data/annotation")
)

pkgs = c("clusterProfiler", "GOSemSim", "tidyverse", "venn")
pkgs.To.Install = ! pkgs %in% installed.packages()
if (any(pkgs.To.Install)) install.packages(pkgs[pkgs.To.Install])
for (curPkg in pkgs) library(curPkg, character.only = T) 

#' Gene sets
#' Genes up/down-regulated
#'	- Arabidopsis-Mlp52166 
#'	- Arabidopsis-Mlp72983
#'	- Poplar-Mlp52166 
#'	- Poplar-Mlp72983 
#' Genes up/down-regulated 
#'	- Mlp72983 (both species)
#'	- Mlp52166 (both species)
#' Genes concordantly deregulated in between all lines

#### Common variables ####
SPECIES = c("Arabidopsis", "Poplar")
CSEPs = c("Mlp52166", "Mlp72983")
excludeSets = CSEPs
#### Input files ####
load("RData_files/GO_annotation.RData")
load("RData_files/clusterProfiler_enrichmentResults.RData")
# The script below contains the functions:
# getTable_enrichedTerms, get_PWSim, get_longPW_Sim, get_graphGOs
source("https://gist.github.com/KarenGoncalves/95b6b81ef391d85b4fb5f8866cb448c3") 
Ontologies = unique(GO_term_Info$Ontology)

# Find significantly enriched terms
significant.enriched.set = 
	sapply(list("Down-regulated" = enrichment_down, 
		    "Up-regulated" = enrichment_up),
	       simplify = F, \(resultList) {
	       sapply(Ontologies,
	              simplify = F, \(onto) {
	              	getTable_enrichedTerms(resultList[[onto]])
	              }) %>% do.call(what = "rbind") %>% 
	       		as.data.frame
	       })

GroupNames = unique(significant.enriched.set[[1]]$SetName)

termEnriched_heatmap =
	sapply(c("Up-regulated", "Down-regulated"), 
	       simplify = F, \(dereg) {
	       	colorHigh = ifelse(grepl("Up", dereg),
	       			   "red", "blue")
	       	
	       	(significant.enriched.set[[dereg]] %>%
	       			filter(!SetName %in% excludeSets
	       			       ) %>%
	       			mutate(Ordered_names = 
	       			       	fct_reorder(Name, richFactor),
	       			       Ratio_GenesInTerm = Count/SetSize) %>%
	       			ggplot(aes(x = factor(SetName, levels = GroupNames),
	       				   y = Ordered_names,
	       				   size = Ratio_GenesInTerm,
	       				   color = richFactor)) +
	       			geom_point() +
	       			scale_color_gradient(low = "black",
	       					     high = colorHigh
	       			) +
	       			theme_minimal(base_size = 11) +
	       			labs(x = "", y = "", 
	       			     title = dereg,
	       			     size = "Genes in term\n(proportion)",
	       			     color = "Enrichment\nfactor"))
})

pdf(file = "plots/dotplot_enrichedGO.pdf", 
    width = 13, height = 9)
for (i in termEnriched_heatmap) {
	print(i)
}
dev.off()


########### FOR EFFECTORS ###############
termEnriched_heatmap_cseps =
	significant.enriched.set$`Down-regulated` %>%
	filter(SetName %in% c("Mlp72983",
			      "Mlp52166")
	) %>%
	mutate(Ordered_names = 
	       	fct_reorder(Name, richFactor),
	       Ratio_GenesInTerm = Count/SetSize) %>%
	ggplot(aes(x = factor(SetName, levels = GroupNames),
		   y = Ordered_names,
		   size = Ratio_GenesInTerm,
		   color = richFactor)) +
	geom_point() +
	scale_color_gradient(low = "black",
			     high = "blue"
	) +
	theme_minimal(base_size = 11) +
	labs(x = "", y = "", 
	     title = "Down-regulated",
	     size = "Genes in term\n(proportion)",
	     color = "Enrichment\nfactor")
ggsave("plots/dotplot_CSEPs_GO.pdf",
       height = 5.1, width = 7.3)



significant.enriched.set$`Down-regulated` %>%
	mutate(Ordered_names = 
	       	fct_reorder(Name, richFactor),
	       Ratio_GenesInTerm = Count/SetSize) %>%
	ggplot(aes(x = factor(SetName, levels = GroupNames),
		   y = Ordered_names,
		   size = Ratio_GenesInTerm,
		   color = richFactor)) +
	geom_point() +
	scale_color_gradient(low = "black",
			     high = "blue"
	) +
	theme_minimal(base_size = 8) +
	theme(axis.text.x = element_text(angle = 30,
					 vjust = 1,
					 hjust = 1)) +
	labs(x = "", y = "", 
	     title = "Down-regulated",
	     size = "Genes in term\n(proportion)",
	     color = "Enrichment\nfactor")
ggsave("plots/dotplot_allSets_GO.pdf",
       height = 7.1, width = 7.3)
