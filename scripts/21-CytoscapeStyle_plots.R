# Install and load packages
options("repos" = c(CRAN = "https://cran.rstudio.com/",
		    BioCsoft = "https://bioconductor.org/packages/3.17/bioc/",
		    BioCann = "https://bioconductor.org/packages/3.17/data/annotation")
)

pkgs = c("clusterProfiler", "GOSemSim", 
	 "tidyverse", "igraph", "DOSE",
	 "org.At.tair.db", "ggraph", "ggrepel", "scatterpie")
pkgs.To.Install = ! pkgs %in% installed.packages()
if (any(pkgs.To.Install)) install.packages(pkgs[pkgs.To.Install])
for (curPkg in pkgs) library(curPkg, character.only = T) 

#### Common variables ####
SPECIES = c("Arabidopsis", "Poplar")
CSEPs = c("Mlp52166", "Mlp72983")
geneSets = c("Arabidopsis_Mlp52166",
	     "Arabidopsis_Mlp72983",
	     "Poplar_Mlp52166",
	     "Poplar_Mlp72983"#,
	     #"Mlp52166",
	     #"Mlp72983"
	     )

colorSets = c("#e07a5f",
	   "#3d405b",
	   "#81b29a",
	   "#f2cc8f"#,
	   # "red",
	   # "blue"
	   )

similarity_measure = "Lin"

#### Input files ####
load("RData_files/GO_annotation.RData")
load("RData_files/semantic_GO_At.Rdata")
load("RData_files/clusterProfiler_enrichmentResults.RData")

Ontologies = unique(GO_term_Info$Ontology)
select = dplyr::select

source("scripts/FUNCTIONS_GOplots.R")

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

pairwise_similarity = list()
for (deregulation in names(significant.enriched.set)) {
	# y is the data frame holding the data about the GO terms deregulated in each group

	pairwise_similarity[[deregulation]] = 
		sapply(Ontologies, simplify = F, \(onto) {
			get_PWSim(
				GO_Ontology_df = significant.enriched.set[[deregulation]], 
				#' Table must have at least an ID and an Ontology column
				ontology = onto, 
				semData = semData_At[[onto]],
				# We use the semantic data from Arabidopsis, 
				# since poplar does not have one
				simMeasure = similarity_measure
			)
		})
}	

longer_PW_similarity = list()
for (deregulation in names(significant.enriched.set)) {
	longer_PW_similarity[[deregulation]] =
		sapply(Ontologies, simplify = F, \(onto) {
			get_longPW_Sim(pairWS_sim = pairwise_similarity[[deregulation]][[onto]], 
				       GO_Ontology_df = significant.enriched.set[[deregulation]], 
				       ontology = onto, 
				       similarity_measure = similarity_measure)
		})
}

graph_GOs= list()
for (deregulation in names(significant.enriched.set)) {
	graph_GOs[[deregulation]] = 
		sapply(Ontologies, simplify = F, \(onto) {
			plotTitle = paste(deregulation, "-", 
					  gsub("_", " ", onto)
			)
			
			get_graphGOs(Ontology = onto,
				     GO_Ontology_df = significant.enriched.set[[deregulation]],
				     longDF_PWSim = longer_PW_similarity[[deregulation]][[onto]],
				     geneSetNames = geneSets,
				     color_vector = colorSets) +
				labs(title = plotTitle,
				     fill = "")
		})
}

pdf("plots/CytoscapeStyle_GOenrichment.pdf", 
    height = 8, width = 8)

for (deregulation in names(significant.enriched.set)) {
	for (onto in Ontologies) { 
		print(graph_GOs[[deregulation]][[onto]])
		}
}

dev.off()
	
