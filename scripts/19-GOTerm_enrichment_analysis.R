### GO enrichment with clusterProfiler ###
# Install and load packages
options("repos" = c(CRAN = "https://cran.rstudio.com/",
		    BioCsoft = "https://bioconductor.org/packages/3.17/bioc/",
		    BioCann = "https://bioconductor.org/packages/3.17/data/annotation")
)

pkgs = c("clusterProfiler", "tidyverse", "org.At.tair.db")
pkgs.To.Install = ! pkgs %in% installed.packages()
if (any(pkgs.To.Install)) install.packages(pkgs[pkgs.To.Install])
for (curPkg in pkgs) library(curPkg, character.only = T) 


## Load RData files 
load("RData_files/Gene_sets.RData")
load("RData_files/GO_annotation.RData")
SPECIES = c("Arabidopsis", "Poplar")
Ontologies = names(annotation.by.ontology)
deregulation = paste0(c("Up", "Down"), "-regulated")
similarityGO_cutoff = 0.3

##### GO enrichment
getEnrichmentsByGroups = function(deregulation, ontology) {
	set.names = names(geneSETs[[deregulation]])
	all_results = 
		sapply(set.names, simplify = F, \(SN) {
		if (grepl('_', SN)) {
			genes.in.set = geneSETs[[deregulation]][[SN]]
			species = str_split_i(SN, "_", 1)
		} else {
			species = "Poplar"
			genes.in.set = geneSETs[[deregulation]][[SN]]$Poplar_gene
		}
		annotation = annotation.by.ontology[[ontology]][[species]]
		
		organism = ifelse(species == "Poplar",
				  "Populus trichocarpa",
				  "Arabidopsis thaliana")
		result = enricher(gene = genes.in.set,
			 pvalueCutoff = 0.1,
			 pAdjustMethod = "BH",
			 TERM2GENE = annotation$TERM2GENE,
			 TERM2NAME = annotation$TERM2NAME
		)
		result@organism = organism
		result@ontology = gsub("^(\\w)\\w+_(\\w)\\w+$", "\\1\\2",
				       ontology) %>% toupper()
		
		clusterProfiler::simplify(
			x = result,
			cutoff = similarityGO_cutoff
		) %>% filter(p.adjust < .05, 
			     qvalue < 0.2, 
			     Count >= 2) %>%
			mutate(richFactor = 
			       	Count / as.numeric(sub("/\\d+", "", BgRatio)),
			       SetSize = sub("\\d+/", "", GeneRatio) %>%
			       	as.numeric
			       )
	})
	return(all_results)
}
	
enrichment_up = sapply(simplify = F, Ontologies, \(onto) {
	getEnrichmentsByGroups(
		deregulation = "Up-regulated", 
		ontology = onto
	)
})

enrichment_down = sapply(simplify = F, Ontologies, \(onto) {
	getEnrichmentsByGroups(
		deregulation = "Down-regulated", 
		ontology = onto
	)
})

save(enrichment_down, enrichment_up,
     file = "RData_files/clusterProfiler_enrichmentResults.RData")

