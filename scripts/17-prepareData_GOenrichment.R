# start by creating a vector with all the packages you need
## Set the websites from where R will download the packages
options("repos" = c(CRAN = "https://cran.rstudio.com/",
		    BioCsoft = "https://bioconductor.org/packages/3.17/bioc/",
		    BioCann = "https://bioconductor.org/packages/3.17/data/annotation")
	
)
pkgs = c("biomaRt", "tidyverse")
pkgs.To.Install = ! pkgs %in% installed.packages()
if (any(pkgs.To.Install)) install.packages(pkgs[pkgs.To.Install])
for (curPkg in pkgs) library(curPkg, character.only = T) 
SPECIES = c("Arabidopsis", "Poplar")

biomartDatasets = paste0(c("athaliana", "ptrichocarpa"),
			 "_eg_gene")
names(biomartDatasets) = SPECIES

expressedGenes = sapply(SPECIES, simplify = F, \(species) {
	read.delim(paste0(species, "/counts_filtered.tsv"),
		   header = T, row.names = 1) %>%
		rownames()
})

biomart = sapply(simplify = F, SPECIES, \(species) {
	useMart(biomart = "plants_mart",
		host = "https://plants.ensembl.org",
		dataset = biomartDatasets[[species]],
		version = "Ensembl Plants Genes 57"
	)
})

features = c("ensembl_gene_id", "description", 
	     "go_id", "name_1006", "namespace_1003")
renameCols = c("GENE", "description", 
	       "GO", "Name", "Ontology")

geneAnnotation = sapply(SPECIES, simplify = F, \(species) {
	tableData = getBM(attributes = features,
			  mart = biomart[[species]])
	names(tableData) = renameCols
	tableData = tableData %>%
		filter(GO != "",
		       GENE %in% expressedGenes[[species]])
	
	list(TERM2GENE = unique(tableData[, c("GO", "GENE")]),
	     TERM2NAME = unique(tableData[, c("GO", "Name")]),
	     TERM2ONTO = unique(tableData[, c("GO", "Ontology")])
	)
})

for (species in SPECIES) {
	missingInfo = which(geneAnnotation[[species]]$TERM2NAME$Name == "")
	geneAnnotation[[species]]$TERM2NAME[missingInfo,] =
		lapply(missingInfo, \(GO_id_r) {
			GO_id = geneAnnotation[[species]]$TERM2NAME$GO[GO_id_r]
			termName = go2term(GO_id)$Term
			
			data.frame(GO = GO_id,
				   Name = 
				   	ifelse(length(termName) == 0,
				   	       "", termName)
				   )
		}) %>% do.call(what = "rbind") %>%
		as.data.frame
	
	geneAnnotation[[species]]$TERM2NAME = 
		geneAnnotation[[species]]$TERM2NAME[
			which(geneAnnotation[[species]]$TERM2NAME$Name != ""),]
	
	
	missingInfo = which(geneAnnotation[[species]]$TERM2ONTO$Ontology == "")
	geneAnnotation[[species]]$TERM2ONTO[missingInfo,] =
		lapply(missingInfo, \(GO_id_r) {
			GO_id = geneAnnotation[[species]]$TERM2ONTO$GO[GO_id_r]
			onto = go2ont(GO_id)$Ontology
			
			data.frame(GO = GO_id,
				   Ontology = 
				   	ifelse(length(onto) == 0,
				   	       "", onto))
		}) %>% do.call(what = "rbind") %>%
		as.data.frame
	
	geneAnnotation[[species]]$TERM2ONTO = 
		geneAnnotation[[species]]$TERM2ONTO[
			which(geneAnnotation[[species]]$TERM2ONTO$Ontology != ""),]
}

GO_term_Info = sapply(SPECIES, simplify = F,
		      \(species) {
		      	merge(geneAnnotation[[species]]$TERM2NAME,
		      	      geneAnnotation[[species]]$TERM2ONTO,
		      	      by = "GO")
		      }) %>% do.call(what = "rbind") %>%
	unique() %>% filter(GO != "")

GO_term_Info$Ontology = ifelse(GO_term_Info$Ontology == "BP", "biological_process", 
			       ifelse(GO_term_Info$Ontology == "MF",
			              "molecular_function", GO_term_Info$Ontology)
			       )
Ontologies = GO_term_Info$Ontology |> unique()

annotation.by.ontology = sapply(Ontologies, simplify = F, \(onto) {
	sapply(SPECIES, simplify = F, \(species) {
		annotation_list = geneAnnotation[[species]]
		GO_onto = (annotation_list$TERM2ONTO %>%
			filter(Ontology == onto))$GO
		list(TERM2GENE = annotation_list$TERM2GENE %>%
		     	filter(GO %in% GO_onto),
		     TERM2NAME = annotation_list$TERM2NAME %>%
		     	filter(GO %in% GO_onto),
		     TERM2ONTO = annotation_list$TERM2ONTO %>%
		     	filter(GO %in% GO_onto))
	})
})

save(geneAnnotation, GO_term_Info, annotation.by.ontology,
     file = "RData_files/GO_annotation.RData")

source("https://gist.github.com/KarenGoncalves/e1b37865d5f6bdeff1c8c490ba8df963") # load function getSemData
getSemData(outputPath = "RData_files/semantic_GO_At.RData") 
