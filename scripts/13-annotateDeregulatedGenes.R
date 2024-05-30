#' annotateGenes_deregulatedByEffector.R
#' 
#' This script filters the degs_tables to get the genes that are either 
#' up or down-regulated in both lines of each effector,
#' then finds the annotation from ensembl, including:
#' gene name, gene description, 
#' GO ID (for biological process, molecular function and cellular compartment).
#'
#' Outputs: 
#' - Annotation files:
#'	- Arabidopsis/Annotation_commonDEGs.txt
#'	- Arabidopsis/DEGs_annotation.txt
#'	- Poplar/Annotation_commonDEGs.txt
#'	- Poplar/DEGs_annotation.txt
#' - Homologs (deregulated genes in one species and their homologs in the other) 
#'	- Arabidopsis/DEGs_homologs.txt
#'	- Poplar/DEGs_homologs.txt


#### run the script to prepare the input datasets and load packages
source("scripts/loadAndPrepareDEGs_Files.R")

#### Variables ####
effectors = c("Mlp72983", "Mlp52166")
biomartDatasets = paste0(c("athaliana", "ptrichocarpa"),
			 "_eg_gene")
names(biomartDatasets) = SPECIES

deregulated_genes = 
	sapply(simplify = F, SPECIES,
       \(species) {
       	data = degs_table[[species]] %>%
       		filter(abs(log2FoldChange) > 2,
       		       padj < 0.01 &
       		       	plantLine != discarded_line)
       	
        csep1 = data %>% filter(SampleName == effectors[1]) %>%
        	pivot_wider(id_cols = Gene,
        		    names_from = plantLine,
        		    values_from = log2FoldChange)
        csep1_common = data.frame(
        	csep = apply(csep1[,-1], 1, \(rowSelected) {
        		case_when(all(rowSelected > 2) ~ "Up-regulated",
        			  all(rowSelected < -2) ~ "Down-regulated",
        			  .default = "Not_deregulated")}),
        	Gene = csep1$Gene)
        names(csep1_common)[1] = paste0(species, "_", effectors[1])
        csep2 = data %>% filter(SampleName == effectors[2]) %>%
        	pivot_wider(id_cols = Gene,
        		    names_from = plantLine,
        		    values_from = log2FoldChange)
        csep2_common = data.frame(
        	csep = apply(csep2[,-1], 1, \(rowSelected) {
        		case_when(all(rowSelected > 2) ~ "Up-regulated",
        			  all(rowSelected < -2) ~ "Down-regulated",
        			  .default = "Not_deregulated")}),
        	Gene = csep2$Gene)
        names(csep2_common)[1] = paste0(species, "_", effectors[2])
        
        output = full_join(csep1_common, csep2_common, 
        	  by = "Gene")[, c(2, 1, 3)]
        output[, 2] = ifelse(is.na(output[, 2]),
        		     "Not_deregulated",
        		     output[, 2])
        output[, 3] = ifelse(is.na(output[, 3]),
        		     "Not_deregulated",
        		     output[, 3])
        write.table(x = output, 
        	    file = paste0(species, "/Deregulated_genesCSEPs.txt"), 
        	    append = F, quote = F, sep = "\t", 
        	    row.names = F, col.names = T)
        
        dereg = apply(output[, -1], 1, \(rowSelected)
        	      any(rowSelected != "Not_deregulated"))
        output[dereg, ]
        })


#### Create biomart objects 
biomart = sapply(simplify = F,
		 SPECIES,
		 \(species) {
		 	useMart(biomart = "plants_mart",
		 		host = "https://plants.ensembl.org",
		 		dataset = biomartDatasets[[species]]
		 	)
		 })

geneLists = sapply(simplify = F,
		   SPECIES,
		   \(species) {
		   	deregulated_genes[[species]]$Gene
		   })

attributes = sapply(SPECIES, simplify = F,
		    \(species) {
		    	features = c("ensembl_gene_id", "description", 
		    		     "go_id", "name_1006", "namespace_1003")
		    	homologs = c(
		    		"ensembl_gene_id",
		    		"eg_homolog_ensembl_gene",
		    		"eg_homolog_perc_id"
		    	)
		    	if (species == "Arabidopsis") {
		    		feature = c(features, 
		    			"external_gene_name")
		    		homologs[2:3] = paste0("ptrichocarpa_", homologs[2:3])
		    	} else { 
		    		homologs[2:3] = paste0("athaliana_", homologs[2:3])
		    	}
		    	list(features, homologs)
		    })

annotateGenes = 
	sapply(simplify = F,
	       SPECIES,
	       \(species) {
	       	annotation_feature =
	       		getBM(filters = "ensembl_gene_id",
	       		      values = geneLists[[species]], 
	       		      attributes = attributes[[species]][[1]],
	       		      mart = biomart[[species]]) 
	       	
		selectedCols = !grepl("external", names(annotation_feature))
	       	names(annotation_feature)[selectedCols] = c("Gene_id", "description", 
	       						  "GO_id", "GO_name", "GO_class")
	       	
	       	annotation_homolog = 
	       		getBM(filters = "ensembl_gene_id",
	       		      values = geneLists[[species]], 
	       		      attributes = attributes[[species]][[2]],
	       		      mart = biomart[[species]])
	       	
	       	newNames = paste0(
	       		(gsub("athaliana.+", "Arabidopsis", 
	       		      names(annotation_homolog)[2:3]) %>%
	       		 	gsub(pattern = "ptrichocarpa.+", replacement = "Poplar")), 
	       		c("_gene", "_perc_id"))
	       	names(annotation_homolog) = c(paste0(species, "_gene"),
	       					 newNames)
	       	
	       	write.table(annotation_homolog, file = paste0(species, "/DEGs_homologs.txt"),
	       		    sep = "\t", quote = F, append = F, row.names = F)
	       	write.table(annotation_feature, file = paste0(species, "/DEGs_annotation.txt"),
	       		    sep = "\t", quote = F, append = F, row.names = F)
	       	list(annotation = annotation_feature,
	       	     homologs = annotation_homolog)
	       	
	       })
