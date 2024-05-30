### Prepare gene set RData file ###
# Install and load packages
options("repos" = c(CRAN = "https://cran.rstudio.com/",
		    BioCsoft = "https://bioconductor.org/packages/3.17/bioc/",
		    BioCann = "https://bioconductor.org/packages/3.17/data/annotation")
)

pkgs = c("clusterProfiler", "tidyverse", "venn")
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
minGeneCount = 5
gene_set_names = c(
	paste0(rep(SPECIES, each = 2), "_", CSEPs),
	CSEPs,
	"common_all")

#### Input files ####
deregulated_genes = sapply(SPECIES, simplify = F,
			   \(species) {
			   	read.delim(paste0(species, "/Deregulated_genesCSEPs.txt"))
			   })	

homologs = sapply(SPECIES, simplify = F,
		  \(species) {
		  	read.delim(paste0(species, "/DEGs_homologs.txt"))
		  })	

homologs_joined = full_join(
	homologs$Arabidopsis[, -3],
	homologs$Poplar[, -3],
	by = c("Arabidopsis_gene", "Poplar_gene")
)

#### Get gene sets

## Get sets of genes deregulated in each species*effector combination
#' Genes up/down-regulated
#'	- Arabidopsis-Mlp52166 
#'	- Arabidopsis-Mlp72983
#'	- Poplar-Mlp52166 
#'	- Poplar-Mlp72983 

geneSets_CSEPs_species = 
	sapply(c("Up-regulated", "Down-regulated"), 
	       simplify = F,
	       \(dereg) {
	       	sapply(gene_set_names[1:4], simplify = F,
	       	       \(species_csep) {
	       	       	species = str_split_i(species_csep, "_", 1)
	       	       	
	       	       	deregulated_genes[[species]]$Gene[
	       	       		deregulated_genes[[species]][, species_csep] == dereg
	       	       	] %>% unique
	       	       })
	       })


## Get sets of genes deregulated in both species by each effector
#' Genes up/down-regulated 
#'	- Mlp72983 (both species)
#'	- Mlp52166 (both species)

geneSets_CSEPs = sapply(
	c("Up-regulated", "Down-regulated"), simplify = F, \(dereg) { 
		sapply(CSEPs, simplify = F, \(csep) {
			geneSet_select = grep(
				csep, value = T,
				names(geneSets_CSEPs_species[[dereg]])
			)
			species = str_split_i(geneSet_select, "_", 1)
			
			perSpecies = sapply(1:2, simplify=F, \(x) {
				GS = geneSets_CSEPs_species[[dereg]] [[ geneSet_select[x] ]] %>%
					as.data.frame
				names(GS) = paste0(species[x], "_gene")
				left_join(GS, homologs_joined,
					  by = names(GS))
			})
			
			inner_join(perSpecies[[1]],
				   perSpecies[[2]],
				   by =  c("Arabidopsis_gene", "Poplar_gene")) %>%
				unique
		})
	})

common_all = sapply(
	c("Up-regulated", "Down-regulated"), simplify = F, \(dereg) { 
		#' If the number of genes deregulated by either effector in both species is < 10,
		#' No need to do GO enrichment for DEGs common between all effectors/species
		
		if ( all( sapply(geneSets_CSEPs[[dereg]], nrow) > minGeneCount ) ) {
			set1 = geneSets_CSEPs[[dereg]][[ CSEPs[1] ]]
			set2 = geneSets_CSEPs[[dereg]][[ CSEPs[2] ]]
			inner_join(set1, set2, by = c("Arabidopsis_gene", "Poplar_gene")) %>%
				unique
		}
	})

geneSETs = sapply(
	c("Up-regulated", "Down-regulated"), simplify = F, \(dereg) { 
		sets = list(
			geneSets_CSEPs_species[[dereg]],
			geneSets_CSEPs[[dereg]],
			common_all[[dereg]]) %>% unlist(recursive = F)
		notEmpty = sapply(names(sets), \(SN) {
			if (class(sets[[SN]]) == "character") {
				length(sets[[SN]]) >= minGeneCount
			} else {
				nrow(sets[[SN]]) >= minGeneCount
			}
		})
		
		sets[ gene_set_names %in% names(sets)[notEmpty] ]
	})

save(geneSETs, file = "RData_files/Gene_sets.RData")
