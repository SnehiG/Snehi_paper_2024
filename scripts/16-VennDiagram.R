pkgs = c("tidyverse", "cluster", "ggpubr", "venn", "ggvenn", "ggVennDiagram")
pkgs.To.Install = ! pkgs %in% installed.packages()
if (any(pkgs.To.Install)) install.packages(pkgs[pkgs.To.Install])
for (curPkg in c(pkgs) ) library(curPkg, character.only = T) 

theme_bw() %>% theme_set()


#### VARIABLES ####
SPECIES = c("Arabidopsis", "Poplar")
CSEPs = c("Mlp52166", "Mlp72983")

annotateGenes = sapply(SPECIES, simplify = F,
		       \(species) {
		       	read.delim(paste0(species, "/DEGs_homologs.txt"))
		       })
deregulated_genes = sapply(SPECIES, simplify = F,
			   \(species) {
			   	read.delim(paste0(species, "/Deregulated_genesCSEPs.txt"))
			   })

#### Venn diagram ####
venn_data =
	full_join(annotateGenes$Arabidopsis[, -3], 
		  annotateGenes$Poplar[, -3],
		  by = c("Arabidopsis_gene", "Poplar_gene")
	) %>% 
	full_join(deregulated_genes$Arabidopsis, 
		  by = join_by(Arabidopsis_gene == Gene))  %>%
	full_join(deregulated_genes$Poplar, 
		  by = join_by(Poplar_gene == Gene))  %>%
	mutate(geneID = paste0(Arabidopsis_gene, " - ", Poplar_gene))

colsSelect = sapply(CSEPs, simplify = F, \(effector)
		    grep(effector, 
		         names(venn_data), 
		         value = T)
) %>% unlist %>% unname

for (i in colsSelect) {
	venn_data[, i] =
		venn_data[, i] %>%
		replace_na(replace = "Not_deregulated") %>% 
		factor(levels = c("Up-regulated",
				  "Down-regulated",
				  "Not_deregulated")
		)
}

upRegulated_genes = apply(venn_data, 1, \(DEG) {
	# get rows that have columns with value == "Up-regulated,
	# And no columns with "down-regulated"
	any(DEG == "Up-regulated")
}) %>% which
upRegulated_vennData = venn_data[upRegulated_genes,]

upRegulated_Lists = sapply(colsSelect, simplify = F,
			   \(species_csep) {
			   	selectRows = (
			   		upRegulated_vennData %>% dplyr::select(
			   			any_of(species_csep)) == 
			   			"Up-regulated")
			   	upRegulated_vennData$geneID[selectRows[, 1]]
			   })

### Venn down-regulated genes ###

downRegulated_genes = apply(venn_data, 1, \(DEG) {
	# get rows that have columns with value == "Down-regulated,
	# And no columns with "down-regulated"
	any(DEG == "Down-regulated")
}) %>% which
downRegulated_vennData = venn_data[downRegulated_genes,]

downRegulated_Lists = sapply(colsSelect, simplify = F,
			     \(species_csep) {
			     	selectRows = (
			     		downRegulated_vennData %>% 
			     			dplyr::select(
			     			any_of(species_csep)) == 
			     			"Down-regulated")
			     	genes = downRegulated_vennData$geneID[selectRows[, 1]]
			     	
			     })

#### Make plot

vennD = (Venn(upRegulated_Lists) %>% process_data())
vennD@region[, c("Up-regulated", "Down-regulated")] = 
	sapply(c("Up-regulated", "Down-regulated"), \(dereg) {
		N_dereg_genes = 
			sapply(SPECIES, \(speciesG) {
			geneNames = venn_data[,paste0(speciesG, "_gene")]
			
			sapply(vennD@region$name, \(vennZones) {
				columns_Selected = str_split_1(
					vennZones, 
					pattern = "\\.\\."
				)
				colsExcluded = colsSelect[which(!colsSelect %in% columns_Selected)]
				
				if (length(columns_Selected) == 1) {
					KEEP_ROWS = (venn_data[, columns_Selected] == dereg) &
						apply(venn_data[, colsExcluded], 1, \(rowV) {
							all(rowV != dereg)
						})
				} else if (length(colsExcluded) == 0) {
					KEEP_ROWS = apply(venn_data[, columns_Selected], 1, \(rowV) {
						all(rowV == dereg)
					})
				} else if (length(colsExcluded) == 1) {
					KEEP_ROWS = (venn_data[, colsExcluded] != dereg) &
						apply(venn_data[, columns_Selected], 1, \(rowV) {
							all(rowV == dereg)
						})
				} else {
					KEEP_ROWS = 
						apply(venn_data[, colsExcluded], 1, \(rowV) {
							all(rowV != dereg)
						}) & apply(venn_data[, columns_Selected], 1, \(rowV) {
							all(rowV == dereg)
						})
				}
				
				if (length(KEEP_ROWS) == 0) {
					""
				} else {
					geneNames[KEEP_ROWS] %>%
						unique %>% length
				}
			})
		})
		
		apply(N_dereg_genes, 1, \(x) paste(x, collapse = ", "))
		
	})
	
vennCorrected = read.csv("plots/Venn_data.csv", header = T)

vennD@region$`Up-regulated` = vennCorrected$Up.regulated %>%
	gsub(pattern = ", ", replacement = "  /  ")
vennD@region$`Down-regulated` = vennCorrected$Down.regulated  %>%
	gsub(pattern = ", ", replacement = "  /  ")

ggplot() +
	geom_sf(fill = NA, data = vennD@region, linewidth = .8) +
	geom_sf_text(aes(label = `Up-regulated`), 
		     color = "red",  nudge_y = 0.01, data = vennD@region,
		     size = 5) +
	geom_sf_text(aes(label = `Down-regulated`), 
		     color = "blue", nudge_y = -0.02, data = vennD@region,
		     size = 5) +
	geom_sf_text(aes(label = paste0(gsub("_", "\n", name), "\n")),
		     data = vennD@setLabel, nudge_y = -0.01,
		     size = 5) + 
	labs(caption = 
	     	"Genes are presented with their homologs in the format:\n(arabidopsis / poplar).") +
	theme_void() +
	theme(plot.caption = element_text(size = 12))

ggsave("plots/Venn_diagram_deregulatedGenes_homologs.pdf",
       height = 7, width = 9)
