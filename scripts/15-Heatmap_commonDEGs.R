pkgs = c("tidyverse", "cluster", "ggpubr", "venn")
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
##### HEATMAP HOMOLOGS #####

# left_join finds the rows from the first table present in the second:
# Poplar homologs from arabidopsis DEGs that are deregulated in poplar
# Keeps all rows from the first table not found in the second.
geneState_DEGsArabidopsis_wide =
	inner_join(annotateGenes$Arabidopsis[, -3], 
		   deregulated_genes$Poplar, by = join_by(Poplar_gene == Gene)
	) %>% filter(
		Poplar_Mlp72983 != "Not_deregulated" |
			Poplar_Mlp52166 !=  "Not_deregulated"
	) %>%
	left_join(deregulated_genes$Arabidopsis, 
		  by = join_by(Arabidopsis_gene == Gene))  %>%
	mutate(geneID = paste0(Arabidopsis_gene, " - ", Poplar_gene))


# Using grep to get the columns with effector names
colsSelect = sapply(CSEPs, simplify = F, \(effector)
		    grep(effector, names(geneState_DEGsArabidopsis_wide), 
		         value = T)
) %>% unlist %>% unname

# looping through the vector above to replace NAs with "Not_deregulated"
# and transform the column from character to factor
for (i in colsSelect) {
	geneState_DEGsArabidopsis_wide[, i] =
		geneState_DEGsArabidopsis_wide[, i] %>%
		replace_na(replace = "Not_deregulated") %>% 
		factor(levels = c("Up-regulated",
				  "Down-regulated",
				  "Not_deregulated")
		)
}

gene_clusters = geneState_DEGsArabidopsis_wide[,3:ncol(geneState_DEGsArabidopsis_wide)] %>%
	data.frame(row.names = "geneID") %>%
	cluster::daisy(metric = "gower") %>%
	cluster::diana(diss = T)

geneState_DEGsArabidopsis_long = geneState_DEGsArabidopsis_wide %>%
	pivot_longer(cols = all_of(colsSelect),
		     names_to = "plant_effector", 
		     values_to = "Gene_state"
	)

geneState_DEGsArabidopsis_long$Gene_state = 
	geneState_DEGsArabidopsis_long$Gene_state %>%
	gsub(pattern = "Not_deregulated", replacement = "Not deregulated") %>%
	factor(levels = c("Up-regulated",
			  "Down-regulated")
	)

geneState_DEGsArabidopsis_long$geneID = factor(
	geneState_DEGsArabidopsis_long$geneID,
	levels = gene_clusters$order.lab
)
geneState_DEGsArabidopsis_long = 
	geneState_DEGsArabidopsis_long %>%
	mutate(Species = str_split_i(plant_effector, pattern = "_", i = 1),
	       CSEP = str_split_i(plant_effector, pattern = "_", i = 2))

geneState_DEGsArabidopsis_long %>%
	ggplot(aes(y = geneID,
		   x = CSEP,
		   fill = Gene_state)) + 
	geom_tile() + 
	scale_fill_manual(values = c("red","blue"),
			  breaks = c("Up-regulated",
			  	   "Down-regulated"),
			  na.value = "white",
			  name = "") + 
	labs(x = "", y = "") +
	facet_wrap(~Species) + 
	theme_minimal() + 
	# theme(axis.ticks = element_blank(),
	#       axis.text.y = element_blank(),
	#       legend.position = "bottom")
	theme(axis.text.x = element_text(size = 8),
	      axis.ticks = element_blank(),
	      axis.text.y = element_text(size = 6),
	      panel.grid = element_blank())

ggsave(filename = "plots/Heatmap_homologs_deregulated.pdf",
       height = 10, width = 10)

#### Venn diagram ####
venn_data =
	full_join(annotateGenes$Arabidopsis[, -3], 
		   deregulated_genes$Poplar, by = join_by(Poplar_gene == Gene)
	) %>% 
	full_join(deregulated_genes$Arabidopsis, 
		  by = join_by(Arabidopsis_gene == Gene))  %>%
	mutate(geneID = paste0(Arabidopsis_gene, " - ", Poplar_gene))

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
			   		downRegulated_vennData %>% dplyr::select(
			   			any_of(species_csep)) == 
			   			"Down-regulated")
			   	downRegulated_vennData$geneID[selectRows[, 1]]
			   })
pdf("plots/SimpleVenn_DEGs.pdf")
par(mfrow = c(2, 1),
    cex = 1.3)
venn(upRegulated_Lists %>% unname, 
     box = F, ellipse = T,  par = F, 
     zcolor = "style", opacity = 0.1, plotsize = 30,
     snames = gsub("_", " - ", names(upRegulated_Lists)))
title(main = "Up-regulated genes", line = -16)

venn(downRegulated_Lists %>% unname, 
     box = F, ellipse = T, par = F, 
     zcolor = "style", opacity = 0.1, plotsize = 30,
     snames = gsub("_", " - ", names(downRegulated_Lists)))
title(main = "Down-regulated genes", line = -16)

dev.off()
