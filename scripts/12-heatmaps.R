#### run the script to prepare the input datasets and load packages
source("scripts/loadAndPrepareDEGs_Files.R")

#### Check files
ls()
str(degs_table)

wideTable = 
	sapply(SPECIES, simplify = F,
	       \(species) {
	       	filtered = degs_table[[species]] %>%
	       		filter(abs(log2FoldChange) > 2 &
	       		       	padj < 0.01)
	       	columns = c("Gene", "log2FoldChange", "plantLine")
	       	wide = filtered[, columns] %>%
	       		pivot_wider(id_cols = Gene, 
	       			    names_from = plantLine,
	       			    values_from = log2FoldChange, 
	       			    values_fill = 0) %>%
	       		as.data.frame
	       	rownames(wide) = wide$Gene
	       	wide
	       })

genes.dendrogram = 
	sapply(SPECIES, simplify = F,
	       \(species) {
	       	
	       	data.set = wideTable[[species]][,-1]
	       	
	       	dist(data.set, method = "euclidian") %>%
	       		hclust()
	       })

plotTable = 
	sapply(SPECIES, simplify = F, \(species) {
		geneOrder = genes.dendrogram[[species]]$labels[
			genes.dendrogram[[species]]$order
		]
		data = wideTable[[species]] %>%
			pivot_longer(cols = !Gene,
				     names_to = "plantLine",
				     values_to = "FoldChange")
		data$Gene = 
			factor(data$Gene,
			       levels = geneOrder)
		data
	})
	
# remove grid from next plots
theme_set(theme_bw() +
	theme(panel.grid = element_blank(),
	      axis.ticks = element_blank()))

plotName = paste0("plots/Heatmaps_", Sys.Date(), ".pdf")

heatmaps.species = sapply(
	SPECIES, simplify = F, \(species) {
		totalGenesDeregulated = length(genes.dendrogram[[species]]$order)
		plotTable[[species]] %>%
			ggplot(aes(x = Gene, y = plantLine, fill = FoldChange)) +
			geom_tile() +
			scale_fill_gradient2(high = "red2",
					     low = "blue",
					     mid = "white",
					     na.value = "white",
					     name = "Fold change") +
			labs(y = "", x ="", title = species,
			     subtitle = paste0("Genes deregulated: ", 
			     		totalGenesDeregulated)) +
			theme(axis.text.x = element_blank())
	})

pdf(plotName, width = 8, height = 6)
ggarrange(plotlist = heatmaps.species,
	  ncol = 1, nrow = 2,
	  align = "v"
	  )
dev.off()
 
