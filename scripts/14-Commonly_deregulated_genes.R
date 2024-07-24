#' Commonly_deregulated_genes.R
#' 
#' This script generates the scatter plots of genes deregulated in different lines of the same effector,
#' tables with these data and barplot of number of genes deregulated in each line and
#' 
#' Outputs: 
#' - Plots:
#'	- plots/Arabidopsis_scatterPlot_DEGs\*DATE\*.pdf
#'	- plots/Poplar_scatterPlot_DEGs\*DATE\*.pdf
#'	- plots/Arabidopsis_numberDEGs\*DATE\*.pdf
#'	- plots/Poplar_numberDEGs\*DATE\*.pdf
#' - DEGs tables 
#'	- Arabidopsis/Mlp72983_deregulatedGenes.csv
#'	- Arabidopsis/Mlp52166_deregulatedGenes.csv
#'	- Poplar/Mlp72983_deregulatedGenes.csv
#'	- Poplar/Mlp52166_deregulatedGenes.csv

#### run the script to prepare the input datasets and load packages
source("scripts/loadAndPrepareDEGs_Files.R")
barplotTicks = list(Arabidopsis = seq(0, 1250, 250),
		    Poplar = seq(0, 4000, 1000))
#### Plot ####
for (species in SPECIES) {
	# Make one plot without the undesired plant line if necessary
	if (discarded_line %in% design[[species]]$plantLine) {
		# Remove the line from the table
		filter(degs_table[[species]],
		       plantLine != discarded_line) %>%
			# create a bar plot (do not put y)
			ggplot(aes(y = plantLine, fill = SampleName)) +
			geom_bar() +
			scale_x_continuous(breaks = barplotTicks[[species]]) +
			# geom_bar counts the number of times each plant line appears
			labs(x = "Number of deregulated genes",
			     y = "") +
			# Split the plot by the direction of the deregulation
			facet_wrap(~Deregulation, nrow = 2) + 
			# Tilt the text in the x
			 theme(legend.position = "bottom",
			 	# axis.text.x = element_text(angle = 45,
			 	#       		     vjust = 1,
			 	#       		     hjust = 1)
			 	)
		ggsave(paste0("plots/", species, "_numberDEGs_without_",
			      discarded_line, "_",
			      Sys.Date(), ".pdf"))
	}
	
	ggplot(degs_table[[species]], 
	       aes(y = plantLine, fill = SampleName)) +
		geom_bar() +
		labs(x = "Number of deregulated genes",
		     y = "") +
		scale_x_continuous(breaks = barplotTicks[[species]]) +
		facet_wrap(~Deregulation, nrow = 2) + 
		 theme(legend.position = "bottom",
		 	#axis.text.x = element_text(angle = 45,
		# 				 vjust = 1,
		# 				 hjust = 1))
		 	)
	 ggsave(paste0("plots/", species, "_numberDEGs_", Sys.Date(), ".pdf"))
	
}

#### Pivot wider by Sample ####
columnsToSelect = c("Gene", "plantLine", "log2FoldChange", "padj")
wideTable_DEGs_effector =
	sapply(simplify = F, SPECIES, \(species) {
		DEGS_TABLE = degs_table[[species]]
	sapply(unique(DEGS_TABLE$SampleName), simplify = F, \(effector) {
	       	data = filter(DEGS_TABLE, 
	       		      SampleName == effector)[, columnsToSelect]
	       
	       	output = 
	       		filter(data, padj < 0.01 &
	       		       	abs(log2FoldChange) > 2) %>%
	       		pivot_wider(id_cols = Gene,
	       			    names_from = plantLine,
	       			    values_from = log2FoldChange)
	       	outputFileName = paste0(
	       		species, "/", effector, "_deregulatedGenes.csv"
	       	)
	       	write.csv(output, file = outputFileName,
	       		  quote = F,row.names = F
	       	)
	       
	       	output	
	       })
	})


#### Function to draw scatter plot comparing lines expressing same effector ####
scatterPlot_DEGS_plantLines <- function(wideTable.DEGs) {
	# First we get the plantLine names
	plantLines = names(wideTable.DEGs)[2:3]
	# we filter for the genes that are deregulated in both lines
	common_degs = 
		!apply(wideTable.DEGs, 1, function(x) any(is.na(x)))
	wideTable.CommonDEGs = wideTable.DEGs[common_degs,]
	
	# Then we add columns called plantLine1 and plantLine2 to the tables
	# This will allow their selection in ggplot
	wideTable.CommonDEGs = 
		mutate(wideTable.CommonDEGs,
		       plantLine1 = wideTable.CommonDEGs[[plantLines[1]]],
		       plantLine2 = wideTable.CommonDEGs[[plantLines[2]]]
		) %>% mutate(
			color = case_when(
				plantLine1 > 0 &
					plantLine2 > 0 ~ "Up-regulated in both",
				plantLine1 < 0 &
					plantLine2 < 0 ~ "Down-regulated in both",
				.default = "Discordant deregulation"
			) %>% as.factor
		)
	# We also define if the gene is downregulated in both lines, 
	# up in both or shows discordant deregulation.
	if (length(levels(wideTable.CommonDEGs$color)) == 3) {
		pointColors = c("black", "blue", "red")
	} else {
		pointColors = c("blue", "red")
	}
	
	
	return(wideTable.CommonDEGs %>%
	       	ggplot(aes(x = plantLine1, y = plantLine2,
	       		   color = color)) +
	       	geom_point() +
	       	labs(x = plantLines[1],
	       	     y = plantLines[2]) +
	       	scale_color_manual(values = pointColors,
	       			   name = "")
	)
}

for (species in SPECIES) {
	pdf(paste0("plots/", species, "_scatterPlot_DEGs_", Sys.Date(), ".pdf"))
	lapply(names(wideTable_DEGs_effector[[species]]), \(effector) {
		cat(species, effector, "\n")
		(wideTable_DEGs_effector[[species]][[effector]] %>%
		scatterPlot_DEGS_plantLines() +
			theme(legend.position = "bottom") ) %>%
			print()
	})
	dev.off()
}
dev.off()
