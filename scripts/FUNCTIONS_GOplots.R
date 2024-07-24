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

# Functions

getTable_enrichedTerms = function(enrichmentResultList) {
	set.names = names(enrichmentResultList)
	deregulated_GO = sapply(
		set.names, simplify = F, \(SN) {
			if (nrow(enrichmentResultList[[SN]]@result) > 0) {
				enrichmentResultList[[SN]]@result |>
					select(ID, 
					       GeneRatio,
					       Count, 
					       richFactor,
					       SetSize) |>
					mutate(SetName = SN)
			}
		}) %>% 
		do.call(what = "rbind") %>% as.data.frame
	
	resultTable = left_join(deregulated_GO, GO_term_Info,
				by = join_by("ID" == "GO"))
	return(resultTable)
}

get_PWSim = function(GO_Ontology_df, ontology, semData, simMeasure) {
	GOIDs = filter(y, Ontology == ontology)$ID %>%
		unique
	if (length(GOIDs) == 0) { NULL }
	# do not continue if there are no GO terms for the ontology selected
	pairWS_sim = sapply(GOIDs, \(GOID1) {
		sapply(GOIDs, \(GOID2)
		       # Use the same vector in nested loops
		       # to get the paiwise similarity matrix
		       goSim(GOID1, GOID2, 
		             semData = semData, 
		             measure = simMeasure)
		)
	}) %>% rbind
	#' the first loop gives the rows of the matrix, 
	#' so we bind them row-wise to get the matrix
	
	return(pairWS_sim)
}

get_longPW_Sim = function(pairWS_sim, 
			  GO_Ontology_df, 
			  ontology, 
			  similarity_measure) {
	colSimilarity = paste0(similarity_measure, "_sim")
	dataGO = GO_Ontology_df[, c("Name", "ID")] |> unique()
	longPW_Sim = 
		pairWS_sim |>
		data.frame(GOID1 = rownames(pairWS_sim)) |>
		pivot_longer(cols = !GOID1,
			     names_to = "GOID2",
			     values_to = colSimilarity, 
			     values_drop_na = T) |>
		mutate(GOID2. = gsub("\\.", ":", GOID2)) |>
		select(GOID1, GOID2., colSimilarity) |>
		rename(GOID2. = "GOID2") |>
		mutate(Description1 = 
		       	sapply(GOID1, simplify = T, \(x) dataGO$Name[dataGO$ID == x]) |>
		       	as.character(),
		       Description2 = 
		       	sapply(GOID2, simplify = T, \(x) dataGO$Name[dataGO$ID == x]) |>
		       	as.character()
		)
	return(longPW_Sim)
}


get_graphGOs <- function(Ontology, GO_Ontology_df, longDF_PWSim, geneSetNames, color_vector = NULL, legendName = "", edgeColour = NULL) {
	edgeColour = ifelse(is.null(edgeColour), 'darkgrey', edgeColour)
	
	GOs_inSets = 
		filter(GO_Ontology_df,
		       Ontology == Ontology,
		       SetName %in% geneSetNames)$Name |>
		unique()
	
	longDF_PWSim_filtered = 
		longDF_PWSim |>
		filter(Description1 %in% GOs_inSets,
		       Description2 %in% GOs_inSets)
	
	g <- graph.data.frame(longDF_PWSim_filtered[, 4:5], 
			      directed=FALSE)
	
	E(g)$width <- sqrt(longDF_PWSim_filtered[, 3] * 5)
	E(g)$weight <- longDF_PWSim_filtered[, 3]
	g <- delete.edges(g, E(g)[longDF_PWSim_filtered[, 3] < 0.1])
	pieChart_nodes = as_data_frame(g, "vertices")
	pieChart_nodes[, c("x", "y")] = 
		layout_with_graphopt(g, charge = 0.25)
	pieChart_nodes[, geneSetNames] <- 
		sapply(geneSetNames, \(SN) {
			enrichedGOs = filter(GO_Ontology_df, 
					     Ontology == Ontology,
					     SetName == SN)$Name
			c(pieChart_nodes$name %in% enrichedGOs) %>%
				unname %>%
				as.numeric()
			
		}) %>% as.data.frame()
	
	V(g)$x = pieChart_nodes$x;
	V(g)$y = pieChart_nodes$y
	
	p = ggraph(g, "manual", x = V(g)$x, y = V(g)$y) 
	
	if (length(E(g)$width) > 0) {
		p <- p + geom_edge_link(alpha=.8,
					colour=edgeColour)
	}
	
	p <- p + geom_scatterpie(cols = geneSetNames, 
				 data = pieChart_nodes,
				 colour = NA, pie_scale = 1)
	if (!is.null(color_vector)) {
		p <- p + scale_fill_manual(values = color_vector,
					   breaks = geneSetNames)
	}
		
	p <- p + geom_node_text(aes(label=name), 
				repel=TRUE, 
				bg.color = "white",
				force = 1.5)
	p + theme_void()
}
