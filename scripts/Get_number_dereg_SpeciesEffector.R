library(tidyverse)

links = c(Arabidopsis_521 = "Arabidopsis/Mlp52166_deregulatedGenes.csv",
	  Arabidopsis_729 = "Arabidopsis/Mlp72983_deregulatedGenes.csv",
	  Poplar_729 = "Poplar/Mlp52166_deregulatedGenes.csv",
	  Poplar_521 = "Poplar/Mlp72983_deregulatedGenes.csv"
)

sapply(links, \(x) {
	read.csv(x) %>% 
		apply(MARGIN = 1, FUN = \(rowN) !any(is.na(rowN))) %>%
		which %>%
		length
}) 
