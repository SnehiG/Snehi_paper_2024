getSemData = function(orgDB = "org.At.tair.db", keytype = "TAIR", outputPath) {
	if (is.null(orgDB)) { orgDB = "org.At.tair.db"}
	if (is.null(keytype)) { keytype = "TAIR"}
	if (is.null(outputPath)) {
		orgname = gsub("org.([A-Za-z\\.]).db", "\\1", orgDB)
		outputPath = paste0("./", gsub("\\.", "_", orgDB), ".RData")
	}
	# start by creating a vector with all the packages you need
	## Set the websites from where R will download the packages
	options("repos" = c(CRAN = "https://cran.rstudio.com/",
			    BioCsoft = "https://bioconductor.org/packages/3.17/bioc/",
			    BioCann = "https://bioconductor.org/packages/3.17/data/annotation")
		
	)
	pkgs = c("GOSemSim", "tidyverse", orgDB)
	pkgs.To.Install = ! pkgs %in% installed.packages()
	if (any(pkgs.To.Install)) install.packages(pkgs[pkgs.To.Install])
	for (curPkg in pkgs) library(curPkg, character.only = T) 
	
	Ontologies = c("BP", "CC", "MF")
	semData = 
		sapply(Ontologies, simplify = F, \(onto) {
			ontology = sub("^(\\w)\\w+_(\\w)\\w+$",
				       "\\1\\2", onto) |> toupper()
			
			godata(orgDB, keytype = keytype, ont = ontology)
		})
	save(semData, file = outputPath)
}
