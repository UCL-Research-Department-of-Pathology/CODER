library(CNomplexity)
CNomplexity.path = path.package('CNomplexity')
datapath = file.path(CNomplexity.path, 'data')
datafiles = list.files(datapath,full.names=TRUE)
datafile = datafiles[grep("bedpe",datafiles)]


chromothripsisRes = chromothripsis(
		bedpeFile=datafile,
			bedpeChromCol1=1, # bedpe chrom1 col
			bedpePosCol1=2, # bedpe pos1 col
			direction1col=7,
			bedpeChromCol2=4, # bedpe chrom2 col
			bedpePosCol2=5, # bedpe pos2col
			direction2col=8,
			bedpeSampleCol=9,
			sepbedpe=FALSE,
			bedpeHead=TRUE
			) 
