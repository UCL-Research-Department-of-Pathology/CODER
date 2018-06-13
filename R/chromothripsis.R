# ======================================================================================
#				TEST FUNCTIONS
# ======================================================================================

# function to test whether a single chromosome has exponential segment lengths from seg file
segLengthsExponential = function(seg,startCol=3,endCol=4,verbose=FALSE,diag=FALSE) 
  {
  breakpoints = sapply(1:(nrow(seg)-1),
                       FUN=function(i) (seg[i+1,startCol]+seg[i,endCol])/2)
  seglengths = diff(breakpoints) # segment lengths
  if(diag) return(seglengths)
  # are segment lengths exponentially distributed?
  test = ks.test(seglengths, "pexp", 1/mean(seglengths)) # p>0.05 indicates that segLengths fit exponential distr
  return(test$p.value) # p<0.05 indicates not exponential  (suggesting chromothripsis)
  }

# function to test whether a single chromosome has exponential breakpoint distances from bedpe
breakpointsExponential = function(bedpe,chrom,chromCol1=1,posCol1=2,chromCol2=4,posCol2=5,verbose=FALSE,diag=FALSE) 
  {
  # breakpoints
  breakpoints1 = bedpe[which(bedpe[,chromCol1]==chrom),posCol1]
  breakpoints2 = bedpe[which(bedpe[,chromCol2]==chrom),posCol2]
  breakpoints=c(breakpoints1,breakpoints2)
  # order breakpoints
  breakpoints = sort(breakpoints)
  diffs = diff(breakpoints)
  if(verbose)
	{
	plot(x=jitter(c(rep(1,length(diffs)),rep(2,length(diffs))),factor=0.5),
		y=c(diffs,rexp(length(diffs),1/mean(diffs))),
		col=c(rep("black",length(diffs)),rep("gray",length(diffs))),
		xlab="Obs/Exp",ylab="Distance between breakpoints")
	legend("topright",legend=c("Obs.","Exp"),pch=1,col=c("black","gray"))
	}
  if(diag) return(diffs)
  # are segment lengths exponentially distributed?
  test = ks.test(diffs, "pexp", 1/mean(diffs)) # p>0.05 indicates that segLengths fit exponential distr
  return(test$p.value<0.05) # low p (<0.05) indicates not exponential  (suggesting chromothripsis)
  }


# randomness of DNA fragment joins
# counts of +/+ -/- +/- -/+ should be random (1/4,1/4,1/4,1/4)
randomJoins = function(bedpe,direction1col=9,direction2col=10,svclassCol=NULL,verbose=FALSE,pThresh=0.8,diag=FALSE)
  {
if(!is.null(svclassCol))
	{
	counts = table(bedpe[,svclassCol])
	if(length(counts)<4) counts = c(counts,rep(0,4-length(counts)))
	} else {
  	joins = paste0(bedpe[,direction1col],bedpe[,direction2col])
	counts = table(joins)
	if(length(counts)<4) 
		{
		saveNames = names(counts)
		counts = c(counts,rep(0,4-length(counts)))
		#possNames = c("++","+-","-+","--")
		#names(counts) = c(saveNames,possNames[which(!possNames%in%names(counts))])
		}
	}
  if(diag) return(counts)
  if(verbose)
	{
	barplot(rbind(counts,rep(sum(counts)/4,4)),beside=TRUE,legend.text=c("Obs.","Exp."))
	print(counts)	
	}
  # goodness of fit test to multinomial
  test = chisq.test(counts,p=rep(0.25,4)) # p>0.05 indicates that counts fit multinomial distr
  return(test$p.value>pThresh) # high p (>0.05) indicates multinomial  (suggesting chromothripsis)
  }

# randomness of DNA fragment order
# two sides of each breakpoint should be random draws from all breakpoint positions
randomOrder = function(bedpe,chromCol1=1,posCol1=2,chromCol2=4,posCol2=5,nSims=1000,pThresh=0.8)
  {
  # breakpoints
  breakpoints1 = bedpe[,c(chromCol1,posCol1)]
  breakpoints2 = bedpe[,c(chromCol2,posCol2)]
  colnames(breakpoints1)=colnames(breakpoints2)=c("chrom","pos")
  breakpoints=rbind(breakpoints1,breakpoints2)
  # order breakpoints
  breakpoints = breakpoints[order(breakpoints[,1],breakpoints[,2]),]
  breakpoints = unique(breakpoints)
  # indices
  indices = apply(bedpe,MARGIN=1,FUN=function(x)
    {
    x = gsub(" ","",x)
      c(
      which(breakpoints[,1]==unlist(x[chromCol1])&
              as.numeric(breakpoints[,2])==as.numeric(x[posCol1])),
      which(breakpoints[,1]==unlist(x[chromCol2])&
              as.numeric(breakpoints[,2])==as.numeric(x[posCol2]))
      )})
  indicesScore = mean(abs(indices[2,]-indices[1,]))
  # monte carlo simulations
  sims = replicate(nSims,abs(diff(sample(nrow(breakpoints),2))))
  # p value
  pVal = sum(sims>indicesScore)/nSims # p>0.05 indicates random draw (suggesting chromothripsis)
  return(pVal>pThresh) # high p (>0.05) indicates random draw (suggesting chromothripsis)
  }

# ability to walk chromosome
walkChrom = function()
  {

  }

# ======================================================================================
#				SINGLE WRAPPER
# ======================================================================================

# combining p-values with fishers method
fishersMethod = function(Ps)
  {
  test=-2*sum(log(Ps))
  pchisq(test,df=2*length(Ps),lower.tail=FALSE)
  }

# run a single test
runSingle = function(bedpe,
	direction1col=9,direction2col=10,
	chromCol1=1,posCol1=2,
	chromCol2=4,posCol2=5,nSims=1000,
	startCol=3,endCol=4,svclassCol=NULL,pThresh=0.8,diag=FALSE)
	{
	dobedpe = nrow(bedpe)>0
	if(dobedpe)
		{
		# check for random joins
		P2 = randomJoins(bedpe,
			direction1col=direction1col,
			direction2col=direction2col,
			svclassCol=svclassCol,pThresh=pThresh,diag=diag)
		# check for random selection of breakpoints
		P3 = randomOrder(bedpe,
			chromCol1=chromCol1,
			posCol1=posCol1,
			chromCol2=chromCol2,
			posCol2=posCol2,
			nSims=nSims,pThresh=pThresh)
		if(diag) return(list(P2,P3))
		return(c(P2,P3))
		} else {
		return(NA)
		}
	}


# ======================================================================================
#				WINDOW WRAPPER
# ======================================================================================

# split into windows, then check for chromothripsis
splitWindow = function(bedpe,chrom,size=3e7,gap=1e6,
	chromCol=2,startCol=3,endCol=4,chromCol1=1,posCol1=2,
	chromCol2=4,posCol2=5,direction1col=9,direction2col=10,
	breaksLimit=30,pThresh=0.8,chromStart=0,chromEnd=2e8,
	svclassCol=NULL,diag=FALSE)
	{
	if(nrow(bedpe)<breaksLimit) return(NA) # lower limit on number of fusions
	# p value for exponential distribution of breakpoints
	P1 = breakpointsExponential(bedpe,
			chrom=chrom,
			chromCol1=chromCol1,
			posCol1=posCol1,
			chromCol2=chromCol2,
			posCol2=posCol2,diag=diag) 
	# get windows of size
	split = seq(from=min(chromStart),to=max(chromEnd)-size,by=gap)
	if((max(split)+size)<max(chromEnd)) 
		{
		split = c(split,split[length(split)]+gap)
		names(split)[length(split)] = split[length(split)]
		}
	# get Granges
	bedpeGrange1 = as(paste0("chr",
			bedpe[,chromCol1],":",
			bedpe[,posCol1],"-",
			bedpe[,posCol1]),
		"GRanges")
	bedpeGrange2 = as(paste0("chr",
			bedpe[,chromCol2],":",
			bedpe[,posCol2],"-",
			bedpe[,posCol2]),
		"GRanges")
	# check for chromothripsis in each window
	res = sapply(1:length(split),FUN=function(i)
		{
		checkGrange = as(paste0("chr",chrom,":",
					split[i],"-",
					split[i]+size),
				"GRanges")
		bedpeIndex1 = findOverlaps(bedpeGrange1,checkGrange)
		bedpeIndex2 = findOverlaps(bedpeGrange2,checkGrange)
		bedpeIndex = unique(bedpeIndex1@from,bedpeIndex2@from)
		P = runSingle(bedpe=bedpe[bedpeIndex,],
			startCol=startCol,
			endCol=endCol,
			chromCol1=chromCol1,
			posCol1=posCol1,
			chromCol2=chromCol2,
			posCol2=posCol2,
			direction1col=direction1col,
			direction2col=direction2col,
			svclassCol = svclassCol,	
			pThresh=pThresh,diag=diag)
		if(diag) return(P)
		if(!any(is.na(P))) 
			{
			#return(fishersMethod(c(P1,P))) #fishers method
			return(sum(c(P1,P))>1) # cutoff method - 2 tests must pass cutoff
			} else {
			return(NA)
			}
		})
	if(diag) return(res)
	if(all(is.na(res))) return(P1)
	res[which(is.na(res))] = FALSE
	names(res) = split
	return(res)
	}


# ======================================================================================
#				OVERALL WRAPPER
# ======================================================================================

# function to combine overlapping regions
combineRegions = function(regions,sampleCol=1,chromCol=2,startCol=3,endCol=4)
	{
	if(nrow(regions)<2) return(regions)
	doOverlap=TRUE
	while(doOverlap)
		{
		seqStrings = paste0("chr",regions[,chromCol],":",regions[,startCol],"-",regions[,endCol])
		seqGranges = as(seqStrings,"GRanges")
		overlaps = findOverlaps(seqGranges,seqGranges)
		doOverlap = length(overlaps@from)!=nrow(regions)
		if(!doOverlap) break
		overlapInfo = sapply(1:length(seqStrings),FUN=function(x) overlaps@to[which(overlaps@from==x)],simplify=FALSE)
		newRegions = t(sapply(overlapInfo,FUN=function(x) c(unique(regions[x,sampleCol]),
					unique(regions[x,chromCol]),
					min(as.numeric(regions[x,c(startCol,endCol)])),
					max(as.numeric(regions[x,c(startCol,endCol)]))
					)))
		regions = unique(newRegions)
		}
	return(regions)
	}

# function to get runs and clean output
getRuns = function(chromScores,chrom,samp,size)
	{
	if(length(chromScores)>1)
		{
		#chromBool = chromScores<0.05
		chromBool = chromScores
		if(!any(chromBool)) return(NULL)
		runs = rle(chromBool)
		ends = cumsum(runs$lengths)
		starts = c(1,ends[-length(ends)]+1)
		windowStarts = as.numeric(names(chromBool)[starts[which(runs$values==TRUE)]])
		windowEnds = as.numeric(names(chromBool)[ends[which(runs$values==TRUE)]])
		windowEnds = windowEnds+size
		windows = cbind(samp,chrom,windowStarts,windowEnds)
		windows = combineRegions(windows) 
		return(windows)
		} else {
		#return(cbind(samp,chrom))
		# only return results if have looked at translocations
		return(NULL)
		}
	}




# function to run whole chromothripsis analysis
chromothripsis = function(bedpeFile=NULL, # directory of separate bedpes, or single bedpe file
			size=3e7, # window size
			gap=1e6, # gap between sliding windows
			bedpeChromCol1=1, # bedpe chrom1 col
			bedpePosCol1=2, # bedpe pos1 col
			direction1col=7, # orientation of first partner
			bedpeChromCol2=4, # bedpe chrom2 col
			bedpePosCol2=5, # bedpe pos2col
			direction2col=8, # orientation of second partner
			bedpeSampleCol=9, # bedpe sample col
			svclassCol=NULL,
			doParallel=FALSE, # parallel computation
			nCores = NULL,	# number of cores
			samplesToRun = NULL, # which samples to run
			chromsToRun = NULL, # which chromosomes to run
			sepbedpe=TRUE, # Are bedpes separate
			bedpeHead=FALSE, # does bedpe file have header
			bedpeEnding=".brass.annot.bedpe.gz",
			breaksLimit=30, #  minimum number of breakpoints on chromosomes
			pThresh=0.8, # p value threshold for tests
			cytoFile=NULL,
			diag=FALSE
			) 
	{
	if(doParallel&is.null(nCores)) nCores = detectCores()
	# read in bedpe
	if(!sepbedpe)
		{
		allbedpe = readFile(bedpeFile,bedpeHead)
		samples = unique(allbedpe[,bedpeSampleCol])
		} else {
		samples = list.files(bedpeFile)
		}
	# get chromomsome start/ends
	chromInfo = getChromInfo(cytoFile)
	# run analysis per sample per chromosome
	# loop over samples
	Ps = sapply(samples,FUN=function(y)
		{
		print(y)
		# load bedpe for this sample
		if(sepbedpe)
			{
			bedpe = readFile(paste0(bedpeFile,"/",y),bedpeHead)
			} else {
			bedpe = allbedpe[which(paste0(allbedpe[,bedpeSampleCol])==paste0(y)),]
			}
		if(is.null(chromsToRun))
			{
			chromosomes = unique(c(paste0(bedpe[,bedpeChromCol1]),
						paste0(bedpe[,bedpeChromCol2])))
			} else {
			chromosomes = chromsToRun
			} 
		# chromothripsis calculation
		if(doParallel)
			{
			# loop over chromosomes
			res = mclapply(chromosomes,FUN=function(x)
				{
				print(x)
				# keep bedpe rows for this chms
				index1=paste0(bedpe[,bedpeChromCol1])==paste0(x)
				index2=paste0(bedpe[,bedpeChromCol2])==paste0(x)
				indexBedpe = which(index1|index2)
				if(length(indexBedpe)==0) return(NULL)
				# check for chromothripsis
				chromScores = splitWindow(bedpe=bedpe[indexBedpe,],
					chrom=paste0(x),
					size=size,
					chromCol=chromCol,
					startCol=startCol,
					endCol=endCol,
					chromCol1=bedpeChromCol1,
					posCol1=bedpePosCol1,
					chromCol2=bedpeChromCol2,
					posCol2=bedpePosCol2,
					direction1col=direction1col,
					direction2col=direction2col,
					breaksLimit=breaksLimit,pThresh=pThresh,
					chromStart=chromInfo[1,paste0(x)],
					chromEnd=chromInfo[2,paste0(x)])
				# just output regions that are chromothriptic
				getRuns(chromScores,paste0(x),paste0(y),size)},mc.cores=nCores)	
			} else {
			# loop over chromosomes
			res = sapply(chromosomes,FUN=function(x)
				{
				print(x)
				# keep bedpe rows for this chms
				index1=paste0(bedpe[,bedpeChromCol1])==paste0(x)
				index2=paste0(bedpe[,bedpeChromCol2])==paste0(x)
				indexBedpe = which(index1|index2)
				if(length(indexBedpe)==0) return(NULL)
				# check for chromothripsis
				chromScores = splitWindow(bedpe=bedpe[indexBedpe,],
					chrom=paste0(x),
					size=size,
					chromCol=chromCol,
					startCol=startCol,
					endCol=endCol,
					chromCol1=bedpeChromCol1,
					posCol1=bedpePosCol1,
					chromCol2=bedpeChromCol2,
					posCol2=bedpePosCol2,
					direction1col=direction1col,
					direction2col=direction2col,
					breaksLimit=breaksLimit,pThresh=pThresh,
					chromStart=chromInfo[1,paste0(x)],
					chromEnd=chromInfo[2,paste0(x)],
					diag=diag)
				if(diag) return(chromScores)
				# just output regions that are chromothriptic
				getRuns(chromScores,paste0(x),paste0(y),size)},simplify=FALSE)
			}
		names(res) = chromosomes
		return(res)
		},simplify=FALSE)
	names(Ps) = samples
	return(Ps)
	}



# get chromosome lengths
getChromInfo = function(cytoFile=NULL)
	{
	data = loadCytoBand(cytoFile)
	data[,1] = gsub("chr","",data[,1])
	sapply(c(1:22,"X","Y"),FUN=function(x) 
		c(min(data[which(data[,1]==x),2]),
		max(data[which(data[,1]==x),3])))
	}


# ======================================================================================
#				SIMULATE CHROMOTHRIPSIS
# ======================================================================================

# function to create a list of SV breakpoints from simulated chromothripsis
chromothripsisSim = function(chromLength=1000,chrom="22",
			pLoss=0.2, # probability of losing a segment
			poisMean=50 # poisson mean for number of breakpoints
			)
	{
	# get breakpoints
	nBreaks = rpois(1,poisMean)
	# clustered breakpoints
	# oscillating probability
	#nPeaks = rpois(1,1)+1
	#prob = sin(seq(from=0,to=2*pi*nPeaks,length.out=chromLength))+1
	#perturbation = sample(2:chromLength,1)
	#prob = prob[c(perturbation:chromLength,1:perturbation-1)]
	# single chromothriptic region
	size = chromLength/10
	start = sample(chromLength-size,1)
	prob = rep(1,length=chromLength)
	prob[start:(start+size)] = 10
	breakpoints = sort(sample(chromLength,nBreaks,prob=prob,replace=FALSE))
	# make segments	
	segStarts = c(0,breakpoints)
	segEnds = c(breakpoints,chromLength)
	seg = data.frame(segStarts,segEnds,"+","-",stringsAsFactors=FALSE)
	# lose some segments
	loss = rbinom(nrow(seg),1,pLoss)	
	if(sum(loss)>0) seg = seg[-which(loss==1),]
	# orientation of segments for stitching
	orientation = rbinom(nrow(seg),1,0.5)
	seg[which(orientation==1),1:2] = seg[which(orientation==1),2:1]
	seg[which(orientation==1),3:4] = seg[which(orientation==1),4:3]
	# new order of segments
	newOrder = sample(nrow(seg),replace=FALSE)
	seg = seg[newOrder,]
	if(nrow(seg)>1)
		{
		# fusions
		out = cbind(chrom,seg[-nrow(seg),c(2,4)],chrom,seg[-1,c(1,3)])
		colnames(out) = c("chrom1","pos1","direction1","chrom2","pos2","direction2")
		return(out)
		} else {
		return(NA)
		}
	}



# perform chromothripsis simulation once
singleSimTest = function(chrom="22",chromLength=50818468,pLoss=0.2,poisMean=50)
	{
	bedpe = chromothripsisSim(chromLength,chrom,pLoss,poisMean)
	while(is.na(bedpe))
		{
		bedpe = chromothripsisSim(chromLength,chrom,pLoss,poisMean)
		}
	any(splitWindow(bedpe=bedpe,
		chrom=chrom,
		chromCol1=1,
		posCol1=2,
		chromCol2=4,
		posCol2=5,
		direction1col=3,
		direction2col=6,
		chromStart=0,
		chromEnd=chromLength))
	}

# perform chromothripsis simulation multiple times
multiSimTest = function(nSim=10,chrom="22",chromLength=50818468,pLoss=0.2,poisMean=50)
	{
	test = replicate(nSim,singleSimTest(chrom,chromLength,pLoss,poisMean))
	sum(test[which(!is.na(test)&test)])/nSim
	}

