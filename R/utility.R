# ======================================================================================
#			UTILITY FUNCTIONS
# ======================================================================================

# read in a file
readFile = function(file,head)
	{
	tmp = gsub("[.]gz","",file)
	ending = rev(strsplit(tmp,split="[.]")[[1]])[1]
	if(ending=="csv") return(read.csv(file,head=head,as.is=TRUE))
	if(ending%in%c("txt","tsv")) return(read.table(file,sep="\t",head=head,as.is=TRUE))
	return(read.table(file,sep="\t",head=head,as.is=TRUE))	
	}


# load cytoBand file
loadCytoBand = function(file=NULL)
	{
	if(is.null(file))
		{
		tmpEnv = new.env()
		data(list="cytoBand", package='CNomplexity',envir=tmpEnv)
		return(tmpEnv[["cytoBand"]])
		} else {
		return(readFile(file))
		}
	}
