##################################
# Dynamic Conditional Correlation
##################################

#require(methods)

############
# corMat & update
############

cor.update <- function(x, innovations, lambda)
{
	x <- cov.update(x, innovations, lambda)
	x <- cov2cor(x)
	x <- 0.5 * (x + t(x))
	x
}

cov.update <- function(x, data, lambda)
{
	valid.data <- is.finite(data)
	
	if(sum(valid.data) > 0)
	{
		x[valid.data, valid.data] <- lambda * x[valid.data, valid.data] + 
				(1 - lambda) *  crossprod(t(data[valid.data]))
		x <- 0.5 * (x + t(x))		
	}
	x 
}


DCC <- function(data, lambda = 0.94, showMsg=TRUE, initial.matrix = NULL, aggregate.over.days = 1)
{		
	if(is.null(initial.matrix))
	{
		#initialize the correlation matrix, start off as identity
		cor.mat <- diag(1, NCOL(data))
		rownames(cor.mat) <- colnames(cor.mat) <- colnames(data)
	
	}
	if(!is.null(initial.matrix))
	{
		if(NROW(initial.matrix) != NCOL(initial.matrix)) stop("\n\"initial.matrix\" must be square")
		if(NROW(initial.matrix) != NCOL(data)) stop("\n\"initial.matrix\" doesn't have the same number of columns as the input \"data\"")
		
#		if(any(!rownames(initial.matrix)%in%colnames(data))) stop("\nrownames of the initial.matrix should match the colnames of the input \"data\"\n")
#		if(any(!colnames(initial.matrix)%in%colnames(data))) stop("\nrownames of the initial.matrix should match the colnames of the input \"data\"\n")
		
		rownames(initial.matrix) <- colnames(initial.matrix) <- colnames(data)
		cor.mat <- initial.matrix
	}

	if(aggregate.over.days > 1)
	{
		# use the filter method, much faster than rollapply
		data <- .aggregateData(data, aggregate.over.days, method="filter")
		#remove the NA's
		data <- data[-(1:(aggregate.over.days - 1)),]
	}
	
	
	# show a progress window
	if (showMsg){
		tstr = "Calculating DCC"
		pb <- winProgressBar(title=tstr,label=" ",min=1,max=NROW(data))
		#Sys.sleep(0.1)		
	}
	
	#increment over time, updating the correlation matrix estimate for one step ahead
#	pb <- winProgressBar(title=sprintf("Calculating correlation matrix with lambda = %.4f",lambda),label=" ",min=1,max=NROW(data))
	for(i in 1:NROW(data))
	{
		if (showMsg){
			infostr = sprintf("%d%% done",round(100*(i/NROW(data))))
			setWinProgressBar(pb,value=i,label=infostr,title=tstr)
		}
	
		inos <- data[i,]
		cor.mat <- cor.update(cor.mat, innovations = as.numeric(inos), lambda = lambda)
	}
	
	if (showMsg){
		Sys.sleep(0.1)
		close(pb)		
	}
#	close(pb)         
	
	cor.mat
	
}



.aggregateData <- function(data, aggregate.over.days = 5, method=c("filter","rollapply")){
	
	method = match.arg(method)
	
	if (method=="rollapply"){
		data <- rollapply(data,width = floor(aggregate.over.days),
				FUN = sum, align = "right", na.pad=TRUE) / sqrt(aggregate.over.days)		
	}
	if (method=="filter"){
		f <- rep(1,aggregate.over.days)
		data.in <- data
		data <- reclass(filter(data,filter=f,method="convolution",sides=1), data) / 
				sqrt(aggregate.over.days)		
		colnames(data) <- colnames(data.in)
	}
	return(data)
}



#.ewma.na.pad <- function(x, alpha)
#{
#	ewma(x, alpha = alpha)
#}
#
#setOldClass("xts")
#
#setClass("CorMat",representation("matrix", lambda ="numeric"))
#
#setGeneric("update",function(x, innovations) standardGeneric("update") )
#
#setMethod("update", signature(x = "CorMat", innovations = "numeric"), 
#		function(x, innovations)
#		{
#			lambda <- x@lambda
#			x <- cov.update(x, innovations, lambda)
#			x <- cov2cor(x)
#			x <- 0.5 * (x + t(x))
#			new("CorMat", x, lambda = lambda)
#		})
#
#setMethod("update", signature(x = "CorMat", innovations = "xts"), 
#		function(x, innovations)
#		{
#			x <- update(x, as.numeric(innovations))
#		})
#
#setMethod("update", signature(x = "CorMat", innovations = "zoo"), 
#		function(x, innovations)
#		{
#			x <- update(x, as.numeric(innovations))
#		})




