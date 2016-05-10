

ac.cov.update <- function(x, ino, ino.lag, lambda)
{
	valid.ino <- is.finite(ino)
	valid.ino.lag <- is.finite(ino.lag)
	
	valid.data <- valid.ino & valid.ino.lag
	
	if(sum(valid.data) > 0)
	{
		x[valid.data, valid.data] <- lambda * x[valid.data, valid.data] + 
				(1 - lambda) *  ino[valid.data] %*% t(ino.lag[valid.data])#crossprod(t(ino[valid.data]))
		x <- 0.5 * (x + t(x))		
	}
	x 
}


DCAC <- function(data, lambda = 0.94, showMsg=TRUE, initial.ac.matrix = NULL, 
		initial.cov.matrix = NULL, 
		aggregate.over.days = 1, get.auto.cor = FALSE)
{	
	#require data to be an xts
	if(!any(class(data)%in%"xts")) stop("\nin DCAC, data needs to be an \"xts\"\n")

	if(is.null(initial.ac.matrix))
	{
		#initialize the correlation matrix, start off as identity
		acov.mat <- matrix(0, ncol = NCOL(data), nrow = NCOL(data))
		rownames(acov.mat) <- colnames(acov.mat) <- colnames(data)
	}
	
	if(is.null(initial.cov.matrix))
	{
		#initialize the covariance matrix
		cov.mat <-  matrix(0, ncol = NCOL(data), nrow = NCOL(data))
		rownames(cov.mat) <- colnames(cov.mat) <- colnames(data)
		
	}
	
  #check initial matrices are of the correct size   
	if(!is.null(initial.ac.matrix))
	{
		if(NROW(initial.ac.matrix) != NCOL(initial.ac.matrix)) stop("\n\"initial.matrix\" must be square")
		if(NROW(initial.ac.matrix) != NCOL(data)) stop("\n\"initial.matrix\" doesn't have the same number of columns as the input \"data\"")
		
		rownames(initial.ac.matrix) <- colnames(initial.ac.matrix) <- colnames(data)
		acov.mat <- initial.ac.matrix
	}
	
	if(!is.null(initial.cov.matrix))
	{
		if(NROW(initial.cov.matrix) != NCOL(initial.cov.matrix)) stop("\n\"initial.matrix\" must be square")
		if(NROW(initial.cov.matrix) != NCOL(data)) stop("\n\"initial.matrix\" doesn't have the same number of columns as the input \"data\"")
		
		rownames(initial.cov.matrix) <- colnames(initial.cov.matrix) <- colnames(data)
		cov.mat <- initial.cov.matrix
	}
	
	if(aggregate.over.days > 1)
	{
		# use the filter method, much faster than rollapply
		data <- .aggregateData(data, aggregate.over.days, method="filter")
		#remove the NA's
		data <- data[-(1:(aggregate.over.days - 1)),]
	}
	
	#lag data
	data.lag <- lag(data, k = aggregate.over.days, na.pad = FALSE)
	
	#make sure the dates match
	data <- data[index(data.lag)]
	
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
		inos.lag <- data.lag[i,]
		acov.mat <- ac.cov.update(acov.mat, ino = as.numeric(inos), 
					ino.lag = as.numeric(inos.lag), lambda = lambda)
	
    if(get.auto.cor)
    {
      cov.mat <- ac.cov.update(cov.mat, ino = as.numeric(inos), 
                               ino.lag = as.numeric(inos), lambda = lambda)
    }
	}
	
	if (showMsg){
		Sys.sleep(0.1)
		close(pb)		
	}
  output <- acov.mat	
	if(get.auto.cor) 
  {
    output <- sweep(acov.mat, 2, sqrt(diag(cov.mat)), "/")
    output <- sweep(output, 1, sqrt(diag(cov.mat)), "/")
	}
	
  output
}