


## make (dc) auto correlation reference frames using endpoints
dcacRefFrames <- function(rar.norm = NULL,cor.decay = 500, on = "quarters", aggregate.over.days = 2,
                         initial.cor.mat = NULL, prog.bar = FALSE, get.auto.cor = FALSE)
{
  if(is.null(rar.norm)) stop("\n\"rar.norm\" is NULL, please provide")
  
  #set the initial correlation matrix, if it's not provided
  if(is.null(initial.cor.mat)) initial.cor.mat <- diag(NCOL(rar.norm))
  
  if(NCOL(rar.norm) != NCOL(initial.cor.mat)) stop("\nnumber of columns in initial.cor.mat does not match number of columns in rar.norm\n")
  
  ###### using "normal" risk adjusted returns to get the correlation out
  eps <- endpoints(rar.norm, on = on)
  eps <- eps[-length(eps)]
  
  #bit of hack here.. removing eps
  if(eps[2] - eps[1] < aggregate.over.days) eps <- eps[-1]
  
  if(prog.bar) 
  {
    tstr = "Calculating DCC for reference frames"
    pb <- winProgressBar(title = tstr, label = " ", min = 1, 
                         max = length(eps) - 1)
  }
  
  #generate the reference frames
  ref.cor.mat <- lapply(1:(length(eps) - 1), function(i)
  {
    if(prog.bar) 
    {
      infostr = sprintf("%d%% done", round(100 * (i/(length(eps) - 1))))
      setWinProgressBar(pb, value = i, label = infostr, 
                        title = tstr)
    }
    #dates (from endpoints) specifying interval
    date.1 <- eps[i] + 1
    date.2 <- eps[i + 1]
    
    #specify the interval of time
    use.interval <- seq(date.1, date.2, 1)
    
    #evolve the correlation matrix over the intervals
    output <- DCAC(rar.norm[use.interval,], lambda = d.to.lambda(cor.decay), 
                  initial.cov.matrix = initial.cor.mat,
                  aggregate.over.days = aggregate.over.days, showMsg = FALSE,
                   get.auto.cor = get.auto.cor)
    
    #update the initial correlation matrix, for the next interval
    initial.cor.mat <<- output
    output
  })
  
  if (prog.bar) close(pb)
  
  #get the dates which each reference frame represents
  cor.mat.dates <- index(rar.norm[eps[2:length(eps)]])	
  
  #give names (reference dates) to the list (reference frames)
  names(ref.cor.mat) <- as.character(cor.mat.dates)
  
  ref.cor.mat
}

getDcacMat <- function(ref.cor.mat = NULL, date = NULL, rar.norm = NULL, 
                      usable.markets = NULL, 
                      cor.decay = 500, aggregate.over.days = 2,
                       get.auto.cor = FALSE)
{
  if(is.null(ref.cor.mat)) stop("\n\"ref.cor.mat\" not provided\n")
  if(is.null(date)) stop("\n\"date\" not provided\n")
  if(is.null(rar.norm)) stop("\n\"rar.norm\" not provided\n")
  
  if(is.null(usable.markets)) usable.markets <- colnames(rar.norm)
  
  #set date to be of class POSIX
  date <- as.POSIXct(date)
  
  #get the dates of the reference correlation matrices
  cor.mat.dates <- names(ref.cor.mat)
  
  #get the most recent reference correlation matrix to use
  use.ref.cor.mat <- max(which(date >= as.POSIXct(cor.mat.dates, tz = "UTC")))
  
  #think about this step here, currently uses the previous reference frame if
  #aggreagte.over.days is > 1, which is to avoid the case of not having enough data 
  #when trying to aggregate over multiple days
  #REMOVED
  if(aggregate.over.days > 1) use.ref.cor.mat <- use.ref.cor.mat - 1
  
  cm.date <- as.POSIXct(cor.mat.dates[use.ref.cor.mat], tz = "UTC")
  
  #the date to start updating the reference correlation matrix				
  start.date <- index(rar.norm)[which(index(rar.norm)%in%cm.date) + 1]
  #	start.date <- index(rar.norm)[which(index(rar.norm)%in%cm.date) + 2 - aggregate.over.days]
  
  
  #the interval to update the correlation matrix over
  date.interval <- paste(start.date, date, sep = "::")
  initial.cor.mat <- ref.cor.mat[[use.ref.cor.mat]]
  
  #if the var.date lands on the "reference date", just use that correlation
  #matrix
  if(date%in%cm.date)
  {
    acor.mat <- initial.cor.mat[usable.markets, usable.markets]
  }
  if(!date%in%cm.date)
  {				
    acor.mat <- DCAC(rar.norm[date.interval,usable.markets], showMsg = FALSE,
                   lambda = d.to.lambda(cor.decay),
                   initial.cov.matrix = initial.cor.mat[usable.markets, usable.markets],
                   aggregate.over.days = aggregate.over.days,
                   get.auto.cor = get.auto.cor)
  }
  
  acor.mat
  
}
