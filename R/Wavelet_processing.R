#'@title Internal wavelet processing
#'@description  Internal function, perform a wavelet transform and thresholding
#'@param Time01 rescalled measurment positions (ordered)
#'@param y the signal
#'@param lev_res the maximum level of resolution needed

wavproc <- function(Time01,y,lev_res)
{
  #Kovac and Silvermann 2000
  mygrid <- wavethresh::makegrid(t=Time01,y=y)
  LDIRWD <- irregwd(mygrid,filter.number=1)
  class(LDIRWD) <- "wd"
 
  myThres <- threshold(LDIRWD,policy = "universal",type="soft",dev = madmad,levels = 1:(LDIRWD$nlevels-1))


 

  res <- c()
  for( i in 0: lev_res)
  {
    res <- c(res, accessD( LDIRWD,lev = i) )
  }

  res
}
