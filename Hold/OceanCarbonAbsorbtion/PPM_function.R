austinPetzold_1986 <- function(lambda, K490){
  
  # -- ORIGINAL DATA FROM AUSTIN AND PETZOLD (1986) -- #
  wave = seq(350,700,10)
  M = c(2.1442, 2.0504, 1.9610, 1.8772, 1.8009, 1.7383, 1.7591, 1.6974, 1.6108, 1.5169, 1.4158, 1.3077,
        1.1982, 1.0955, 1.0000, 0.9118, 0.8310, 0.7578, 0.6924, 0.6350, 0.5860, 0.5457, 0.5146, 0.4935,
        0.4840, 0.4903, 0.5090, 0.5380, 0.6231, 0.7001, 0.7300, 0.7301, 0.7008, 0.6245, 0.4901, 0.2891)
  
  Kdw <- c(0.0510,0.0405,0.0331,0.0278,0.0242,0.0217,0.0200,0.0189,0.0182,0.0178,0.0176,0.0176,
           0.0179,0.0193,0.0224,0.0280,0.0369,0.0498,0.0526,0.0577,0.0640,0.0723,0.0842,0.1065,
           0.1578,0.2409,0.2892,0.3124,0.3296,0.3290,0.3559,0.4105,0.4278,0.4521,0.5116,0.6514)
  
  # -- INTERPOLATE TO WAVELENGTH OF INTEREST-- #
  M_l   <- approx(x=wave, y=M, xout=lambda)$y
  Kdw_l <- approx(x=wave, y=Kdw, xout=lambda)$y
  
  # -- GET REFERENCE WAVELENGTH (=490 FOR NOW) AND APPLY MODEL-- #
  ref= which(wave==490)
  
  Kd = (M_l/M[ref]) * (K490 - Kdw[ref]) + Kdw_l
  
  return(Kd)
}


# this needs to be changed accroding to the MLD grid size and also the resoltuion of
# the sat data (e.g. 4km or 9km)

MLD_stretch_larger <- function(MLD){
  
  MLD.new <- array(NA,c(4320,2160))
  for (i in 1:4320){
    for (j in 1:2160){
      MLD.new[((i-1)*2+1):((i-1)*2+2),((j-1)*2+1):((j-1)*2+2)] <- rep(MLD[i,j],2)  
    }
  }
  return(MLD.new)
  
  
}

read_hdf_files <- function(x, files.order){

# load in raster/hdf file


file <- raster(x[files.order[20]])

# set the coords limits and system

 extent(file) <- extent(-180, 180, -90, 90)
 crs(file) <- "+proj=longlat +proj=eck4"

# file[file<0]   <- NA
# plot(file)

# convert to a dataframe and then a matrix
## FOR DEVELOPMENT PURPOSES coarsen grid to 1 of a degree from 0.1
# file$x <-  trunc(file$x*10^1)/10^1
# file$y <-trunc(file$y*10^0)/10^0

file <- data.frame(rasterToPoints(file))
colnames(file)[3] <- "z"

# file$x <-  round_any(file$x, 10)
# file$y <-  round_any(file$y, 10)
# file$z[file$z<0]   <- NaN  
# file <- aggregate(z ~ x+y, file, FUN = mean, na.action = na.pass)

file <- acast(file, rev(y)~x, value.var="z")

return(file)                                   

}

daylength <- function(doy,lat){
  
  # calculates daylength.  required inputs are day of year (yd)
  # and latitude (lat)
  
  dec <- 23.5*cos((2*pi*(doy - 172))/365);       #declination of sun
  tmp <- -tan(lat*pi/180)*tan(dec*pi/180)
  tmp[tmp>1]  <- 1                              #Check for daylengths less than 24 hours
  tmp[tmp<(-1)] <- -1                           #Check for daylengths greater than 24 hours
  dl  <- acos(tmp)*180/pi/15*2                  #sunrise hour angle / 15deg/hr x 2 (sr +ss)
  return(dl)                                    #Hours
}


PPM <- function(bbp, chl, irr, k490, MLD, doy, month){
  
  #Remove negative fill values
  bbp[bbp<0]   <- NA  
  chl[chl<0]   <- NA
  irr[irr<0]   <- NA
  k490[k490<0] <- NA
  MLD[MLD<0] <- NA
  
  
  #Stretch MLD and ZNOx
  
  # this can be used if the dimensions of the MLD data grid do not match. Check to see whether
  # the MLD data are compatible with 4km or 9km sat data
  
  
  #Isolate Ocean Cells with data
  ocean <- which(is.finite(bbp) & is.finite(k490) & is.finite(chl) 
                 & is.finite(irr) & is.finite(MLD), 
                 arr.ind = T)
  
  bbp   <- bbp[ocean]  
  chl   <- chl[ocean]
  irr   <- irr[ocean]
  k490  <- k490[ocean]
  MLD   <- MLD[ocean]
  

  #############################################
  # 2) Declare Spectral Values
  lambda      <- c(400, 412, 443, 490, 510, 555, 625, 670, 700)
  parFraction <- c( 0.0029, 0.0032, 0.0035, 0.0037, 0.0037, 0.0036, 0.0032, 0.0030, 0.0024)
  X           <- c(.11748, .122858, .107212, .07242, .05943, .03996, .04000, .05150, .03000)
  e           <- c(.64358, .653270, .673358, .68955, .68567, .64204, .64700, .69500, .60000) 
  Kw          <- c(.01042, .007932, .009480, .01660, .03385, .06053, .28400, .43946, .62438)
  
  #############################################
  # 3) Initalize Values
  
  y0   = 0.0003
  z    = c(0:200)
  r    = 0.1
  uMax = 2.0
  
  Klambda <- array(NA, c(length(chl), 9))
  E       <- array(NA, c(length(chl), 9))
  Kdif    <- array(NA, c(length(chl), 9))
  PAR     <- array(0, c(length(chl)))
  
  #Create Daylength matrix
  
  
  dl      <- array(0, c(4320, 2160))
  lat     <- seq(-89.958336, by=0.083333336, length.out=2160)
  dlin    <- rev(daylength(doy[1],lat)) 
  for (i in 1:4320) {dl[i,] = rep(dlin[i], 2160)}
  dl      <- dl[ocean[,2]]
  
  # Need to increase the deimensions of the dl box from 2160 x 4320 to 4320 x 8640 if you 
  # are using 4km resolution sat data
  
  #Compute median PAR in MLD
  
  for (i in 1:9) {
    Klambda[,i]  <- austinPetzold_1986(lambda[i],k490)
    E[,i]        <-  irr * parFraction[i] * 0.975 * exp(-Klambda[,i] * MLD / 2.0)}
  for (i in 1:8) {
    PAR <- PAR + (lambda[i+1] - lambda[i]) * (E[,i] + E[,i+1]) / 2}
  
  # PAR     <- PAR / dl  # Changed from above
  irr     <- irr / dl 
  PAR_kd <- irr^0.45/k490
  
  # PPM model coefficients
  a = -16.8
  b= 1.55
  c= 47.027
  d = 0.0125
  
  # Photoacclimation terms 
    dm <- 19 *exp(0.038*PAR_kd)
  sm <- ((1+exp(-0.15*irr))/(1+exp(-3*PAR)))
  sm[sm < 1] = 1
  dm[dm > 1000] = 1000
  theta = (dm*sm)  
  
  # use either bbp or model C
    carbon = 13000 * (bbp - 0.00035)
  # carbon = (dm*sm) * chl 
  
  # Light limited function for depth 
  IgFunc <- (1- exp(-5 * irr))
  
  # Light limitation function  
  lightmu <-  (I(1 / dm * a) + b) 
  
  # Nutrient limitation function  
  nutmu <-  (I(1 / (dm*sm) * c) + d) 
  
  # Estimate phytoplankton growth 
  mu <- (lightmu*nutmu) *  IgFunc
  mu[mu > uMax] = uMax
  
  for (i in 1:9) {  #Need to find out E just above MLD for later on
    Kdif[,i] <- Klambda[,i] - (Kw[i] +  X[i] * chl^e[i])
    E[,i]    <-  irr * parFraction[i] * 0.975 * exp(-Klambda[,i] * floor(MLD))}  
  
  #############################################
  # 3) Compute NPP through depth
  NPP     <- array(NA, c(length(irr), 201))
  
  #Parameters dependent on level z-1
  mu.z     <- mu
  chl.z    <- chl
  E.z      <- E
  
  #save(mu, file = paste(getwd(), "/Model_output/mu/", month,".RData", sep=""))  
  #save(carbon, file = paste(getwd(), "/Model_output/carbon/", month,".RData", sep=""))  
  #save(MLD, file = paste(getwd(), "/Model_output/MLD/", month,".RData", sep=""))  
  
  #Fill in NPP where data is present. Will overwrite NPP for  data beneath the MLD
  for (k in 1:201){NPP[,k] <- mu * carbon} 
  
  #Start at shallowest MLD  
  surface <- ceiling(min(MLD))+1
  
  # rm(Klambda);rm(Kdif);
  # rm(E); rm(E.z)
  # rm(bbp);rm(chl);rm(dl);rm(theta)
    # gc()
  
  
  for (k in 1:201){
   
    # k= 10 
    # deep <- which(MLD<=z[k], arr.ind = T)
    
    irr.z        =  irr * exp(-k490*z[k])
    
    IgFunc      = (1 - exp(-5 * irr.z))
    
    lightmu[k] <-  (I(1 / dm[k] * a) + b) 
    nutmu[k] <-  (I(1 / (dm[k]*sm[k]) * c) + d) 
    
    mu[k] <- (lightmu[k]*nutmu[k]) *  IgFunc[k]
    
    mu[k][ mu[k]>uMax] <- uMax
    
    ind <- which(mu.z <  r & MLD<=z[k], arr.ind=T)

    # chl[k]    = chl_c[k] * carbon[k]
    NPP[,k]  = round(mu[k] * carbon,3)
    
    #Parameters dependent on level z-1
    mu.z[k]   <- mu[k]
    # chl.z[k]  <- chl[k]
    # NPP.z[,k]  = round(mu.z[k] * carbon,3)
    # for (l in 1:9){E.z[k,l] <-  E[k,l]}
    
    print(k)
    
  }
  
  
  # Need to parse between above and below the mixed layer 

  #Remove some data if necessary 
  #rm(E.z); rm(deltaZ); rm(IgFunc);
  #rm(chl.z); rm(mu.z); rm(nutTempFunc);
  
  NPP.surface <- array(NA,c(2160,4320))   
  NPP.all     <- array(NA,c(2160,4320)) 
  NPP.deep    <- array(NA,c(2160,4320))
  
  #For integration, truncate MLD to 200 m
  MLD[MLD>200] = 200
  
  for (i in 1:length(irr)){
    #Integrate
    fn <- splinefun(z,NPP[i,])       
    if(class(try(integrate(fn,0,MLD[i]), silent = TRUE))!="try-error"){
      NPP.surface[ocean[i,1],ocean[i,2]] <- round(integrate(fn,0,MLD[i])$value,digits=2)}
    if(class(try(integrate(fn,0,200), silent = TRUE))!="try-error"){
      NPP.all[ocean[i,1],ocean[i,2]]     <-  round(integrate(fn,0,200)$value,digits=2)}
  }
  
  # NPP.all   <- NPP.all[ocean]  
  # NPP.surface   <- NPP.surface[ocean]  
  NPP.deep <- NPP.all - NPP.surface
  
  #Save Rdata
  save(NPP.surface, file = paste(getwd(), "/Model_output/NPP Surface/", month,".RData", sep="")) 
  save(NPP.deep, file = paste(getwd(), "/Model_output/NPP Deep/", month,".RData", sep=""))  
  save(NPP.all, file = paste(getwd(), "/Model_output/NPP/", month,".RData", sep=""))  
  
}




