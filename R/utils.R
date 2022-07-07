#RETURNS THE PDF OF THE GIVEN DISTRIBUTION
.getDistPDF <- function(pars, dist){
  if (dist=="normal")
    return(function(x) dnorm(x, mean=pars["mean"], sd=pars["sd"]))
  else if (dist=="lognormal")
    return(function(x) dlnorm(x, meanlog=pars["meanlog"], sdlog=pars["sdlog"]))
}

#RETURNS THE RANDOM NUMBER GENERATOR OF THE GIVEN DISTRIBUTION
.getDistRNG <- function(pars, dist){
  if (dist=="normal")
    return(function(x) rnorm(x, mean=pars["mean"], sd=pars["sd"]))
  else if (dist=="lognormal")
    return(function(x) rlnorm(x, meanlog=pars["meanlog"], sdlog=pars["sdlog"]))
}

#RETURNS THE DESCRIPTION OF A GIVEN DISTRIBUTION
.getDistDesc <- function(pars, dist) {
  if (dist=="normal")
    return(paste("normal(mean=", round(pars["mean"], 2), ", ", "sd=", round(pars["sd"], 2),  ")", sep=""))
  else if (dist=="lognormal")
    return(paste("lognormal(meanlog=", round(pars["meanlog"], 2), ", ", "sdlog=", round(pars["sdlog"], 2), ")", sep=""))
}

##RETRIEVES THE MEAN OF THE PRIOR FOR MOSQUITO BITING FREQUENCY
.getPriorMosqBitingFreqMean <- function(pars, dist) {
  if (dist=="normal")
    return(pars["mean"])
  else if (dist=="lognormal")
    return(exp(pars["meanlog"]+pars["sdlog"]^2/2))
}

#SMOOTHS A TIME SERIES WITH +- N STEPS AVERAGE PER POINT
.smoothUDSeries<- function(series, n){
  return(as.numeric(stats::filter(series, rep(1/n, n), sides=2, circular=TRUE)))
}

##CALCULATES MOSQ biting RATE COMPONENT DEPENDENT ON HUMIDITY IN TIME
.mvse_hum_effect_aV<- function(U,meanU){
  ef<- (U-meanU)/sqrt(1+(U-meanU)^2)
  return(ef)
  ##no need to trim the function for bio meaning
}

##CALCULATES MOSQ DEATH RATE COMPONENT DEPENDENT ON TEMPERATURE IN TIME
.mvse_temp_effect_muV <- function(temp){
  ef <- (0.8692-0.1599*temp+(0.01116*(temp^2))-0.0003408*(temp^3)+0.000003809*(temp^4))
  ef[which(1/ef>180 | ef < 0)]<- 1/180 # fix negative numbers and biologically impossible numbers to 1/180 (180 days)
  return(ef)
}

##CALCULATES MOSQ DEATH RATE COMPONENT DEPENDENT ON HUMIDITY IN TIME
.mvse_hum_effect_muV <- function(U, meanU){
  ef <- -(U-meanU)/sqrt(1+(U-meanU)^2) ##v1.0.3
  return(ef)
  ##no need to trim the function for bio meaning
}

##CALCULATES EPSVH (PROB TRANS VECTOR TO HUMAN) DEPENDENT ON TEMPERATURE IN TIME
.mvse_temp_epsVH <- function(temp){ ##testing
  ep <- 0.001044*temp*(temp-12.286)*((32.461-temp)^(0.5))
  ep[which(ep<0)] <- 0 # fix negative numbers as bio -> trim to zero
  ep[which(is.nan(ep))] <- 0
  return(ep)
}

##CALCULATES EPSHV (PROB TRANS HUMAN TO VECTOR) DEPENDENT ON TEMPERATURE IN TIME
.mvse_temp_epsHV <- function(temp){ ##testing
  ep <- 0.0729*temp-0.9037
  ep[temp<12.4] <- 0
  ep[temp>26.1] <- 1
  return(ep)
}

##CALCULATES MOSQ INCUBATION TO INFECTION RATE DEPENDENT ON TEMPERATURE IN TIME
.mvse_temp_effect_gammaV <- function(temp){
  Tk <- temp+273.15
  R <- 1.987 # cal/(mol*K)
  ef <- (24.0*(0.003359*(Tk/298.)*exp((15000./R)*(1/298.-1./Tk))/(1.+ exp((6.203*(10^21)/R)*(1./(-2.176*(10^30))-1./Tk)))))
  ef[which(ef<0)]<- 0 #fix negative numbers as bio -> trim to zero
  return(ef)
}

##CALCULATES MOSQ aquatic mortality rate COMPONENT DEPENDENT ON TEMP IN TIME
.mvse_temp_effect_muA <- function(temp){
  ef<- (2.13-0.3797*temp+0.02457*temp^2-0.0006778*temp^3+0.000006794*temp^4)
  return(ef)
}

##CALCULATES MOSQ aquatic dev rate COMPONENT DEPENDENT ON TEMP IN TIME
.mvse_temp_effect_epsA <- function(temp){
  ef<-  0.131-0.05723*temp+0.01164*temp^2-0.001341*temp^3+0.00008723*temp^4-0.000003017*temp^5+5.153*10^(-8)*temp^6-3.42*10^(-10)*temp^7
  return(ef)
}

##CALCULATES MOSQ oviposition rate COMPONENT DEPENDENT ON TEMP IN TIME
.mvse_temp_effect_theta <- function(temp){
  ef<-  -5.4+1.8*temp-0.2124*temp^2+0.01015*temp^3-0.0001515*temp^4
  return(ef)
}

##CALCULATES MOSQ egg success rate COMPONENT DEPENDENT ON TEMP IN TIME
.mvse_temp_effect_C<- function(temp){
  ef<-  (-184.8+27.94*temp-0.9254*temp^2+0.009226*temp^3)/100
  return(ef)
}

##CALCULATES MOSQ egg success rate effect COMPONENT DEPENDENT ON RAINFALL IN TIME
.mvse_rain_effect_C <- function(R, meanR){
  ef <- (R-meanR)/sqrt(1+(R-meanR)^2)
  return(ef)
}
