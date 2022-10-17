##########################################################################################################
#' Class representing MVSE model
#'
#' A \code{mvsemodel} object contains both climate data and prior distributions
#' for  human/vector/virus biological parameters. The \code{sampling}
#' method defined in this class may be used to draw samples of the ecological scaling coefficients
#' (rho, eta, alpha) from this model.
#'
#' Instances of \code{mvsemodel} are created by calling the function \code{\link{mvse_model}}.
#'
#' @slot model_name The model name, an object of type \code{character}.
#' @slot model_category The model category, an object of type \code{character}.
#' @slot climate_data The climate data, an object of type \code{dataframe}.
#' @slot priors The prior distributions of the human/vector/virus biological parameters,
#' an object of type \code{list}.
#'
#' @section Methods:
#' \describe{\item{\code{show}}{Print the default summary for the model.}}
#' \describe{\item{\code{print}}{Print all information contained in the model.}}
#'
#' @name mvsemodel-class
#' @rdname mvsemodel
#' @include utils.R
#' @import dplyr tidyr ggplot2 RColorBrewer
#' @importFrom bayesplot bayesplot_grid
#' @importFrom methods new
#' @importFrom stats dlnorm dnorm rlnorm rnorm
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom rlang .data 
setClass("mvsemodel",
         slots = c(
           model_name = "character",
           model_category = "character",
           climate_data = "data.frame",
           priors = "list"
         ),
         prototype = list(
           model_name = character(0),
           model_category = character(0),
           climate_data = data.frame(),
           priors = list()
         )
)

##########################################################################################################
#' Construct a MVSE model
#'
#' Construct an instance of S4 class \code{mvsemodel}. A \code{mvsemodel} object can be used to draw samples
#' of the ecological scaling coefficients (i.e. rho, eta, alpha) along with the transmission potential measure (i.e. index P) 
#' from the defined model.
#'
#' @param model_name A character string naming the model; defaults to \code{"anon_model"}.
#' @param model_category A character string naming the model category. If \code{"denv_aegypti"} is specified,
#' the informed probability distributions of human/vector/virus biological parameters will automatically be set.
#' Otherwise, it defaults to \code{user-defined} and the user must specify the informed probability distributions for these parameters.
#' @param climate_data A data frame for the climate data. The data frame should have four named columns; \code{date} for dates in the form
#' \code{"\%Y-\%m-\%d"}, \code{T} for temperature in degrees Celsius and \code{H} for relative humidity (\%).
#' @param priors A named list of lists specifying the priors for the host/vector/virus biological parameters. 
#' (see \strong{Details} section for further information).
#' @param warning An optional logical specifying whether warning messages should be displayed. 
#'
#'
#' @details The \code{priors} argument requires data in a specific format. The list must contain five lists \code{mosq_life_exp}, \code{mosq_inc_per},
#' \code{mosq_biting_freq}, \code{human_life_exp} and \code{human_inc_per} which specify the distributions of mosquito life expectancy (days),
#' mosquito incubation period (days), mosquito biting frequency (bites/female/day), human life expectancy (years), and
#' human incubation period (days) respectively. The list for each of the five distributions
#' are specified by \code{pars}, a named numeric vector specifying the parameters of the distribution, and \code{dist}, a string for the name of the distribution.
#' For example, a normally distributed mosquito incubation period with a mean of 5 days and a standard deviation 3 days can be
#' specified as \code{mosq_inc_per=list(dist="normal", pars=c(mean=5, sd=3))}. Currently only the normal distribution with parameters \code{mean}
#' and \code{sd} and the lognormal distribution with parameters \code{meanlog} and \code{sdlog} are available for use. 
#'
#' @return An instance of S4 class \code{\linkS4class{mvsemodel}} that can be passed to \code{\link{sampling}}.
#' 
#' @examples 
#' # Daily temperature, relative humidity and precipitation data from 1981 in 2019 in the 
#' # city of Feira de Santana in Bahia, Brazil.
#' data(climateFSA)
#' 
#' # User-defined model
#' priors <- list()
#' ### mosquito life expectancy (days)
#' priors$mosq_life_exp <- list(dist="normal", pars=c(mean=12, sd=2))
#' ### mosquito-virus incubation period (days)
#' priors$mosq_inc_per <- list(dist="normal", pars=c(mean=7, sd=2))
#' ### mosquito biting frequency (bites/female/day)
#' priors$mosq_biting_freq <- list(dist="normal", pars=c(mean=0.25, sd=0.01))
#' ### human life expectancy (years)
#' priors$human_life_exp <- list(dist="normal", pars=c(mean=71.1, sd=2))
#' ### human-virus incubation period (days)
#' priors$human_inc_per <- list(dist="normal", pars=c(mean=5.8, sd=1))
#' ### human-virus infectious period (days)
#' priors$human_inf_per <- list(dist="normal", pars=c(mean=5.9, sd=1))
#' user_model <- mvse_model(model_name="my_mvse_model", climate_data=climateFSA, priors=priors)
#' 
#' # Pre-defined model with priors for dengue virus transmission 
#' # by Aedes aegypti mosquitoes in human hosts. 
#' aegypti_model <- mvse_model(model_name="my_aegypti_mvse_model", 
#' model_category="denv_aegypti", climate_data=climateFSA)
#' 
#' @export
mvse_model <- function(model_name="anon_model", model_category="user-defined", climate_data=data.frame(), priors=list(), warning=TRUE) {
  if (model_category=="denv_aegypti") {
    if (warning) print("Model category is 'denv_aegypti' so user-defined prior distributions will be ignored")
    priors <- .get_aegypti_priors()
  }
  new("mvsemodel", model_name=model_name, model_category=model_category, climate_data=climate_data, priors=priors)
}

##########################################################################################################
#' Print the default summary for the model
#'
#' @param object An object of class \code{\linkS4class{mvsemodel}}.
#'
#' @return NULL
#' @export
setMethod("show", "mvsemodel",
          function(object) {
            cat("Class:", class(object), "\n")
            cat("Model name:", object@model_name, "\n")
            cat("Model category:", object@model_category, "\n")
            nrowShow <- min(10, nrow(object@climate_data))
            trans_climate_data <- object@climate_data[, c("date", "T", "H")]
            cat("Climate data (limited to first 10 rows): \n")
            print(trans_climate_data[1:nrowShow, ], quote=FALSE)
            cat("Priors: \n")
            cat("\t Mosquito life expectancy (days) \t\t:", .getDistDesc(pars=object@priors$mosq_life_exp$pars, dist=object@priors$mosq_life_exp$dist), "\n")
            if (object@model_category=="denv_aegypti") {
              cat("\t Mosquito incubation period (days) \t\t \n")
              for (ii in seq_along(object@priors$mosq_inc_per$temps)) {
                temp <- object@priors$mosq_inc_per$temps[ii]
                if (ii==1) temp <- paste0("  <", temp, "\u00B0C", sep="")
                else if (ii==length(object@priors$mosq_inc_per$temps)) temp <- paste("  >", object@priors$mosq_inc_per$temps[ii-1], "\u00B0C", sep="")
                else temp <- paste(object@priors$mosq_inc_per$temps[ii-1], "-", temp, "\u00B0C", sep="")
                description <- .getDistDesc(pars=object@priors$mosq_inc_per$dists[[ii]]$pars, dist=object@priors$mosq_inc_per$dists[[ii]]$dist)
                cat("\t\t\t\t\t", temp, "\t:", description, "\n")
              }
            }
            else {
              cat("\t Mosquito incubation period (days) \t\t:", .getDistDesc(pars=object@priors$mosq_inc_per$pars, dist=object@priors$mosq_inc_per$dist),  "\n")
            }
            cat("\t Mosquito biting frequency (bites/female/day) \t:", .getDistDesc(pars=object@priors$mosq_biting_freq$pars, dist=object@priors$mosq_biting_freq$dist),  "\n")
            cat("\t Human life expectancy (years) \t\t\t:", .getDistDesc(pars=object@priors$human_life_exp$pars, dist=object@priors$human_life_exp$dist),  "\n")
            cat("\t Human incubation period (days) \t\t:", .getDistDesc(pars=object@priors$human_inc_per$pars, dist=object@priors$human_inc_per$dist),  "\n")
            cat("\t Human infectious period (days) \t\t:", .getDistDesc(pars=object@priors$human_inf_per$pars, dist=object@priors$human_inf_per$dist),  "\n")
          }
)

##########################################################################################################
#' Print the default summary of model (with all information). 
#'
#' @param x An object of class \code{\linkS4class{mvsemodel}}.
#' @param ... Additional arguments passed to the \code{print} method for \code{mvsemodel} objects.
#' 
#' @return NULL
#' @export
setMethod("print", "mvsemodel",
          function(x, ...) {
            cat("Class:", class(x), "\n")
            cat("Model name:", x@model_name, "\n")
            cat("Model category:", x@model_category, "\n")
            climate_data <- x@climate_data[, c("date", "T", "H")]
            cat("Climate data: \n")
            print(climate_data, quote=FALSE)
            cat("Priors: \n")
            cat("\t Mosquito life expectancy (days) \t\t", .getDistDesc(pars=x@priors$mosq_life_exp$pars, dist=x@priors$mosq_life_exp$dist), "\n")
            if (x@model_category=="denv_aegypti") {
              cat("\t Mosquito incubation period (days) \t\t \n")
              for (ii in seq_along(x@priors$mosq_inc_per$temps)) {
                temp <- x@priors$mosq_inc_per$temps[ii]
                if (ii==1) temp <- paste0("  <", temp, "\u00B0C", sep="")
                else if (ii==length(x@priors$mosq_inc_per$temps)) temp <- paste("  >", x@priors$mosq_inc_per$temps[ii-1], "\u00B0C", sep="")
                else temp <- paste(x@priors$mosq_inc_per$temps[ii-1], "-", temp, "\u00B0C", sep="")
                description <- .getDistDesc(pars=x@priors$mosq_inc_per$dists[[ii]]$pars, dist=x@priors$mosq_inc_per$dists[[ii]]$dist)
                cat("\t\t\t\t\t", temp, "\t:", description, "\n")
              }
            }
            else {
              cat("\t Mosquito incubation period (days) \t\t:", .getDistDesc(pars=x@priors$mosq_inc_per$pars, dist=x@priors$mosq_inc_per$dist),  "\n")
            }
            cat("\t Mosquito biting frequency (bites/female/day) \t:", .getDistDesc(pars=x@priors$mosq_biting_freq$pars, dist=x@priors$mosq_biting_freq$dist),  "\n")
            cat("\t Human life expectancy (years) \t\t\t:", .getDistDesc(pars=x@priors$human_life_exp$pars, dist=x@priors$human_life_exp$dist),  "\n")
            cat("\t Human incubation period (days) \t\t:", .getDistDesc(pars=x@priors$human_inc_per$pars, dist=x@priors$human_inc_per$dist),  "\n")
            cat("\t Human infectious period (days) \t\t:", .getDistDesc(pars=x@priors$human_inf_per$pars, dist=x@priors$human_inf_per$dist),  "\n")
          }
)

##########################################################################################################
#' Set (or replace) the probability distribution of a human/vector/virus biological parameter.
#'
#' @param object An object of class \code{\linkS4class{mvsemodel}}.
#' @param name A character string specifying the parameter whose distribution is to be set (or replaced).
#' Valid strings include \code{"mosq_inc_per"}, \code{"mosq_life_exp"},
#' \code{"mosq_biting_freq"}, \code{"human_life_exp"}, \code{"human_inc_per"} and
#' \code{"human_inf_per"}.
#' @param dist A character string specifying the distribution name. Currently only normal and log-normal distributions
#' are supported.
#' @param pars A numeric vector specifying the parameters of the distribution. These must be consistent
#' with the indicated distribution. \code{mean} and \code{sd} for the normal distribution and \code{meanlog} and \code{sdlog}
#' for the log-normal distribution.
#'
#' @details The following parameters can be modified: mosquito life expectancy in days (\code{mosq_life_exp}),
#' mosquito incubation period in days (\code{mosq_inc_per}),
#' mosquito biting frequency in bites/female/day (\code{mosq_biting_freq}), human life expectancy in years (\code{human_life_exp}),
#' and human incubation period in days (\code{human_inc_per}) respectively. The in-build probability distributions for
#' the `denv_aegypti` model category cannot be modified.
#'
#' @return An instance of S4 class \code{\linkS4class{mvsemodel}} with the probability distribution updated for
#' the specified human/vector/virus biological parameter.
#'
#' @examples
#' # obtain climate data
#' data(climateFSA)
#'
#' # define mvse model without specifying priors
#' user_model <- mvse_model(model_name="my_mvse_model", climate_data=climateFSA, priors=list())
#'
#' # set priors of mvse model
#' # (e.g. Mosquito life expectancy to normal distribution with 
#' # mean 12 days and standard deviation 2 days)
#' user_model <- set_prior(object=user_model, name="mosq_life_exp", 
#' dist="normal", pars=c(mean=12, sd=2))
#'
#' # replace previously set prior
#' # (e.g. Mosquito life expectancy to normal distribution 
#' # with mean 15 days and standard deviation 4 days)
#' user_model <- set_prior(object=user_model, name="mosq_life_exp", 
#' dist="normal", pars=c(mean=15, sd=4))
#'
#' @usage NULL
#' @name set_prior
#' @export
setGeneric("set_prior",
           def=function(object, name, dist, pars)
           {standardGeneric("set_prior")})

#' @rdname set_prior
setMethod("set_prior", "mvsemodel",
          function(object, name, dist, pars){
            if (object@model_category=="denv_aegypti")
              stop("Model category is 'denv_aegypti': cannot change priors")
            if (name=="mosq_inc_per") {
              object@priors$mosq_inc_per <- list(dist=dist, pars=pars)
            } else if (name=="mosq_life_exp") {
              object@priors$mosq_life_exp <- list(dist=dist, pars=pars)
            } else if (name=="mosq_biting_freq") {
              object@priors$mosq_biting_freq <- list(dist=dist, pars=pars)
            } else if (name=="human_life_exp") {
              object@priors$human_life_exp <- list(dist=dist, pars=pars)
            } else if (name=="human_inc_per") {
              object@priors$human_inc_per <- list(dist=dist, pars=pars)
            } else if (name=="human_inf_per") {
              object@priors$human_inf_per <- list(dist=dist, pars=pars)
            } else {
              stop(paste0(name, " is not a valid parameter!"))
            }
            return(object)
          }
)

##########################################################################################################
#' Draw samples from the model defined by class \code{mvsemodel}
#'
#' An MCMC routine (Metropolis-Hastings algorithm) is run to draw samples of factors rho,
#' eta and alpha from the model defined by class \code{mvsemodel}. The resulting posterior distributions
#' of these factors are then used to draw samples of the time series of Index P.
#'
#' @param object An object of class \code{\linkS4class{mvsemodel}}.
#' @param iter A positive integer specifying the number of MCMC steps to run the routine for.
#' @param warmup A real value between 0 and 1 specifying the proportion of MCMC steps to be considered warmup.
#' The warmup samples are discarded and are not used for inference.
#' @param seed The seed for random number generation. The default is generated from 1 to the maximum integer supported by \code{R} on the machine.
#' @param init A named list specifying the initial guesses for the factors in the MCMC routine. Elements include
#' \code{rho}, \code{eta} and \code{alpha} representing the initial guesses of the factors rho, eta and alpha respectively.
#' \code{alpha} does not need to be specified for MVSE models of category \code{"denv_aegypti"}.
#' @param gauJump A named list specifying the standard deviations of the Gaussian proposal distributions for the factors in
#' the MCMC routine. Elements include \code{rho}, \code{eta} and \code{alpha} representing the standard deviations
#' of the Gaussian proposal distributions of rho, eta and alpha respectively. \code{alpha} does not need to be specified
#' for MVSE models of category \code{"denv_aegypti"}.
#' @param smoothing An optional positive integer used for smoothing the climate data. For example,
#' 7 will smooth the time series using +- 7 time points per existing point (up and down that point).
#' @param samples A positive integer specifying the number of samples to draw from the estimated posterior distributions of rho, eta and alpha
#' to estimate the time series of Index P.
#' @param verbose \code{TRUE} or \code{FALSE}: flag indicating whether to print intermediate output from the MCMC routine
#' to the console, which may be helpful for determining progress of algorithm. The default is \code{FALSE}.
#'
#' @return An object of S4 class \code{mvsefit} representing the sampled results for Index P along with additional relevant data.
#'
#' @examples
#' # obtain climate data
#' data(climateFSA)
#'
#' # define a mvse model
#' priors <- list()
#' ### mosquito life expectancy (days)
#' priors$mosq_life_exp <- list(dist="normal", pars=c(mean=12, sd=2))
#' ### mosquito-virus incubation period (days)
#' priors$mosq_inc_per <- list(dist="normal", pars=c(mean=7, sd=2))
#' ### mosquito biting frequency (bites/female/day)
#' priors$mosq_biting_freq <- list(dist="normal", pars=c(mean=0.25, sd=0.01))
#' ### human life expectancy (years)
#' priors$human_life_exp <- list(dist="normal", pars=c(mean=71.1, sd=2))
#' ### human-virus incubation period (days)
#' priors$human_inc_per <- list(dist="normal", pars=c(mean=5.8, sd=1))
#' ### human-virus infectious period (days)
#' priors$human_inf_per <- list(dist="normal", pars=c(mean=5.9, sd=1))
#' user_model <- mvse_model(model_name="my_mvse_model", climate_data=climateFSA, priors=priors)
#'
#' # run the MCMC sampling procedure to estimate the epi-entomological parameters as well Index P
#' user_fit <- sampling(object=user_model, iter=10^5, warmup=0.2, seed=123, 
#' init=c(rho=1, eta=10, alpha=3), samples=1000)
#' user_fit
#'
#' @usage NULL
#'
#' @seealso \code{\linkS4class{mvsemodel}} and \code{\linkS4class{mvsefit}}
#' @name sampling
#' @export
setGeneric(
  name="sampling",
  def=function(object, iter=10^5, warmup=0.5, seed=sample.int(.Machine$integer.max, 1),
               init=c(rho=1, eta=10, alpha=3), gauJump=c(rho=0.5, eta=2, alpha=2),
               smoothing=NULL, samples=1000, verbose=FALSE)
  {standardGeneric("sampling")}
)

#' @rdname sampling
setMethod(f="sampling",
          signature="mvsemodel",
          definition=function(object, iter=10^5, warmup=0.5, seed=sample.int(.Machine$integer.max, 1),
                              init=c(rho=1, eta=10, alpha=3), gauJump=c(rho=0.5, eta=2, alpha=2),
                              smoothing=NULL, samples=1000, verbose=FALSE) {
            set.seed(seed)
            mvse_data <- object@climate_data

            # prepare climate time series
            if (!is.null(smoothing)) {
              mvse_data$H <- (.smoothUDSeries(mvse_data$H, smoothing))
              mvse_data$T <- (.smoothUDSeries(mvse_data$T, smoothing))
            }
            mvse_data$oH<- mvse_data$H
            mvse_data$H<-  (mvse_data$H-min(mvse_data$H, na.rm=TRUE))/(max(mvse_data$H, na.rm=TRUE)-min(mvse_data$H, na.rm=TRUE)) #normalize

            # essential for P
            mvse_data$hum_aV <- .mvse_hum_effect_aV(mvse_data$H, mean(mvse_data$H))
            mvse_data$temp_muV <- .mvse_temp_effect_muV(mvse_data$T)
            mvse_data$hum_muV <- .mvse_hum_effect_muV(mvse_data$H, mean(mvse_data$H))
            mvse_data$temp_epsVH <- .mvse_temp_epsVH(mvse_data$T)
            mvse_data$temp_epsHV <- .mvse_temp_epsHV(mvse_data$T)
            mvse_data$temp_gammaV <- .mvse_temp_effect_gammaV(mvse_data$T)

            # estimate rho, eta and alpha
            priors <- object@priors
            prior_mosq_biting_freq_mean=.getPriorMosqBitingFreqMean(pars=priors$mosq_biting_freq$pars,
                                                                    dist=priors$mosq_biting_freq$dist)
            mcmc_output <- .estimateFactorsRhoEtaAlpha(model_category=object@model_category,
                                                       iter=iter, warmup=warmup, init=init, gauJump=gauJump,
                                                       mvse_data=mvse_data, priors=priors,
                                                       prior_mosq_biting_freq_mean=prior_mosq_biting_freq_mean)
            if (object@model_category=="denv_aegypti") posterior_list <- mcmc_output[c("rho", "eta")]
            else posterior_list <- mcmc_output[c("rho", "eta", "alpha")]

            # estimate time series for index P
            gen_quantities <- .sample_generated_quantities(posterior_list=posterior_list, mvse_data=mvse_data,
                                                           mvsemodel=object, n=samples, verbose=verbose)
            
            # return mvsefit object
            this_mvsefit <- .mvse_fit(model_name=object@model_name, model_category=object@model_category,
                                      mvsemodel=object, mvse_args = list(iter=iter, warmup=warmup, init=init, gauJump=gauJump, seed=seed),
                                      sim=c(list(accepted=mcmc_output$count_accepted, mvse_data=mvse_data), gen_quantities, posterior_list))
            return(this_mvsefit)
          }
)

##########################################################################################################
#' Plot climate data for \code{mvsemodel} objects
#'
#' The default plot shows the time series for temperature and humidity.
#'
#' @param object An object of class \code{\linkS4class{mvsemodel}}.
#' @param vars A vector of strings specifying the climatic variables to plot. Possible arguments include
#' \code{"temperature"}, \code{"humidity"}, or the pair of the two strings.
#' @param smoothing An optional positive integer used for smoothing the climate data. For example,
#' 7 will smooth the time series using +- 7 time points per existing point (up and down that point).
#' @param filename The name of the file where the plot in PNG format is saved (optional).
#'
#' @return A ggplot object that can be further customized using the \bold{ggplot2} package.
#'
#' @examples
#' #define a mvse model
#' data(climateFSA)
#' user_model <- mvse_model(model_name="my_mvse_model", climate_data=climateFSA, priors=list())
#'
#' # plot climate data
#' plot_climate(object=user_model) # temperature and humidity
#' plot_climate(object=user_model, vars=c("temperature"))
#'
#' @usage NULL
#'
#' @seealso \code{\linkS4class{mvsemodel}}
#' @name plot_climate
#' @export
setGeneric(
  name="plot_climate",
  def=function(object, vars=c("temperature", "humidity"), smoothing=NULL, filename=NULL)
  {standardGeneric("plot_climate")}
)

#' @rdname plot_climate
setMethod(f="plot_climate",
          signature="mvsemodel",
          definition=function(object, vars=c("temperature", "humidity"), smoothing=NULL, filename=NULL) {
            data <- object@climate_data
            if (!is.null(smoothing)) {
              data$H <- (.smoothUDSeries(data$H, smoothing))
              data$T <- (.smoothUDSeries(data$T, smoothing))
              data$R <- (.smoothUDSeries(data$R, smoothing))
            }
            data$date <-as.Date(data$date, format='%Y-%m-%d')
            lab_dates <- seq(data$date[1], tail(data$date, 1),
                             by=paste0(round(tail(data$date, 1)-data$date[1])/5, " days"))

            # plot for a single variable
            if (length(vars)==1) {
              if (vars=="temperature") {
                p <- ggplot(data) +
                  theme_bw() +
                  geom_line(aes_string(x="date", y="T"), color='cadetblue4') +
                  scale_y_continuous(limits=c(0, 50), breaks=seq(0, 50, 10)) +
                  labs(y="Temperature")
              }
              else if (vars=="humidity") {
                p <- ggplot(data) +
                  theme_bw() +
                  geom_line(aes_string(x="date", y="H"), color='magenta') +
                  scale_y_continuous(limits=c(0, 100), breaks=seq(0, 100, 20)) +
                  labs(y="Relative humidity")
              }
            }

            # plot for pairs of parameters
            if (length(vars)==2) {
              if (all(sort(vars)==c("humidity", "temperature"))) {
                data <- data %>%
                  mutate(`T` = 2*`T`) %>%
                  gather(key="variable", value="value", c("H", "T")) %>%
                  filter(.data$variable %in% c("H", "T"))
                p <- ggplot(data) +
                  theme_bw() +
                  geom_line(data=data, aes_string(x="date", y="value", color="variable")) +
                  scale_y_continuous(limits=c(0, 100), breaks=seq(0, 100, 20), name="Humidity",
                                     sec.axis = sec_axis(trans=~./2, name="Temperature")) +
                  scale_color_manual(breaks=c("H", "T"), values=c('magenta', 'cadetblue4'),
                                     labels=c("Relative humidity", "Temperature"), guide="none") +
                  theme(legend.position="bottom", legend.title=element_blank(),
                        axis.title.y.left=element_text(color='magenta'),
                        axis.title.y.right=element_text(color='cadetblue4'))
              }
            }
            p <- p +
              scale_x_date(breaks=lab_dates, date_labels="%b %y")
            if (!is.null(filename)) {
              png(filename, width = 1024, height = 768, res=150)
              print(p)
              dev.off()
            }
            return(p)
          }
)

##########################################################################################################
#' Plot the informed probability distributions of the host/vector/virus biological parameters.
#'
#' The default plot shows the probability density function for all parameters.
#'
#' @param object An object of class \code{\linkS4class{mvsemodel}}.
#' @param pars If not NULL, a character vector indicating which parameters to include in the plots. By
#' default, all parameters are included.
#' @param filename The name of the file where the plot in PNG format is saved (optional).
#'
#' @return If only one parameter is specified, a single ggplot object is returned. Otherwise, many ggplot
#' objects organized into a grid via \code{\link{bayesplot_grid}} is returned.
#'
#' @examples
#' #define a mvse model
#' data(climateFSA)
#' priors <- list()
#' ### mosquito life expectancy (days)
#' priors$mosq_life_exp <- list(dist="normal", pars=c(mean=12, sd=2))
#' ### mosquito-virus incubation period (days)
#' priors$mosq_inc_per <- list(dist="normal", pars=c(mean=7, sd=2))
#' ### mosquito biting frequency (bites/female/day)
#' priors$mosq_biting_freq <- list(dist="normal", pars=c(mean=0.25, sd=0.01))
#' ### human life expectancy (years)
#' priors$human_life_exp <- list(dist="normal", pars=c(mean=71.1, sd=2))
#' ### human-virus incubation period (days)
#' priors$human_inc_per <- list(dist="normal", pars=c(mean=5.8, sd=1))
#' ### human-virus infectious period (days)
#' priors$human_inf_per <- list(dist="normal", pars=c(mean=5.9, sd=1))
#' user_model <- mvse_model(model_name="my_mvse_model", climate_data=climateFSA, priors=priors)
#'
#' # plot prior distributions of epi-entomological parameters
#' plot_priors(object=user_model)
#' plot_priors(object=user_model, pars=c("mosq_life_exp", "mosq_inc_per", "human_inf_per"))
#'
#' @usage NULL
#'
#' @seealso \code{\linkS4class{mvsemodel}}
#' @name plot_priors
#' @export
setGeneric(
  name="plot_priors",
  def=function(object, pars=NULL, filename=NULL)
  {standardGeneric("plot_priors")}
)

#' @rdname plot_priors
setMethod(f="plot_priors",
          signature="mvsemodel",
          definition=function(object, pars=NULL, filename=NULL) {
            this_priors <- object@priors
            mosq_life_exp_pdf <- .getDistPDF(pars=this_priors$mosq_life_exp$pars, dist=this_priors$mosq_life_exp$dist)
            if (object@model_category!="denv_aegypti")
              mosq_inc_per_pdf <- .getDistPDF(pars=this_priors$mosq_inc_per$pars, dist=this_priors$mosq_inc_per$dist)
            else
              mosq_inc_per_pdf <- this_priors$mosq_inc_per
            mosq_biting_freq_pdf <- .getDistPDF(pars=this_priors$mosq_biting_freq$pars, dist=this_priors$mosq_biting_freq$dist)
            human_life_exp_pdf <- .getDistPDF(pars=this_priors$human_life_exp$pars, dist=this_priors$human_life_exp$dist)
            human_inc_per_pdf <- .getDistPDF(pars=this_priors$human_inc_per$pars, dist=this_priors$human_inc_per$dist)
            human_inf_per_pdf <- .getDistPDF(pars=this_priors$human_inf_per$pars, dist=this_priors$human_inf_per$dist)
            prior_pdf_list <- list(mosq_life_exp=mosq_life_exp_pdf, mosq_inc_per=mosq_inc_per_pdf,
                                   mosq_biting_freq=mosq_biting_freq_pdf, human_life_exp=human_life_exp_pdf,
                                   human_inc_per= human_inc_per_pdf, human_inf_per=human_inf_per_pdf)
            all_pars <- c("mosq_life_exp", "mosq_inc_per", "mosq_biting_freq", "human_life_exp", "human_inc_per", "human_inf_per")
            labels <- c("Mosquito life expectancy \n (days)", "Mosquito incubation period \n (days)",
                        "Mosquito biting frequency \n (bites/female/day)",
                        "Human life expectancy \n (years)", "Human incubation period \n (days)", "Human infectious period \n (days)")
            names(labels) <- all_pars
            pal <- brewer.pal(length(all_pars), "Set3")
            tmp <- pal[5]; pal[5] <- pal[2]; pal[2] <- tmp
            names(pal) <- all_pars

            get_values <- function(dist, pars) {
              if (dist=="normal") {
                mean <- pars["mean"]
                sd <- pars["sd"]
                values <- seq(mean-3*sd, mean+3*sd, mean/10^3)
                return(values)
              } else if (dist=="lognormal") {
                meanlog <- pars["meanlog"]
                sdlog <- pars["sdlog"]
                mean <- exp(meanlog+sdlog^2/2)
                sd <- sqrt((exp(sdlog^2)-1)*exp(2*meanlog+sdlog^2))
                start <- max(0, mean-3*sd)
                values <- seq(start, mean+3*sd, mean/10^3)
                return(values)
              } else {
                return(0)
              }
            }

            i = 1
            if (is.null(pars)) pars <- all_pars
            plot_list <- list()
            for (par in pars) {
              if (!(par %in% pars)) next
              if (par=="mosq_inc_per" && object@model_category=="denv_aegypti") {
                values_list <- list()
                densities_list <- list()
                temp_list <- list()
                for (ii in seq_along(prior_pdf_list[[par]]$temps)) {

                  values_list[[ii]] <- get_values(dist=prior_pdf_list[[par]]$dists[[ii]]$dist,
                                                  pars=prior_pdf_list[[par]]$dists[[ii]]$pars)
                  this_pdf <- .getDistPDF(dist=prior_pdf_list[[par]]$dists[[ii]]$dist,
                                          pars=prior_pdf_list[[par]]$dists[[ii]]$pars)
                  densities_list[[ii]] <- this_pdf(values_list[[ii]])
                  temp_list[[ii]] <- rep(prior_pdf_list[[par]]$temps[ii], length(values_list[[ii]]))
                }
                values <- do.call(`c`, values_list)
                densities <- do.call(`c`, densities_list)
                temps <- do.call(`c`, temp_list)
                density_df <- do.call(cbind, list(values, densities, temps)) %>%
                  as.data.frame() %>%
                  setNames(c("value", "density", "temp")) %>%
                  mutate(temp=as.factor(.data$temp)) %>%
                  filter(.data$value < 100)
                blue_colors <- colorRampPalette(brewer.pal(8, "Blues"))(length(values_list))
                p <- ggplot(data=density_df) +
                  theme_bw() +
                  geom_ribbon(aes_string(x="value", ymax="density", ymin=0, fill="temp"), alpha=0.8) +
                  geom_line(aes_string(x="value", y="density", color="temp")) +
                  labs(title=labels[par], fill="Temperature", color="Temperature") +
                  scale_fill_manual(breaks=prior_pdf_list[[par]]$temps,
                                    values=blue_colors,
                                    labels=c("    <15\u00B0C",
                                             "15-17.5\u00B0C",
                                             "17.5-20\u00B0C",
                                             "20-22.5\u00B0C",
                                             "22.5-25\u00B0C",
                                             "25-27.5\u00B0C",
                                             "27.5-30\u00B0C",
                                             "30-32.5\u00B0C",
                                             "  >32.5\u00B0C")) +
                  scale_color_manual(breaks=prior_pdf_list[[par]]$temps,
                                     values=blue_colors,
                                     labels=c("    <15\u00B0C",
                                              "15-17.5\u00B0C",
                                              "17.5-20\u00B0C",
                                              "20-22.5\u00B0C",
                                              "22.5-25\u00B0C",
                                              "25-27.5\u00B0C",
                                              "27.5-30\u00B0C",
                                              "30-32.5\u00B0C",
                                              "  >32.5\u00B0C")) +
                  theme(plot.title=element_text(size=8), axis.title=element_blank(),
                        legend.position=c(0.6, 0.5), legend.key.size = unit(0.3, "cm"),
                        legend.title=element_blank())
              }
              else {
                values <- get_values(dist=this_priors[[par]]$dist, pars=this_priors[[par]]$pars)
                density <- prior_pdf_list[[par]](values)
                density_df <- data.frame(x=values, y=density)
                p <- ggplot(data=density_df) +
                  theme_bw() +
                  geom_ribbon(aes_string(x="x", ymax="y", ymin=0), fill=pal[par], alpha=0.8) +
                  geom_line(aes_string(x="x", y="y"), color="black") +
                  labs(title=labels[par]) +
                  theme(plot.title=element_text(size=8), axis.title=element_blank())
              }
              plot_list[[i]] <- p
              i=i+1
            }
            if (!is.null(pars) & length(pars)==1) p <- plot_list[[1]]
            else p <- bayesplot_grid(plots=plot_list, grid_args=list(nrow=ceiling(length(pars)/3)))

            if (!is.null(filename)) {
              png(filename, width = 1024, height = 768, res=150)
              print(p)
              dev.off()
            }
            return(p)
          }
)


#####################################################################
#### INTERNAL FUNCTIONS
#####################################################################

## WRAPPER TO CALL THE MCMC ON RHO, ETA and ALPHA ESTIMATION, ALSO CLIPS OUTPUT ## ACCORDING TO BURNIN PERIOD
.estimateFactorsRhoEtaAlpha <- function(model_category, iter, warmup, init, gauJump, mvse_data, priors,
                                        prior_mosq_biting_freq_mean){
  # set up prior distribution information
  if (model_category=="denv_aegypti") {
    prior_dists <- c(priors$mosq_life_exp$dist, priors$mosq_biting_freq$dist)
    prior_pars <- list(priors$mosq_life_exp$pars, priors$mosq_biting_freq$pars)
  } else {
    prior_dists <- c(priors$mosq_life_exp$dist, priors$mosq_biting_freq$dist, priors$mosq_inc_per$dist)
    prior_pars <- list(priors$mosq_life_exp$pars, priors$mosq_biting_freq$pars, priors$mosq_inc_per$pars)
  }

  # order the factor variables
  factor_order <- c("rho", "eta", "alpha")
  init <- init[order(match(names(init), factor_order))]
  gauJump <- gauJump[order(match(names(gauJump), factor_order))]

  # run sampling procedure
  results <- .MVSErunMCMC_RhoEtaAlpha(model_category=model_category, iter=iter, init=init, gauJump=gauJump,
                                     mvse_data=mvse_data, prior_dists=prior_dists, prior_pars=prior_pars,
                                     prior_mosq_biting_freq_mean=prior_mosq_biting_freq_mean)

  # collate results
  RHOSacc <- results[1,]; RHOSacc <- RHOSacc[!is.na(RHOSacc)]
  ETASacc <- results[2,]; ETASacc <- ETASacc[!is.na(ETASacc)]
  RHOS <- RHOSacc[(length(RHOSacc)*warmup):length(RHOSacc)]
  ETAS <- ETASacc[(length(ETASacc)*warmup):length(ETASacc)]
  if (model_category!="denv_aegypti") {
    ALPHASacc <- results[3,]; ALPHASacc <- ALPHASacc[!is.na(ALPHASacc)]
    ALPHAS <- ALPHASacc[(length(ALPHASacc)*warmup):length(ALPHASacc)]
    output <- list(rho=RHOS, eta=ETAS, alpha=ALPHAS, count_accepted=length(ALPHASacc))
  }
  else
    output <- list(rho=RHOS, eta=ETAS, count_accepted=length(RHOSacc))
  return(output)
}

## USES THE POSTERIOR DISTRIBUTIONS OF RHO, ETA, ALPHA (OPTIONAL) TO GENERATE TIME SERIES
.sample_generated_quantities <- function(posterior_list, mvse_data, mvsemodel, n, verbose=FALSE){
  ss <- sample(1:length(posterior_list$eta), size=n, replace=TRUE)
  if (mvsemodel@model_category!="denv_aegypti") sMCMC_alpha <- posterior_list$alpha[ss]
  sMCMC_rhos <- posterior_list$rho[ss]
  sMCMC_etas <- posterior_list$eta[ss]

  priors <- mvsemodel@priors
  human_life_exp_rng <- .getDistRNG(pars=priors$human_life_exp$pars, dist=priors$human_life_exp$dist)
  human_inc_per_rng <- .getDistRNG(pars=priors$human_inc_per$pars, dist=priors$human_inc_per$dist)
  human_inf_per_rng <- .getDistRNG(pars=priors$human_inf_per$pars, dist=priors$human_inf_per$dist)
  muHs <- 1/(human_life_exp_rng(n)*365)
  gammaHs <- 1/human_inc_per_rng(n)
  deltaHs <- 1/human_inf_per_rng(n)
  prior_mosq_biting_freq_mean=.getPriorMosqBitingFreqMean(pars=priors$mosq_biting_freq$pars, dist=priors$mosq_biting_freq$dist)
  if (mvsemodel@model_category=="denv_aegypti") {
    mosq_inc_per_rng <- lapply(seq_along(priors$mosq_inc_per$temps),
                               function(x) .getDistRNG(pars=priors$mosq_inc_per$dists[[x]]$pars,
                                                       dist=priors$mosq_inc_per$dists[[x]]$dist))
  }

  if (verbose) print("Simulating empirical index P given distributions of rho, eta and alpha...")

  num_steps <- nrow(mvse_data)
  rep_vals <- rep(NA, n*(num_steps))
  indexP <- matrix(rep_vals, ncol=(num_steps))
  muV_samples <- list()
  aV_samples <- list()
  phiVH_samples <- list()
  phiHV_samples <- list()
  gammaV_samples <- list()
  if (verbose) pb <- txtProgressBar(min=0, max = n-1, style = 3)
  for (ii in 1:n){
    # sample from each of factors' distributions
    eta <- sMCMC_etas[ii]
    rho <- sMCMC_rhos[ii]

    # sample other priors
    muH <- muHs[ii]
    gammaH <- gammaHs[ii]
    deltaH <- deltaHs[ii]
    muV_t <- eta*mvse_data$temp_muV*(1+mvse_data$hum_muV)^rho

    # calculate ento-epi parameters
    muV_t <- eta*mvse_data$temp_muV*(1+mvse_data$hum_muV)^rho
    if (mvsemodel@model_category!="denv_aegypti") {
      alpha <- sMCMC_alpha[ii]
      gammaV_t <- alpha*mvse_data$temp_gammaV
    }
    else gammaV_t <- .sample_gammaV(temps=mvse_data$T, groups=priors$mosq_inc_per$temps,
                                    rngs=mosq_inc_per_rng)
    a_t <- prior_mosq_biting_freq_mean*(1+mvse_data$hum_aV)^rho
    epsilonVH_t<- mvse_data$temp_epsVH
    betaVH <- a_t*epsilonVH_t
    epsilonHV_t <- mvse_data$temp_epsHV;
    betaHV<- a_t*epsilonHV_t

    # correct negative to zero for bio meaning
    muV_t[which(muV_t<0)]<- 0
    gammaV_t[which(gammaV_t<0)]<- 0
    betaHV[which(betaHV<0)]<- 0
    betaVH[which(betaVH<0)]<- 0

    # save sampled biological parameters
    muV_samples[[ii]] <- muV_t
    aV_samples[[ii]] <- a_t
    phiVH_samples[[ii]] <- epsilonVH_t
    phiHV_samples[[ii]] <- epsilonHV_t
    gammaV_samples[[ii]] <- gammaV_t

    # calculate index P
    indexP[ii,]<- (betaVH*betaHV*gammaV_t*gammaH)/(muV_t*(deltaH+muH)*(gammaH+muH)*(gammaV_t+muV_t)) #index P
    indexP[ii, which(indexP[ii,]<0)]<- 0

    if (verbose) setTxtProgressBar(pb, ii)
  }
  if (verbose) close(pb)

  muV_samples <- do.call(cbind, muV_samples)
  aV_samples <- do.call(cbind, aV_samples)
  phiVH_samples <- do.call(cbind, phiVH_samples)
  phiHV_samples <- do.call(cbind, phiHV_samples)
  gammaV_samples <- do.call(cbind, gammaV_samples)
  bio_samples <- list(muV=muV_samples, aV=aV_samples, phiVH=phiVH_samples, phiHV=phiHV_samples, gammaV=gammaV_samples)

  reformat <- function(x) {
    x <- t(x)
    colnames(x) <- 1:ncol(x)
    rownames(x) <- NULL
    x <- as.data.frame(x)
    x <- x %>%
      mutate(date=as.Date(mvse_data$date, format='%Y-%m-%d')) %>%
      select(date, dplyr::everything())
    return(x)
  }
  gen_quantities <- lapply(list(indexP), reformat)
  names(gen_quantities) <- c("indexP")
  gen_quantities <- c(gen_quantities, bio_samples)
  return(gen_quantities)
}

## set the prior distributions for the human/vector ecological/epidemiological parameters
.get_aegypti_priors <- function() {
  # intrinsic incubation period
  human_inc_per <- list(dist="lognormal", pars=c("meanlog"=1.738, "sdlog"=1/sqrt(11.386)))
  # mosquito life expectancy
  mosq_life_exp <- list(dist="normal", pars=c("mean"=10, "sd"=2.55))
  # mosquito biting frequency
  mosq_biting_freq <- list(dist="normal", pars=c("mean"=0.25, "sd"=0.01))
  # human life expectancy
  human_life_exp <- list(dist="normal", pars=c("mean"=70, "sd"=3))
  # human infectious period
  human_inf_per <- list(dist="normal", pars=c("mean"=4, "sd"=0.51))
  # extrinsic incubation period
  mosq_inc_per <- list(temps=c(15, 17.5, 20, 22.5, 25, 27.5, 30, 32.5, 50),
                       dists=list(
                         `15`=list(dist="lognormal", pars=c("meanlog"=6.53, "sdlog"=1.08)),
                         `17.5`=list(dist="lognormal", pars=c("meanlog"=5.21, "sdlog"=0.69)),
                         `20`=list(dist="lognormal", pars=c("meanlog"=4.19, "sdlog"=0.43)),
                         `22.5`=list(dist="lognormal", pars=c("meanlog"=3.33, "sdlog"=0.26)),
                         `25`=list(dist="lognormal", pars=c("meanlog"=2.67, "sdlog"=0.18)),
                         `27.5`=list(dist="lognormal", pars=c("meanlog"=2.15, "sdlog"=0.15)),
                         `30`=list(dist="lognormal", pars=c("meanlog"=1.74, "sdlog"=0.15)),
                         `32.5`=list(dist="lognormal", pars=c("meanlog"=1.41, "sdlog"=0.15)),
                         `50`=list(dist="lognormal", pars=c("meanlog"=1.15, "sdlog"=0.15))
                       ))
  return(list(mosq_life_exp=mosq_life_exp, mosq_inc_per=mosq_inc_per, mosq_biting_freq=mosq_biting_freq,
              human_life_exp=human_life_exp, human_inc_per=human_inc_per, human_inf_per=human_inf_per))
}

## samples the gammaV (rate from viremic blood meal to infectiousness)
.sample_gammaV <- function(temps, groups, rngs) {
  samples <- rep(NA, length(temps))
  ii <- 1
  for (temp in temps) {
    temp_grp_id <- min(which(temp < groups))
    if (length(temp_grp_id)==0) temp_grp_id <- length(groups)
    samples[ii] <- 1/(rngs[[temp_grp_id]](1))
    ii <- ii + 1
  }
  return(samples)
}
