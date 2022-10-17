
# MVSE

MVSE, Mosquito-borne Viral Suitability Estimator, provides methods for
the estimation of an index of climate-based transmission potential for
mosquito-borne viruses. It can be parameterized for any combination of
host, virus and mosquito-species of interest with local climate data.

## Installation

The MVSE package is currently not hosted on CRAN (Comprehensive R
Archive Network). We thus need to install the package directly from the
GithHub repo where it is hosted.

``` r
install.packages("devtools")
install.packages("tidyverse")
require(tidyverse)
require(devtools)

install_github("TaishiNakase/MVSE")
```

## Example

Let’s go through a short example of how one might use MVSE to estimate
Index P for some mosquito-borne virus.

First, we will install some packages.

``` r
library(MVSE)
```

We first need to create a `mvse_model` object, which requires time
series data for climatic variables (i.e. temperature, humidity and
rainfall) as well as epi-entomological prior distributions for the
virus/vector/host system of interest.

``` r
# climate data
data("climateFSA")

# user-defined model
priors <- list(mosq_life_exp=list(pars=c(mean=12, sd=2), dist="normal"),
               mosq_inc_per=list(pars=c(mean=7, sd=2), dist="normal"),
               mosq_biting_freq=list(pars=c(mean=0.25, sd=0.01), dist="normal"),
               human_life_exp=list(pars=c(mean=71.1, sd=2), dist="normal"),
               human_inc_per=list(pars=c(mean=5.8, sd=1), dist="normal"),
               human_inf_per=list(pars=c(mean=5.9, sd=1), dist="normal"))
example_mvse_model <- mvse_model(model_name="user_test", climate_data=climateFSA, priors=priors)
```

Let’s take a look at the input data to be used in the sampling
procedure.

``` r
example_mvse_model
#> Class: mvsemodel 
#> Model name: user_test 
#> Model category: user-defined 
#> Climate data (limited to first 10 rows): 
#>          date        T        H            R
#> 1  1981-01-02 25.13554 73.84253 0.0022509235
#> 2  1981-02-02 25.01051 73.56940 0.0017562567
#> 3  1981-03-02 25.44149 76.77259 0.0091257023
#> 4  1981-04-02 23.51202 84.59996 0.0029805890
#> 5  1981-05-02 22.42301 85.48878 0.0020460677
#> 6  1981-06-02 21.85287 85.16210 0.0020688535
#> 7  1981-07-02 20.71406 83.57096 0.0017366349
#> 8  1981-08-02 20.99485 79.70238 0.0012274350
#> 9  1981-09-02 22.46040 70.62206 0.0005982321
#> 10 1981-10-02 25.42501 65.91237 0.0006568354
#> Priors: 
#>   Mosquito life expectancy (days)        : normal(mean=12, sd=2) 
#>   Mosquito incubation period (days)      : normal(mean=7, sd=2) 
#>   Mosquito biting frequency (bites/female/day)   : normal(mean=0.25, sd=0.01) 
#>   Human life expectancy (years)          : normal(mean=71.1, sd=2) 
#>   Human incubation period (days)         : normal(mean=5.8, sd=1) 
#>   Human infectious period (days)         : normal(mean=5.9, sd=1)
```

Next, let’s perform the MCMC (Markov chain Monte Carlo) sampling
procedure to estimate Index P.

``` r
example_mvse_fit <- fitting(example_mvse_model, verbose=FALSE)
```

Finally, let’s take a look at the estimated distribution of the time
series of Index P.

``` r
indexP_plot <- mcmc_index_dist(example_mvse_fit, index="indexP")
```
