#include <Rcpp.h>
#include <vector>
#include <limits>
using namespace Rcpp;

// function pointer for the pdf functions
typedef double (*fptr)(double x, NumericVector& pars);

// takes a double for the value and NumericVector of parameters
// and returns the log of the normal probability density. 
double normal_log_pdf(double x, NumericVector& pars) {
	return R::dnorm(x, pars[0], pars[1], TRUE);
}

// takes a double for the value and NumericVector of parameters
// and returns the log of the lognormal probability density. 
double lognormal_log_pdf(double x, NumericVector& pars) {
	return R::dlnorm(x, pars[0], pars[1], TRUE);
}

// takes a string for the distribution type and returns 
// its corresponding pdf function. 
fptr log_pdf(String dist) {
	if (dist=="normal") 
		return normal_log_pdf;
	else 
		return lognormal_log_pdf;
}

// [[Rcpp::export]]
NumericMatrix MVSErunMCMC_RhoEtaAlpha(String model_category, unsigned int iter, Rcpp::NumericVector& init, Rcpp::NumericVector& gauJump, 
	                                  Rcpp::DataFrame& mvse_data, Rcpp::StringVector& prior_dists, Rcpp::List& prior_pars, 
	                                  double prior_mosq_biting_freq_mean) {
	// defining enum Factor
	enum Factor {RHO, ETA, ALPHA};
	enum Prior {MLE, MBF, MEIP}; // mosq. life exp., mosq. biting freq. and mosq. EIP

	// proposals
	double pL;     // likelihood of proposal
	double pRho;   // proposed rho factor
	double pEta;   // proposed eta factor
	double pAlpha; // proposed alpha factor

	// current
	double cL = -std::numeric_limits<double>::infinity();
	double cRho = init[RHO];
	double cEta = init[ETA];
	double cAlpha;
	if (model_category!="aegypti") {
		cAlpha = init[ALPHA];
	}

	// standard deviations of kernels
	double rho_sd = gauJump[RHO];
	double eta_sd = gauJump[ETA];
	double alpha_sd;
	if (model_category!="aegypti") {
		alpha_sd = gauJump[ALPHA];
	} 

	// climate effect data
	Rcpp::NumericVector Z = mvse_data["Z"];
	Rcpp::NumericVector Y = mvse_data["Y"];
	Rcpp::NumericVector A = mvse_data["A"];
	Rcpp::NumericVector G = mvse_data["G"];

	// prior distributions
	NumericVector mosq_life_exp_pars = prior_pars[MLE];
	NumericVector mosq_biting_freq_pars = prior_pars[MBF];
	NumericVector mosq_inc_per_pars;
	if (model_category!="aegypti") mosq_inc_per_pars = prior_pars[MEIP];
	fptr mosq_life_exp_pdf = log_pdf(prior_dists[MLE]);
	fptr mosq_biting_freq_pdf = log_pdf(prior_dists[MBF]);
	fptr mosq_inc_per_pdf;
	if (model_category!="aegypti") mosq_inc_per_pdf = log_pdf(prior_dists[MEIP]);

	// intermediate and final outputs
	NumericMatrix accepted(3, iter);
	std::fill(accepted.begin(), accepted.end(), NumericVector::get_na());
	size_t mvse_data_num_rows = Z.length();
	double mu, a, gamma;
	double p1, p2, p3, pAccept;

	// run Metropolis-Hastings algorithm
	for (int ii = 0; ii < iter; ii++) {
		// generate proposals
		while (true) {
			pRho = R::rnorm(cRho, rho_sd);
			if (pRho>0 && pRho<10) break;
		}
		while (true) {
			pEta = R::rnorm(cEta, eta_sd);
			if (pEta>0 && pEta<20) break;
		}

		// compute mean mosq. life exp. and mean mosq. biting rate
		a = 0;
		mu = 0;
		for (size_t jj = 0; jj < mvse_data_num_rows; jj++) {
			mu += 1/(pEta*Z[jj]*pow((1+Y[jj]), pRho));
			a += prior_mosq_biting_freq_mean*pow((1+A[jj]), pRho);
		}
		mu = mu/(1.0*mvse_data_num_rows);
		a = a/(1.0*mvse_data_num_rows);



		if (model_category!="aegypti") {
			// generate proposal
			while (true) {
				pAlpha = R::rnorm(cAlpha, alpha_sd);
				if (pAlpha>0 && pAlpha<10) break;
			}
			// compute mean EIP
			gamma = 0;
			for (size_t jj = 0; jj < mvse_data_num_rows; jj++) {
				gamma += 1/(pAlpha*G[jj]);
			}
			gamma = gamma/(1.0*mvse_data_num_rows);
		}
		
		// compute acceptance probability
		p1 = mosq_life_exp_pdf(mu, mosq_life_exp_pars);
		p2 = mosq_biting_freq_pdf(a, mosq_biting_freq_pars);
		p3 = 0;
		if (model_category!="aegypti") p3 = mosq_inc_per_pdf(gamma, mosq_inc_per_pars);
		pL = p1 + p2 + p3;
		pAccept = exp(pL-cL);

		// update if accepted
		double event = R::runif(0, 1);
		if (event < pAccept) {
			accepted(RHO, ii) = pRho;
			accepted(ETA, ii) = pEta;
			if (model_category!="aegypti") {
				accepted(ALPHA, ii) = pAlpha;
				cAlpha = pAlpha;
			}
			cRho = pRho;
			cEta = pEta;
			cL = pL;
		}
    }

    return accepted;
}