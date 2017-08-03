/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
 * 
 * This file contains routines to serially compute the call and 
 * put price of an Asian option.
 * 
 * Simon Scheidegger -- 06/17.
 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/ 

#include <algorithm>    // Needed for the "max" function
#include <cmath>
#include <iostream>
#include <fstream>
#include <omp.h>

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 A simple implementation of the Box-Muller algorithm, used to 
generate gaussian random numbers; necessary for the Monte Carlo 
method below. */

double gaussian_box_muller(unsigned int s1, unsigned int s2) {
  double x = 0.0;
  double y = 0.0;
  double euclid_sq = 0.0;
//  int tnum = omp_get_thread_num();

  // Continue generating two uniform random variables
  // until the square of their "euclidean distance" 
  // is less than unity
  do {
    //unsigned int s1 = tnum*m*omp_get_wtime()*1000;
    //unsigned int s2 = tnum*tnum*m*omp_get_wtime()*1000 + 123954 ;
    x = 2.0 * rand_r(&s1) / static_cast<double>(RAND_MAX)-1;
    y = 2.0 * rand_r(&s2) / static_cast<double>(RAND_MAX)-1;
    euclid_sq = x*x + y*y;
  } while (euclid_sq >= 1.0);

  return x*sqrt(-2*log(euclid_sq)/euclid_sq);
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Pricing a Asian vanilla call option with a Monte Carlo method

double monte_carlo_call_asia_price(const int& num_sims, const double& S, const double& K, const double& r, const double& v, const double& T, const int& num_m) {
  
  double payoff_sum = 0.0;
  double t = T / (num_m*1.0); 
 
  double S_cur;
  double S_prev;
  double S_path;
  double gauss_bm;
  unsigned int s1;
  unsigned int s2;
  int tnum;

  S_cur = 0.0;
  S_prev = S;
  S_path = 0.0;
 
  // loop over the number of  iterations
  #pragma omp parallel
	{ 
  #pragma omp for private(gauss_bm, tnum, s1, s2)\
	firstprivate(S_cur, S_path, S_prev)\
	reduction(+:payoff_sum) 
  for (int i=0; i<num_sims; i++) {
   
   S_cur = 0.0;
   S_prev = S;
   S_path = 0.0;
   tnum = omp_get_thread_num();
   
    // loop over the number of evaluations/periods  per options
   
    for (int m=0; m<num_m; ++m) {
	s1 = tnum*m*omp_get_wtime()*1000;
	s2 = (tnum*tnum-tnum)* m * 1000 * omp_get_wtime(); 
        gauss_bm = gaussian_box_muller(s1, s2);
        S_cur = S_prev * exp(t*(r-0.5*v*v) + v*sqrt(t)*gauss_bm);
        S_path += S_cur;  // add current value S_cur to the path 
	
	S_prev = S_cur; // update S_prev   
        }
   
 
    payoff_sum += std::max((S_path / num_m) - K, 0.0); // calculate payoff for one call
	}
    } // close parallel section
 return (payoff_sum / static_cast<double>(num_sims)) * exp(-r*T);
}




// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Pricing a Asian vanilla put option with a Monte Carlo method

double monte_carlo_put_asia_price(const int& num_sims, const double& S, const double& K, const double& r, const double& v, const double& T, const int& num_m) {
  
  double payoff_sum = 0.0;
  double t = T / (num_m*1.0);
  double S_adjust = S * exp(t*(r-0.5*v*v)); 
 
  // loop over the number of  iterations
  for (int i=0; i<num_sims; i++) {
   
    double S_cur = 0.0;
    double S_prev = S;
    double S_path = 0.0;
    
    // loop over the number of evaluations/periods  per options
    for (int m=0; m<num_m; ++m) {
 	unsigned int s1 = 1235*m*omp_get_wtime()*1000;
	unsigned int s2 = (m*m -m) * 1000 * omp_get_wtime(); 
        double gauss_bm = gaussian_box_muller(s1, s2);
        S_cur = S_prev * exp(t*(r-0.5*v*v) + v*sqrt(t)*gauss_bm);
        S_path += S_cur;  // add current value S_cur to the path 
	
	S_prev = S_cur; // update S_prev   
        }
    payoff_sum += std::max(K -  (S_path / num_m), 0.0); // calculate payoff for one call
  }

  return (payoff_sum / static_cast<double>(num_sims)) * exp(-r*T);
}


// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
int main(int argc, char **argv) {

  // Parameters                                                                             
  int num_sims = 10000000;   // Number of simulated asset paths                                                       
  int num_m = 10;     // number of evaluations in one simulated asset path

  double S = 100.0;  // Option price                                                                                  
  double K = 100.0;  // Strike price                                                                                  
  double r = 0.05;   // Risk-free rate (5%)                                                                           
  double v = 0.2;    // Volatility of the underlying (20%)                                                            
  double T = 1.0;    // One year until expiry                                                                         

  // Then we calculate the call/put values via Monte Carlo                                                                          
  double call_time_parallel = -omp_get_wtime();
  double call_asia = monte_carlo_call_asia_price(num_sims, S, K, r, v, T, num_m);
  call_time_parallel += omp_get_wtime();

  double put_time_serial = -omp_get_wtime();
  double put_asia = monte_carlo_put_asia_price(num_sims, S, K, r, v, T, num_m);
  put_time_serial += omp_get_wtime();
  //double put = monte_carlo_put_price(num_sims, S, K, r, v, T);

  // Finally we output the parameters and prices                                                                  
  std::cout << "\nPricing an Asian Put and Call option" << std::endl; 
  std::cout << "Number of Paths:   " << num_sims << std::endl;
  std::cout << "# of evals p Path: " << num_m << std::endl;
  std::cout << "Underlying:        " << S << std::endl;
  std::cout << "Strike:            " << K << std::endl;
  std::cout << "Risk-Free Rate:    " << r << std::endl;
  std::cout << "Volatility:        " << v << std::endl;
  std::cout << "Maturity:          " << T << std::endl;

  std::cout << "Call Price:        " << call_asia << std::endl;
  std::cout << "Call time parallel:" << call_time_parallel << std::endl;
  std::cout << "Put Price:         " << put_asia << std::endl;
  std::cout << "Put time serial:   " << put_time_serial << std::endl;
  
  
  return 0;
}
