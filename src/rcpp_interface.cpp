// [[Rcpp::plugins(cpp11)]]

#include <Rcpp.h>
#include <vector>

#include "hostpopulation.hpp"
#include "human.hpp"
#include "mosquito.hpp"


using namespace std;
using namespace Rcpp;
//' @export
//' @useDynLib zikaABM
//[[Rcpp::export]]
NumericMatrix run_simulation(double step, double final_t, int NH, int NM, NumericVector microProbs, double burnin, int seed, NumericVector pars){
  int noRows = final_t/step;
  NumericMatrix timeSeries(noRows,10);
  HostPopulation hpop(NH,NM, microProbs, pars[0], pars[1], pars[2],pars[3],pars[4],pars[5],pars[6],pars[7],pars[8]);
  double t = 0;
  int i = 0;
  int j = 0;
  int tmp;
  Rcpp::Rcout << "Burnin..." << endl;
  hpop.burnin(burnin);
  Rcpp::Rcout << "Burnin done" << endl;
  hpop.printStatus();
  Rcpp::Rcout << endl;
  hpop.seed(seed);
 
  while(i < noRows){
    hpop.stepForward(step);
    hpop.printStatus();
    Rcpp::Rcout << endl;
    t += step;
    timeSeries(i,0) = hpop.countSusceptibleH();
    timeSeries(i,1) = hpop.countExposedH();
    timeSeries(i,2) = hpop.countInfectedH();
    timeSeries(i,3) = hpop.countRecoveredH();
    timeSeries(i,4) = hpop.countSusceptibleM();
    timeSeries(i,5) = hpop.countExposedM();
    timeSeries(i,6) = hpop.countInfectedM();

    timeSeries(i,7) = hpop.getTodayBirths();
    timeSeries(i,8) = hpop.getTodayMicroceph();
    timeSeries(i,9) = hpop.getRisk();
    tmp = hpop.count(Pregnant,Susceptible) + hpop.count(Pregnant, Exposed) + hpop.count(Pregnant,Infected) + hpop.count(Pregnant, Recovered);
    Rcpp::Rcout << "Number pregnant adults: " << tmp << endl << endl;;
    i++;
  }
  return(timeSeries);
}
