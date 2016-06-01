#ifndef HOSTPOPULATION_HPP
#define HOSTPOPULATION_HPP

#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <math.h>
#include <Rcpp.h>

enum State {Susceptible, Exposed, Infected, Recovered, Dead};
enum AgeClass {Child, Adult, Pregnant};

class Human;
class Mosquito;

class HostPopulation{
private:
  // Vector of pointers to all humans
  std::vector<Human*> susceptible_H;
  std::vector<Human*> infected_H;
  std::vector<Human*> exposed_H;
  std::vector<Human*> recovered_H;

  // Pointers to all mosquitoes
  std::vector<Mosquito*> susceptible_M;
  std::vector<Mosquito*> exposed_M;
  std::vector<Mosquito*> infected_M;

  // Pointers to all new individuals changing classes
  std::vector<Human*> new_infected_H;
  std::vector<Human*> new_exposed_H;
  std::vector<Human*> new_recovered_H;
  std::vector<Human*> new_susceptible_H;

  std::vector<Mosquito*> new_infected_M;
  std::vector<Mosquito*> new_exposed_M;
  std::vector<Mosquito*> new_susceptible_M;

  std::vector<Human*> all_dead_H;
  std::vector<Mosquito*> all_dead_M;
  std::vector<Mosquito*> all_births_M;

  Rcpp::NumericVector microProbs;

  int todayBirths;
  int todayMicro;

  double day;
  double risk;
  double LH, LM, Dc, DEM, DEH, b, pHM, pMH, DIH;
  double birthFreq;


public:
  // Constructors
  HostPopulation();
  HostPopulation(int NH, int NM, Rcpp::NumericVector _m);
  HostPopulation(int NH, int NM, Rcpp::NumericVector _m, double _LH, double _LM, double _Dc, double _DEH, double _DEM, double _b, double _pHM, double _DIH, double _pMH);
  ~HostPopulation();

  // Manage population temporal dynamics
  void stepForward(double step);
  void grow(double step);
  void growUp(double step);
  void contact(double step);
  void progress(double step);
  void pregnancies(double step);
  void recoveries(double step);
  void decline(double step);
  void updateCompartments();
  void burnin(double length);
  void incBirths(bool microCeph);
  void seed(int noSeeds);
  // Get properties of HostPopulation
  int getDay();
  int countSusceptibleH();
  int countExposedH();
  int countInfectedH();
  int countRecoveredH();
  int countH();

  int countSusceptibleM();
  int countExposedM();
  int countInfectedM();
  int countM();

  int countDeadM();
  int countBirthsM();
  int count(AgeClass age, State infState);

  int getTodayBirths();
  int getTodayMicroceph();

  double getMicroprob(int week);
  double getRisk();
  // Print out current population status
  void printStatus();
};

#endif
