#ifndef MOSQUITO_HPP
#define MOSQUITO_HPP


#include <vector>
#include <cstdlib>


#include "hostpopulation.hpp"

class Mosquito{
private:

  double age;
public:
  State state;
  HostPopulation* popn;
  
  void upAge(double step);
  void infect(double cur_t);
  void develop();
  void die();
   
  Mosquito(State _state, HostPopulation* _popn);
  Mosquito(State _state, HostPopulation* _popn, double _age);
  ~Mosquito(){}
  bool isInfected();
  bool isExposed();
  bool isSusceptible();
  bool isRecovered();
  bool isDead();

  State getState();
};

#endif
