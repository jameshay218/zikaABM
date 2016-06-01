#ifndef HUMAN_HPP
#define HUMAN_HPP

#include <vector>
#include <cstdlib>

#include "hostpopulation.hpp"

using namespace std;

class Human{
private:

  double age;
  
  int gestWeek;
  double gestDay;
  int infWeek;
  AgeClass ageState;
  double microCephRisk;

public:
  State state;
  HostPopulation* popn;
  
  void getPregnant();
  void getPregnant(int gestWeek);
  bool isPregnant();
  void upAge(double step, double probPregnant);
  void infect(double cur_t);
  void develop();
  void die();
 
  void recover();
  AgeClass getAge();
  
  Human(State _state, HostPopulation* _popn);
  Human(State _state, HostPopulation* _popn, double _age);
  ~Human(){}
  bool isInfected();
  bool isExposed();
  bool isSusceptible();
  bool isRecovered();
  bool isDead();

  State getState();


};

#endif
