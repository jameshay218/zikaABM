#include "human.hpp"
using namespace std;

void Human::getPregnant(){
  // Change to pregnant state and set gestation times to start
  ageState = Pregnant;
  gestWeek = 0;
  gestDay = 0;
}

void Human::getPregnant(int gestWeek){
  ageState = Pregnant;
  gestWeek = gestWeek;
  gestDay = R::unif_rand()*7;
}


bool Human::isPregnant(){
  return(ageState == Pregnant);
}

AgeClass Human::getAge(){
  return(ageState);
}

Human::Human(State _state, HostPopulation* _popn, double _age){
  age = _age;
  state = _state;
  popn = _popn;
  gestWeek = -1;
  gestDay = -1;
  infWeek = -1;
  microCephRisk = 0;
  if(_age >= 18.0) ageState = Adult;
  else ageState = Child;
}



Human::Human(State _state, HostPopulation* _popn){
  age = 0;
  state = _state;
  popn = _popn;
  gestWeek = -1;
  gestDay = -1;
  infWeek = -1;
  microCephRisk = 0;
  ageState = Child;
}

void Human::upAge(double step, double probPregnant){
  // Increase age by a certain amount
  age += step;
  // If age is greater than 18 years, change this Human to be an adult
  if(ageState == Child && age >= 18.0) ageState = Adult;
  if(ageState == Adult && R::unif_rand() <= probPregnant) getPregnant();
  
  // Increase day of gestation by step amount (convert from years to days)
  if(ageState == Pregnant){
    gestDay += step;
    // If at the end of the week, then increase gestation week, and reset gestation day to start of week
    if(gestDay >= 7){
      gestWeek++;
      gestDay = gestDay - 7;
    }
    // If reached the end of the week, then return to being an adult
    // THIS IS WHERE BIRTHS WILL BE EXPLICITLY CALCULATED
    if(gestWeek >= 40){
      popn->incBirths(R::unif_rand() <= microCephRisk);
      ageState = Adult;
      gestWeek = -1;
      gestDay = -1;
    }
  }
}

void Human::infect(double cur_t){
  state = Exposed;
  if(ageState == Pregnant){
    infWeek = gestWeek;
    microCephRisk = popn->getMicroprob(infWeek);
  }
  // Will need to calculate risk of microceph here
  // microCephRisk = ...
}


void Human::recover(){
  state = Recovered;
}


void Human::develop(){
  state = Infected;
}

void Human::die(){
  state = Dead;
}

State Human::getState(){
  return(state);
}

bool Human::isExposed(){
  return(state==Exposed);
}

bool Human::isInfected(){
  return(state == Infected);
}

bool Human::isDead(){
  return(state==Dead);
}

bool Human::isSusceptible(){
  return(state==Susceptible);
}

bool Human::isRecovered(){
  return(state==Recovered);
}



