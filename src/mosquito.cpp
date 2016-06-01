#include "mosquito.hpp"

using namespace std;

Mosquito::Mosquito(State _state, HostPopulation* _popn){
  age = 0;
  state = _state;
  popn = _popn;
}

Mosquito::Mosquito(State _state, HostPopulation* _popn, double _age){
  age = _age;
  state = _state;
  popn = _popn;
}

void Mosquito::infect(double day){
  state = Exposed;
}

void Mosquito::develop(){
  state = Infected;
}

void Mosquito::die(){
  state = Dead;
}

void::Mosquito::upAge(double step){
  age += step;
}

State Mosquito::getState(){
  return(state);
}

bool Mosquito::isExposed(){
  return(state==Exposed);
}

bool Mosquito::isInfected(){
  return(state == Infected);
}

bool Mosquito::isDead(){
  return(state==Dead);
}

bool Mosquito::isSusceptible(){
  return(state==Susceptible);
}

bool Mosquito::isRecovered(){
  return(state==Recovered);
}



