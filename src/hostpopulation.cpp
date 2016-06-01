#include "hostpopulation.hpp"
#include "mosquito.hpp"
#include "human.hpp"

using namespace std;

HostPopulation::HostPopulation(){
  day = 0;
  LH = 70.0;
  LM = 14.0/365.0;
  Dc = 18.0;
  DEH = 5.0/365.0;
  DEM = 5.0/365.0;
  b = 100.0;
  pHM = 0.5;
  DIH = 5.0/365.0;
  pMH = 0.5;
  todayBirths=0;
  todayMicro=0;
  birthFreq = 7.0/365.0;

  for(int i = 0; i < 1000000; ++i){
    Human* firstH = new Human(Susceptible, this, 18.0);
      susceptible_H.push_back(firstH);
  }
 
  for(int i = 0; i < 3000000; ++i){
    Mosquito* firstM = new Mosquito(Susceptible, this);
      susceptible_M.push_back(firstM);
  }
}


HostPopulation::HostPopulation(int NH, int NM, Rcpp::NumericVector _m){

  microProbs = _m;
  day = 0;
  LH = 70.0;
  LM = 14.0/365.0;
  Dc = 18.0;
  DEH = 5.0/365.0;
  DEM = 5.0/365.0;
  b = 100.0;
  pHM = 0.5;
  DIH = 5.0/365.0;
  pMH = 0.5;

  todayBirths=0;
  todayMicro=0;
  birthFreq = 7.0/365.0;

  for(int i = 0; i < NH; ++i){
    Human* firstH = new Human(Susceptible, this,R::unif_rand()*LH);
    susceptible_H.push_back(firstH);
  }
  
  for(int i = 0; i < NM; ++i){
    Mosquito* firstM = new Mosquito(Susceptible, this, R::unif_rand()*LM);
    susceptible_M.push_back(firstM);
  }
}


HostPopulation::HostPopulation(int NH, int NM, Rcpp::NumericVector _m, double _LH, double _LM, double _Dc, double _DEH, double _DEM, double _b, double _pHM, double _DIH, double _pMH){
  microProbs = _m;
  day = 0;
  LH = _LH;
  LM = _LM;
  Dc = _Dc;
  DEH = _DEH;
  DEM = _DEM;
  b = _b;
  pHM = _pHM;
  DIH = _DIH;
  pMH = _pMH;
  risk = 0;

  todayBirths=0;
  todayMicro=0;

  int noChildren = NH * Dc/(LH+Dc);
  int noAdults = NH * LH/(LH+Dc);
  int noPregnant = noAdults * 1/(LH-Dc);
  double probPreg = noPregnant/noAdults;  
  for(int i = 0; i < noChildren; ++i){
    Human* firstH = new Human(Susceptible, this,R::unif_rand()*Dc);
    susceptible_H.push_back(firstH);
  }
  for(int i = 0; i < noAdults; ++i){
    Human* firstH = new Human(Susceptible, this,((R::unif_rand()*LH)+Dc));
    if(R::unif_rand() < probPreg) firstH->getPregnant(R::unif_rand()*39);
    susceptible_H.push_back(firstH);
  }
  
  for(int i = 0; i < NM; ++i){
    Mosquito* firstM = new Mosquito(Susceptible, this, R::unif_rand()*LM);
    susceptible_M.push_back(firstM);
  }
}

void HostPopulation::burnin(double length){
  double step = 1.0;
  while(day < length){
    if(fmod(day,30) < step){
      Rcpp::Rcout << "Burn in year: " << day << endl;
      Rcpp::Rcout << "Number of pregnant: " << count(Pregnant, Susceptible) << endl;
    }
    day += step;
    grow(step);
    //pregnancies(step);
    growUp(step);
    decline(step);
    updateCompartments();
  }
  
}

void HostPopulation::seed(int noSeeds){
  for(int i = 0; i < noSeeds; ++i){
    Human* firstI = new Human(Infected, this, R::unif_rand()*LH);
    infected_H.push_back(firstI);
  }
}

void HostPopulation::stepForward(double step){
  // Current day is changed
  day += step;
  
  todayBirths = 0;
  todayMicro = 0;
 
  //New births of susceptibles and moving around age classes etc
  contact(step);
  grow(step);
  //pregnancies(step);
  growUp(step);
  // Infected population grows from transmission
 

  progress(step);

  // Infected population declines from recovery
  recoveries(step);

  // Deaths from all compartments - MUST BE LAST EVENT
  decline(step);
  
  //cout << "Update" << endl;
  updateCompartments();

}


void HostPopulation::grow(double step){
  // Generate new Human births
  int newBirths = R::rpois((step/LH)*countH());
  for(int i = 0; i < newBirths; ++i){
    Human* h = new Human(Susceptible, this);
    new_susceptible_H.push_back(h);
  }
  
  // Generate new Mosquito births
  newBirths = R::rpois(step*countM()/LM);
  for(int i = 0; i < newBirths; ++i){
    Mosquito* m = new Mosquito(Susceptible, this);
    new_susceptible_M.push_back(m);
  }
}

void HostPopulation::growUp(double step){
  double probPregnant = step/(LH-Dc);
  vector<vector<Human*> > allHumans = {susceptible_H, exposed_H, infected_H, recovered_H};
  vector<vector<Mosquito*> > allM= {susceptible_M, exposed_M, infected_M};
  int max = allHumans.size();
  for(int i = 0; i < max; ++i){
    for(vector<Human*>::iterator it = allHumans[i].begin(); it != allHumans[i].end(); ++it){
      (*it)->upAge(step, probPregnant);
    }
  }
  max = allM.size();
  for(int i = 0; i < max; ++i){
    for(vector<Mosquito*>::iterator it = allM[i].begin(); it != allM[i].end(); ++it){
      (*it)->upAge(step);
    }
  }

}

void HostPopulation::pregnancies(double step){
  // Probability of becoming pregnant at a given time point
  double probPregnant = step/(LH-Dc);
 
  // Put all humans into same vector
  vector<vector<Human*> > allHosts = {susceptible_H, exposed_H, infected_H, recovered_H};
  // Go through all hosts. Each host has a probability of becoming pregnant, assuming that they are a non-pregnant adult
  int max = allHosts.size();
  for(int i = 0; i < max; ++i){
    for(vector<Human*>::iterator it = allHosts[i].begin(); it != allHosts[i].end(); ++it){
      if((*it)->getAge() == Adult && R::unif_rand() <= probPregnant){
	(*it)->getPregnant();
      }
    }
  }
}

void HostPopulation::contact(double step){
  double lambda_H = b*pMH*countInfectedM()*step/countH();
  double lambda_M = b*pHM*countInfectedH()*step/countH();
  int index = 0;
  int noInfected = R::rpois(lambda_H*countSusceptibleH());
  risk = lambda_H;
  Rcpp::Rcout << "b: " << b << endl;
  Rcpp::Rcout << "Risk of infection: " << lambda_H << endl;
  Rcpp::Rcout<< "Incidence: " << double(countInfectedH())/double(countH()) << endl;
  for(int i = 0; i < noInfected; ++i){
    if(countSusceptibleH() > 0){
      index = floor(R::unif_rand()*countSusceptibleH());
      new_exposed_H.push_back(susceptible_H[index]);
      susceptible_H[index]->infect(day);
      susceptible_H[index] = susceptible_H.back();
      susceptible_H.pop_back();
    }
  }
  susceptible_H.insert(susceptible_H.end(), new_exposed_H.begin(),new_exposed_H.end());

  noInfected = R::rpois(lambda_M*countSusceptibleM());
  for(int i = 0; i < noInfected; ++i){
    if(countSusceptibleM() > 0){
      index = floor(R::unif_rand()*countSusceptibleM());
      new_exposed_M.push_back(susceptible_M[index]);
      susceptible_M[index]->infect(day);
      susceptible_M[index] = susceptible_M.back();
      susceptible_M.pop_back();
    }
  }
  susceptible_M.insert(susceptible_M.end(), new_exposed_M.begin(),new_exposed_M.end());
 
  // All susceptible humans have a chance of becoming infected each day
  /*  int max = susceptible_H.size();
  for(int i = 0; i < max; ++i){
    if(R::unif_rand() <= lambda_H){
      susceptible_H[i]->infect(day);
      new_exposed_H.push_back(susceptible_H[i]);
    }
  }
  
  // All susceptible mosquitoes have a chance of becoming infected
  for(vector<Mosquito*>::iterator it = susceptible_M.begin(); it != susceptible_M.end(); ++it){
    if(R::unif_rand() <= lambda_M){
      (*it)->infect(day);
      new_exposed_M.push_back(*it);
    }
  }
  */
}
  
void HostPopulation::progress(double step){
  int index = 0;
  int noProgress = R::rpois(step*countExposedH()/DEH);
  for(int i = 0; i < noProgress; ++i){
    if(countExposedH() > 0){
      index = floor(R::unif_rand()*countExposedH());
      new_infected_H.push_back(exposed_H[index]);
      exposed_H[index]->develop();
      exposed_H[index] = exposed_H.back();
      exposed_H.pop_back();
    }
  }
  exposed_H.insert(exposed_H.end(), new_infected_H.begin(),new_infected_H.end());
 
  noProgress = R::rpois(countExposedM()*step/DEM);
  for(int i = 0; i < noProgress; ++i){
    if(countExposedM() > 0){
      index = floor(R::unif_rand()*countExposedM());
      new_infected_M.push_back(exposed_M[index]);
      exposed_M[index]->develop();
      exposed_M[index] = exposed_M.back();
      exposed_M.pop_back();
    }
  }
  exposed_M.insert(exposed_M.end(), new_infected_M.begin(),new_infected_M.end());
  
 
  /*
    double probProgress = step/DEH;
    for(vector<Human*>::iterator it = exposed_H.begin(); it != exposed_H.end(); ++it){
    if(R::unif_rand() <= probProgress){
    (*it)->develop();
    new_infected_H.push_back(*it);
    }
    }
    probProgress = R::rpois(countExposedM()*step/DEM);
    for(vector<Mosquito*>::iterator it = exposed_M.begin(); it != exposed_M.end(); ++it){
    if(R::unif_rand() <= probProgress){
    (*it)->develop();
    new_infected_M.push_back(*it);
    }
    }*/
}

void HostPopulation::recoveries(double step){
  double probRecover = step/DIH;
  int noRecover = R::rpois(probRecover*countInfectedH());
  int index = 0;
  for(int i = 0; i < noRecover; ++i){
    if(countInfectedH() > 0){
      index = floor(R::unif_rand()*countInfectedH());
      new_recovered_H.push_back(infected_H[index]);
      infected_H[index]->recover();
      infected_H[index] = infected_H.back();
      infected_H.pop_back();
    }
  }
  infected_H.insert(infected_H.end(), new_recovered_H.begin(),new_recovered_H.end());
  
  
  /*
  for(vector<Human*>::iterator it = infected_H.begin(); it != infected_H.end(); ++it){
    if(R::unif_rand() <= probRecover){
      (*it)->recover();
      new_recovered_H.push_back(*it);
    }
    }*/
}

void HostPopulation::decline(double step){
  // HUMANS
  // Go through all classes and a certain number die
  int index = 0;
  double probDeath = step/LH;
  int noDead = R::rpois(countSusceptibleH()*probDeath);
  for(int i = 0; i < noDead; ++i){
    if(countSusceptibleH() > 0){
      index = floor(R::unif_rand()*countSusceptibleH());
      all_dead_H.push_back(susceptible_H[index]);
      susceptible_H[index]->die();
      susceptible_H[index] = susceptible_H.back();
      susceptible_H.pop_back();
    }
  }
  noDead = R::rpois(countExposedH()*probDeath);
  for(int i = 0; i < noDead; ++i){
    if(countExposedH() > 0){
      index = floor(R::unif_rand()*countExposedH());
      all_dead_H.push_back(exposed_H[index]);
      exposed_H[index]->die();
      exposed_H[index] = exposed_H.back();
      exposed_H.pop_back();
    }
  }
  noDead = R::rpois(countInfectedH()*probDeath);
  for(int i = 0; i < noDead; ++i){
    if(countInfectedH() > 0){
      index = floor(R::unif_rand()*countInfectedH());
      all_dead_H.push_back(infected_H[index]);
      infected_H[index]->die();
      infected_H[index] = infected_H.back();
      infected_H.pop_back();
    }
  }
  noDead = R::rpois(countExposedH()*probDeath);
  for(int i = 0; i < noDead; ++i){
    if(countRecoveredH() > 0){
      index = floor(R::unif_rand()*countRecoveredH());
      all_dead_H.push_back(recovered_H[index]);
      recovered_H[index]->die();
      recovered_H[index] = recovered_H.back();
      recovered_H.pop_back();
    }
  }

  // MOSQUITOES
  // Go through all classes and a certain number die
  probDeath = step/LM;
  noDead = R::rpois(countSusceptibleM()*probDeath);
  for(int i = 0; i < noDead; ++i){
    if(countSusceptibleM() > 0){
      index = floor(R::unif_rand()*countSusceptibleM());
      all_dead_M.push_back(susceptible_M[index]);
      susceptible_M[index]->die();
      susceptible_M[index] = susceptible_M.back();
      susceptible_M.pop_back();
    }
  }
  noDead = R::rpois(countExposedM()*probDeath);
  for(int i = 0; i < noDead; ++i){
    if(countExposedM() > 0){
      index = floor(R::unif_rand()*countExposedM());
      all_dead_M.push_back(exposed_M[index]);
      exposed_M[index]->die();
      exposed_M[index] = exposed_M.back();
      exposed_M.pop_back();
    }
  }
  noDead = R::rpois(countInfectedM()*probDeath);
  for(int i = 0; i < noDead; ++i){
    if(countInfectedM() > 0){
      index = floor(R::unif_rand()*countInfectedM());
      all_dead_M.push_back(infected_M[index]);
      infected_M[index]->die();
      infected_M[index] = infected_M.back();
      infected_M.pop_back();
    }
  }

  /*  
      vector<vector<Human*> > allHumans = {susceptible_H, exposed_H, infected_H, recovered_H};
      int max = allHumans.size();
      for(int i = 0; i < max; ++i){
      for(vector<Human*>::iterator it = allHumans[i].begin(); it != allHumans[i].end(); ++it){
      if(R::unif_rand() <= probDeath){
      (*it)->die();
      all_dead_H.push_back(*it);
      }
      }
      }
      probDeath = step/LM;
      vector<vector<Mosquito*> > allMos = {susceptible_M, exposed_M, infected_M};
      max = allMos.size();
      int greb = 0;
      for(int i = 0; i < max; ++i){
      greb = allMos[i].size();
      //    for(vector<Mosquito*>::iterator it = allMos[i].begin(); it != allMos[i].end(); ++it){
      for(int j = 0; j < greb; ++j){
      if(R::unif_rand() <= probDeath){
      allMos[i][j]->die();
      all_dead_H.push_back(*it);
      }
      }
      }*/
}


void HostPopulation::updateCompartments(){
  // Insert all new compartment members into their correct compartment
  susceptible_H.insert(susceptible_H.end(), new_susceptible_H.begin(), new_susceptible_H.end());
  susceptible_M.insert(susceptible_M.end(), new_susceptible_M.begin(), new_susceptible_M.end());

  exposed_H.insert(exposed_H.end(), new_exposed_H.begin(), new_exposed_H.end());
  exposed_M.insert(exposed_M.end(), new_exposed_M.begin(), new_exposed_M.end());

  infected_H.insert(infected_H.end(), new_infected_H.begin(), new_infected_H.end());
  infected_M.insert(infected_M.end(), new_infected_M.begin(), new_infected_M.end());

  recovered_H.insert(recovered_H.end(), new_recovered_H.begin(), new_recovered_H.end());

  susceptible_H.erase(remove_if(susceptible_H.begin(), susceptible_H.end(), [](Human *h){return(!h->isSusceptible());}), susceptible_H.end());
  susceptible_M.erase(remove_if(susceptible_M.begin(), susceptible_M.end(), [](Mosquito *h){return(!h->isSusceptible());}), susceptible_M.end());
 
  exposed_H.erase(remove_if(exposed_H.begin(), exposed_H.end(), [](Human *h){return(!h->isExposed());}), exposed_H.end());
  exposed_M.erase(remove_if(exposed_M.begin(), exposed_M.end(), [](Mosquito *h){return(!h->isExposed());}), exposed_M.end());

  infected_H.erase(remove_if(infected_H.begin(), infected_H.end(), [](Human *h){return(!h->isInfected());}), infected_H.end());
  infected_M.erase(remove_if(infected_M.begin(), infected_M.end(), [](Mosquito *h){return(!h->isInfected());}), infected_M.end());

  recovered_H.erase(remove_if(recovered_H.begin(), recovered_H.end(), [](Human *h){return(!h->isRecovered());}), recovered_H.end());
  
  new_infected_H.clear();
  new_exposed_H.clear();
  new_susceptible_H.clear();
  new_recovered_H.clear();

  new_infected_M.clear();
  new_exposed_M.clear();
  new_susceptible_M.clear();
  
  for(vector<Human*>::iterator it = all_dead_H.begin(); it != all_dead_H.end(); ++it){
    delete *it;
  }
  all_dead_H.clear();
  for(vector<Mosquito*>::iterator it = all_dead_M.begin(); it != all_dead_M.end(); ++it){
    delete *it;
  }
  all_dead_M.clear();

}

HostPopulation::~HostPopulation(){
  // Go through and delete all dead hosts
  vector<vector<Human*> > allHumans =  {susceptible_H, exposed_H, infected_H, recovered_H, all_dead_H, new_susceptible_H, new_exposed_H, new_recovered_H, new_infected_H};
  vector<vector<Mosquito*> > allM = {susceptible_M, exposed_M, infected_M, all_dead_M,new_susceptible_M, new_exposed_M, new_infected_M, all_births_M};
  int max = allHumans.size();
  for(int i = 0; i < max; ++i){
    for(vector<Human*>::iterator it = allHumans[i].begin(); it != allHumans[i].end(); ++it){
      delete *it;
    }
  }
  max = allM.size();
  for(int i = 0; i < max; ++i){
    for(vector<Mosquito*>::iterator it = allM[i].begin(); it != allM[i].end(); ++it){
      delete *it;
    }
  }
}

int HostPopulation::countDeadM(){
  return(all_dead_M.size());
}

int HostPopulation::countBirthsM(){
  return(all_births_M.size());
}

void HostPopulation::printStatus(){
  Rcpp::Rcout << "Current day: " << day << endl;
  Rcpp::Rcout << "Total H: " << countH() << endl;
  Rcpp::Rcout << "Susceptible H: " << susceptible_H.size() << endl;
  Rcpp::Rcout << "Exposed H: " << exposed_H.size() << endl;
  Rcpp::Rcout << "Infected H: " << infected_H.size() << endl;
  Rcpp::Rcout << "Recovered H: " << recovered_H.size() << endl << endl;
  Rcpp::Rcout << "Total M: " << countM() << endl;
  Rcpp::Rcout << "Susceptible M: " << susceptible_M.size() << endl;
  Rcpp::Rcout << "Exposed M: " << exposed_M.size() << endl;
  Rcpp::Rcout << "Infected M: " << infected_M.size() << endl;


}
int HostPopulation::getDay(){return(day);}
int HostPopulation::countSusceptibleH(){return(susceptible_H.size());};
int HostPopulation::countExposedH(){return(exposed_H.size());};
int HostPopulation::countInfectedH(){return(infected_H.size());};
int HostPopulation::countRecoveredH(){return(recovered_H.size());};

int HostPopulation::countH(){
  int all = susceptible_H.size() + exposed_H.size() + infected_H.size() + recovered_H.size();
  return(all);
};

int HostPopulation::countSusceptibleM(){return(susceptible_M.size());};
int HostPopulation::countExposedM(){return(exposed_M.size());};
int HostPopulation::countInfectedM(){return(infected_M.size());};

int HostPopulation::countM(){
  int all = susceptible_M.size() + exposed_M.size() + infected_M.size();
  return(all);
};

int HostPopulation::count(AgeClass age, State infState){
  int number = 0;
  switch(infState){
  case Susceptible:
    number = count_if(susceptible_H.begin(), susceptible_H.end(), [age](Human *h){return(h->getAge()==age);});
    break;
  case Exposed:
    number = count_if(exposed_H.begin(), exposed_H.end(), [age](Human *h){return(h->getAge()==age);});
    break;
  case Infected:
    number = count_if(infected_H.begin(), infected_H.end(), [age](Human *h){return(h->getAge()==age);});
    break;
  case Recovered:
    number = count_if(recovered_H.begin(), recovered_H.end(), [age](Human *h){return(h->getAge()==age);});
    break;
  default:
    number = 0;
  }
  return(number);
}

void HostPopulation::incBirths(bool microCeph){
  todayBirths++;
  if(microCeph) todayMicro++;
}

int HostPopulation::getTodayBirths(){
  return(todayBirths);
}


int HostPopulation::getTodayMicroceph(){
  return(todayMicro);
}

double HostPopulation::getMicroprob(int week){
  return(microProbs[week]);
}

double HostPopulation::getRisk(){
  return(risk);
}
