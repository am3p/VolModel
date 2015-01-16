#include "stdafx.h"

#include <Schedule/Schedule.h>

Schedule::Schedule(boost::shared_ptr<Payoff> payoff, double ttm){
	payoff_ = payoff;
	ttm_ = ttm;

	isBermudan();
}

double Schedule::getttm(){
	return ttm_;
}

boost::shared_ptr<Payoff> Schedule::getPayoff(){
	return payoff_;
}

void Schedule::setttm(double ttm){
	ttm_ = ttm;
}

void Schedule::setPayoff(boost::shared_ptr<Payoff> payoff){
	payoff_ = payoff;
	isBermudan();
}

bool Schedule::isTerminate(std::vector<double> S){
	if (!isBermudan_){
		return true;
	}
	else if(payoff_->BermudanCond(S)){
		return true;
	}
	else{
		return false;
	}
}

void Schedule::isBermudan(){
	boost::shared_ptr<Bermudan> p = boost::dynamic_pointer_cast<Bermudan>(payoff_);
	if (p){
		isBermudan_ = true;
	}
	else{
		isBermudan_ = false;
	}
}

double Schedule::CalculateCF(std::vector<double> S){
	payoff_->CalcPayoff(S);

	std::vector<double> CF = payoff_->getCF();

	std::vector<double>::iterator it = std::min_element(CF.begin(), CF.end());
	return *it;
}