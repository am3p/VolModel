#include "stdafx.h"

#include <Product/Payoff.h>

// Payoff Class
std::vector<double> Payoff::getCF(){
	return CF;
}

// Bermudan Case
// Vanilla type
Vanilla::Vanilla(std::vector<double> strike, double sign, double particip, bool CP){
	strike_ = strike;
	sign_ = sign;
	particip_ = particip;
	CP_ = CP;

	for(unsigned int i = 0; i < strike_.size(); i++){
		CF.push_back(0.0);
	}
}

void Vanilla::CalcPayoff(std::vector<double> S){
	if(CP_){
		for(unsigned int i = 0; i < S.size(); i++){
			CF[i] = (S[i] >= strike_[i]) ? sign_*particip_*(S[i] - strike_[i]) : 0;
		}
	}
	else{
		for(unsigned int i = 0; i < S.size(); i++){
			CF[i] = (S[i] <= strike_[i]) ? sign_*particip_*(strike_[i] - S[i]) : 0;
		}
	}
}		

bool Vanilla::BermudanCond(std::vector<double> S){
	return false;
}


// Digital type
Digital::Digital(std::vector<double> strike, double sign, bool CP, double Coupon){
	strike_ = strike;
	sign_ = sign;
	CP_ = CP;
	Coupon_ = Coupon;

	for(unsigned int i = 0; i < strike_.size(); i++){
		CF.push_back(0.0);
	}
}

void Digital::CalcPayoff(std::vector<double> S){
	if(CP_){
		for(unsigned int i = 0; i < S.size(); i++){
			CF[i] = (S[i] >= strike_[i]) ? Coupon_*sign_ : 0;
		}
	}
	else{
		for(unsigned int i = 0; i < S.size(); i++){
			CF[i] = (S[i] <= strike_[i]) ? Coupon_*sign_ : 0;
		}
	}
}

bool Digital::BermudanCond(std::vector<double> S){
	bool result = true;
	if(CP_){
		for (unsigned int i = 0; i < S.size(); i++){
			result = (result && (S[i] >= strike_[i]));
		}
	}
	else{
		for (unsigned int i = 0; i < S.size(); i++){
			result = (result && (S[i] <= strike_[i]));
		}
	}
	return result;
}


// Terminal Case
// Digital Call + Put case (most common)
DigiCall_Put::DigiCall_Put(std::vector<double> strike, double sign, bool CP, 
						   std::vector<double> strike_put, double Coupon, double Coupon_Dummy){
	strike_ = strike;
	sign_ = sign;
	CP_ = CP;
	strike_put_ = strike_put;
	Coupon_ = Coupon;
	Coupon_Dummy_ = Coupon_Dummy;

	for(unsigned int i = 0; i < strike_.size(); i++){
		CF.push_back(0.0);
	}
}

bool DigiCall_Put::TerminalCond(std::vector<double> MinPrice, std::vector<double> DownBarrier){
	bool result = true;
	for(unsigned int i = 0; i < MinPrice.size(); i++){
		result = result && (MinPrice[i] <= DownBarrier[i]);
	}
	return result;
}

void DigiCall_Put::CalcPayoff(std::vector<double> S){
	bool result = true;
	if(!DownBTouched){
		for (unsigned int i = 0; i < S.size(); i++){
			result = result && (S[i] >= strike_[i]);
		}
		if (result){
			for (unsigned int i = 0; i < S.size(); i++){
				CF[i] = Coupon_*sign_;
			}
		}
		else{
			for (unsigned int i = 0; i < S.size(); i++){
				CF[i] = Coupon_Dummy_*sign_;
			}
		}			
	}
	else{
		for (unsigned int i = 0; i < S.size(); i++){
			result = result && (S[i] >= strike_put_[i]);
		}
		if (result){
			for (unsigned int i = 0; i < S.size(); i++){
				CF[i] = Coupon_*sign_;
			}
		}
		else{
			for (unsigned int i = 0; i < S.size(); i++){
				CF[i] = (S[i] - strike_put_[i])*sign_;
			}
		}
	}

}
