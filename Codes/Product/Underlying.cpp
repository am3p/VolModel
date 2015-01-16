#include "stdafx.h"

#include <Product/Underlying.h>

// Underlying Class
Underlying::Underlying(const double S, 
					   const boost::shared_ptr<InterestRate> rf,
					   const boost::shared_ptr<DividendRate> div,
					   const boost::shared_ptr<Volatility> vol){
	S_ = S;
	rf_ = rf;
	div_ = div;
	vol_ = vol;

	chkInterestType();
	chkDividendType();
	chkVolatilityType();
}

double Underlying::getStock(){
	return S_;
}

double Underlying::getInterest(double t){
	double rf;
	if (!isrfterm){
		rf = rf_->getInterest();
	}
	else{
		rf = rf_->getInterest(t);
	}
	return rf;
}

void Underlying::chkInterestType(){
	boost::shared_ptr<InterestConst> rf_pt = boost::dynamic_pointer_cast<InterestConst>(rf_);
	if (rf_pt){
		isrfterm = false;
	}
	boost::shared_ptr<InterestTerm> rf_term = boost::dynamic_pointer_cast<InterestTerm>(rf_);
	if (rf_term){
		isrfterm = true;
	}
}

double Underlying::getDividend(double t){
	double q;
	if (!isdivterm){
		q = div_->getDividend();
	}
	else{
		q = div_->getDividend(t);
	}
	return q;
}

void Underlying::chkDividendType(){
	boost::shared_ptr<DividendConst> q_pt = boost::dynamic_pointer_cast<DividendConst>(div_);
	if (q_pt){
		isdivterm = false;
		q_pt.reset();
	}
	else{
		boost::shared_ptr<DividendTerm> q_term = boost::dynamic_pointer_cast<DividendTerm>(div_);
		if (q_term){
			isdivterm = true;
			q_term.reset();
		}
	}
}

double Underlying::getVolatility(double t, double S){
	double v;
	if (!isvolterm && !isvolsurf){
		v = vol_->getVolatility();
	}
	else{
		if(isvolterm){
			v = vol_->getVolatility(t);
		}
		else if(isvolsurf){
			v = vol_->getVolatility(t, S);
		}
	}
	return v;
}

void Underlying::chkVolatilityType(){
	boost::shared_ptr<VolatilityConst> vol_pt = boost::dynamic_pointer_cast<VolatilityConst>(vol_);
	if (vol_pt){
		isvolterm = false;
		isvolsurf = false;
		vol_pt.reset();
	}
	else{
		boost::shared_ptr<VolatilityTerm> vol_term = boost::dynamic_pointer_cast<VolatilityTerm>(vol_);
		if (vol_term){
			isvolterm = true;
			isvolsurf = false;
			vol_term.reset();
		}
		else{
			boost::shared_ptr<VolatilitySurf> vol_surf = boost::dynamic_pointer_cast<VolatilitySurf>(vol_);
			if (vol_surf){
				isvolterm = false;
				isvolsurf = true;
				vol_surf.reset();
			}
		}
	}
}