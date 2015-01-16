#include "stdafx.h"

#include <interpolation.h>

#include <Product/Underlying.h>

// Interest Rate Class
InterestConst::InterestConst(double r){
	r_ = r;
}

double InterestConst::getInterest(){
	return r_;
}

void InterestConst::setInterest(double r){
	r_ = r;
}

InterestTerm::InterestTerm(alglib::real_1d_array t, alglib::real_1d_array r){
	t_ = t;
	r_ = r;

	alglib::spline1dbuildlinear(t_, r_, s);
}

double InterestTerm::getInterest(double t){
	return alglib::spline1dcalc(s, t);
}

void InterestTerm::setInterest(alglib::real_1d_array t, alglib::real_1d_array r){
	t_ = t;
	r_ = r;
}

// Dividend Rate Class
DividendConst::DividendConst(double q){
	q_ = q;
}

double DividendConst::getDividend(){
	return q_;
}

void DividendConst::setDividend(double q){
	q_ = q;
}

DividendTerm::DividendTerm(alglib::real_1d_array t, alglib::real_1d_array q){
	t_ = t;
	q_ = q;

	alglib::spline1dbuildlinear(t_, q_, s);
}

double DividendTerm::getDividend(double t){	
	return alglib::spline1dcalc(s, t);
}

void DividendTerm::setDividend(alglib::real_1d_array t, alglib::real_1d_array q){
	t_ = t;
	q_ = q;
}

DividendDiscrete::DividendDiscrete(alglib::real_1d_array t, alglib::real_1d_array q){
	t_ = t;
	q_ = q;
}

// Volatility Class
VolatilityConst::VolatilityConst(double v){
	v_ = v;
}

double VolatilityConst::getVolatility(){
	return v_;
}

void VolatilityConst::setVolatility(double v){
	v_ = v;
}

VolatilityTerm::VolatilityTerm(alglib::real_1d_array t, alglib::real_1d_array v){
	t_ = t;
	v_ = v;

	alglib::spline1dbuildlinear(t_, v_, s);
}

double VolatilityTerm::getVolatility(double t){
	return alglib::spline1dcalc(s, t);
}

void VolatilityTerm::setVolatility(alglib::real_1d_array t, alglib::real_1d_array v){
	t_ = t;
	v_ = v;
}

VolatilitySurf::VolatilitySurf(alglib::real_1d_array t, alglib::real_1d_array K, alglib::real_1d_array v){
	t_ = t;
	K_ = K;
	v_ = v;

	alglib::spline2dbuildbilinearv(t_, t_.length(), K_, K_.length(), v_, 1, s);
}

double VolatilitySurf::getVolatility(double t, double K){	
	return alglib::spline2dcalc(s, t, K);
}

void VolatilitySurf::setVolatility(alglib::real_1d_array t, alglib::real_1d_array K, alglib::real_1d_array v){
	t_ = t;
	K_ = K;
	v_ = v;
}