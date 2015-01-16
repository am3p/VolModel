#include <cmath>
#include <exception>

#include <boost/function.hpp>
#include <boost/bind.hpp>
#include <boost/math/distributions.hpp>

#include <ql/Math/Solver/brent.hpp>

#include <Projects/VolSurface/ImpliedVol/BlackScholesVol.hpp>

double BlackScholesVol::BlackScholesPrice(double S, double rf, double q, double Impv){
	boost::math::normal_distribution<> d(0.0, 1.0);
	
	double sig = (double)CP;
	double dp = (std::log(S/Strike) + (rf-q+0.5*Impv*Impv)*TTM)/(Impv*std::sqrt(TTM));
	double dm = dp - Impv*std::sqrt(TTM);
	double res = sig*std::exp(-rf*TTM) * (S*std::exp((rf-q)*TTM)*cdf(d,sig*dp)-Strike*cdf(d,sig*dm));
	
	return res;
}

double BlackScholesVol::BlackScholesImpvProb(double S, double rf, double q, double price, double impv){
	return (this->BlackScholesPrice(S, rf, q, impv)) - price;
}

double BlackScholesVol::BlackScholesImpv(double S, double rf, double q, double price){
	QuantLib::Brent solver;
	boost::function<double (double Impv)> fcn;
	fcn = boost::bind(&BlackScholesVol::BlackScholesImpvProb, this, S, rf, q, price, _1);

	double res;
	try{
		res = solver.solve(fcn,1e-6,0.1,1e-6,1.0);
	}
	catch (std::exception e){
		res = -1;
	}
	
	return res;
}