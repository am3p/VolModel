#include <cmath>
#include <exception>

#include <boost/function.hpp>
#include <boost/bind.hpp>
#include <boost/math/distributions.hpp>

#include <ql/Math/Solver/brent.hpp>

#include <Projects/VolSurface/ImpliedVol/BlackVol.hpp>

double BlackVol::BlackPrice(double F, double rf, double Impv){
	boost::math::normal_distribution<> d(0.0, 1.0);
	
	double sig = (double)CP;
	double dp = (std::log(F/Strike) + 0.5*Impv*Impv*TTM)/(Impv*std::sqrt(TTM));
	double dm = dp - Impv*std::sqrt(TTM);
	double res = sig*std::exp(-rf*TTM) * (F*cdf(d,sig*dp)-Strike*cdf(d,sig*dm));
	
	return res;
}

double BlackVol::BlackImpvProb(double F, double rf, double price, double impv){
	return (this->BlackPrice(F, rf, impv)) - price;
}

double BlackVol::BlackImpv(double F, double rf, double price){
	QuantLib::Brent solver;
	boost::function<double (double Impv)> fcn;
	fcn = boost::bind(&BlackVol::BlackImpvProb, this, F, rf, price, _1);

	double res;
	try{
		res = solver.solve(fcn,1e-6,0.1,1e-6,1.0);
	}
	catch (std::exception e){
		res = -1;
	}
	
	return res;
}