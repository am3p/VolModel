#include <cmath>

#include <Projects/VolSurface/Model/SABR.hpp>

double SABR::vol (double F, double K, double t){
	double z, xz;

	if (F!=K){
		z = nu/alpha * std::log(F/K);
		xz = std::log((std::sqrt(1-2*rho*z+pow(z,2))+z-rho)/(1-rho));
	}
	else{
		z = 1; xz = 1;
	}

	double impv = alpha * z/xz * (1+(0.25*rho*nu*alpha+(2-3*pow(rho,2))/24*pow(nu,2))*t);
	return impv;	
}