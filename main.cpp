
#include <iostream>
#include <fstream>
#include <vector>

#include <Projects/VolSurface/ImpliedVol/BlackVol.hpp>
#include <Projects/VolSurface/ImpliedVol/BlackScholesVol.hpp>
#include <Projects/VolSurface/Model/SABR.hpp>

#include <boost/math/distributions.hpp>

#include <ql/types.hpp>
#include <ql/Math/Optimization/constraint.hpp>
#include <ql/Math/Optimization/costfunction.hpp>
#include <ql/Math/Optimization/problem.hpp>
#include <ql/Math/Optimization/bfgs.hpp>
#include <ql/Math/Optimization/levenbergmarquardt.hpp>
#include <ql/Math/Interpolations/linearinterpolation.hpp>
#include <ql/Math/Interpolations/cubicinterpolation.hpp>
#include <ql/Math/Interpolations/bilinearinterpolation.hpp>
#include <ql/Math/Interpolations/bicubicsplineinterpolation.hpp>

#include "VBAStruct_impv.h"

class SurfFit : public QuantLib::CostFunction{
public:
	SurfFit(int Size, double F, double t, QuantLib::Array K, QuantLib::Array MarketVol)
		: Size_(Size), F_(F), t_(t), K_(K), MarketVol_(MarketVol) {};
	QuantLib::Disposable<QuantLib::Array> values(const QuantLib::Array& param) const{
		SABR sabr(param);
		QuantLib::Array res(Size_);
		
		for (int i = 0; i < Size_; i++){
			res[i] = (sabr.vol(F_, K_[i], t_) - MarketVol_[i]);
		}
		return res;
	}

	QuantLib::Real value(const QuantLib::Array& param) const{
		QuantLib::Array resArray = values(param);
		QuantLib::Real res = 0;

		for (int i = 0; i < Size_; i++){
			res += resArray[i]*resArray[i];
		}
		return res;
	}

private:
	int Size_;
	double F_;
	double t_;
	QuantLib::Array K_;
	QuantLib::Array MarketVol_;
};


class SurfFit_Adj : public QuantLib::CostFunction{
public:
	SurfFit_Adj(QuantLib::Array smplnum,
				QuantLib::Array F, QuantLib::Array t, QuantLib::Array K, 
				QuantLib::Array MarketVol, 
				int skipind, double nu_a, double nu_b)
		: smplnum_(smplnum), F_(F), t_(t), K_(K), MarketVol_(MarketVol), skipind_(skipind), nu_a_(nu_a), nu_b_(nu_b) {};

	QuantLib::Disposable<QuantLib::Array> values(const QuantLib::Array& param) const{
		QuantLib::Array param_(3,0);
		QuantLib::Array res(MarketVol_.size(), 0);
		int ind = 0;
		for (int i = 0; i < (int)(smplnum_.size()); i++){
			if (smplnum_[i] > 3 && (i+1) > skipind_){
				QuantLib::Real alpha = param[i];
				QuantLib::Real rho = param[smplnum_.size()];
				QuantLib::Real nu = std::exp(nu_a_ * std::log(t_[ind]*365.0) + nu_b_);
				param_[0] = alpha;
				param_[1] = rho;
				param_[2] = nu;
				SABR sabr(param_);
				for (int j = 0; j < (int)(smplnum_[i]); j++){
					QuantLib::Real tmp = sabr.vol(F_[ind], K_[ind], t_[ind]) - MarketVol_[ind];
					res[ind] = tmp;
					ind++;
				}
			}			
			else
			{
				ind += (int)(smplnum_[i]);
			}
		}
		return res;
	}

	QuantLib::Real value(const QuantLib::Array& param) const{
		QuantLib::Array resArray = values(param);
		QuantLib::Real res = 0;
		
		for (int i = 0; i < int(MarketVol_.size()); i++){
			QuantLib::Real tmp = resArray[i];
			res += resArray[i]*resArray[i];
		}
		return res;
	}

private:
	QuantLib::Array smplnum_;
	QuantLib::Array F_;
	QuantLib::Array t_;
	QuantLib::Array K_;
	QuantLib::Array MarketVol_;

	int skipind_;
	double nu_a_; double nu_b_;
};

float rate_interp(struct RateCurve* rc, float t);
float rate_integr(struct RateCurve* rc, float t);


extern "C" __declspec(dllexport) struct SABRparam __stdcall SABRCalc_AllParam(struct SABRData* data){
	QuantLib::Array K(data->size, 0);
	QuantLib::Array MarketVol(data->size, 0);
	for (int i = 0; i < data->size; i++){
		K[i] = (QuantLib::Real)data->opt[i].Strike;
		MarketVol[i] = (QuantLib::Real)data->opt[i].impv;
	}
	QuantLib::Real F = (QuantLib::Real)(data->p.F);
	QuantLib::Real T = (QuantLib::Real)(data->p.TTM / 365.0f);

	QuantLib::EndCriteria endCrit(10000, 1000, 1e-8, 1e-8, 1e-8);

	QuantLib::Array p0_(3, 0);
	p0_[0] = (QuantLib::Real)(data->p.alpha);
	p0_[1] = (QuantLib::Real)(data->p.rho);
	p0_[2] = (QuantLib::Real)(data->p.nu);
	
	QuantLib::Array LB(3,0); LB[0] = 0; LB[1] = -1; LB[2] = 0;
	QuantLib::Array UB(3,0); UB[0] = 1; UB[1] = 1; UB[2] = 100;
	QuantLib::NonhomogeneousBoundaryConstraint bc(LB, UB);

	SurfFit optdata(data->size, F, T, K, MarketVol);

	QuantLib::Problem SABRopt(optdata, bc, p0_);
	QuantLib::LevenbergMarquardt solver;
	QuantLib::EndCriteria::Type solvedCrit = solver.minimize(SABRopt, endCrit);

	QuantLib::Array p1_ = SABRopt.currentValue();

	struct SABRparam res;
	res.alpha = (float)(p1_[0]); res.rho = (float)(p1_[1]); res.nu = (float)(p1_[2]);
	return res;
}

extern "C" __declspec(dllexport) struct SABRparamSeries __stdcall SABRCalc_Adjust(int* ExpSize, int* SkipInd,
																			   struct SABRData* data,
																			   float* nu_a, float* nu_b){
	int Size = 0;
	for (int i = 0; i < (*ExpSize); i++){
		Size += data[i].size;
	}

	int ind = 0;
	QuantLib::Array K(Size, 0);
	QuantLib::Array MarketVol(Size, 0);
	QuantLib::Array T(Size, 0);
	QuantLib::Array F(Size, 0);
	QuantLib::Array p0((*ExpSize) + 1, 0);
	QuantLib::Array LB((*ExpSize) + 1, 0);
	QuantLib::Array UB((*ExpSize) + 1, 0);
	QuantLib::Array smpl_num((*ExpSize), 0);
	for (int i = 0; i < (*ExpSize); i++){
		if (data[i].size > 3 && (i+1) > (*SkipInd)){
			if (data[i].p.alpha > 0)
				p0[i] = (QuantLib::Real)(data[i].p.alpha);
			else
				p0[i] = 0.1;
		}
		LB[i] = 0; UB[i] = 1;
		smpl_num[i] = (QuantLib::Real)(data[i].size);
		for (int j = 0; j < data[i].size; j++){
			K[ind] = (QuantLib::Real)(data[i].opt[j].Strike);
			MarketVol[ind] = (QuantLib::Real)(data[i].opt[j].impv);
			T[ind] = (QuantLib::Real)(data[i].p.TTM / 365.0f);
			F[ind] = (QuantLib::Real)(data[i].p.F);
			ind++;
		}
	}
	LB[(*ExpSize)] = -1; UB[(*ExpSize)] = 1;

	QuantLib::EndCriteria endCrit(10000, 1000, 1e-10, 1e-10, 1e-10);
	QuantLib::NonhomogeneousBoundaryConstraint bc(LB, UB);

	SurfFit_Adj optdata(smpl_num, F, T, K, MarketVol, (*SkipInd), (double)(*nu_a), (double)(*nu_b));

	QuantLib::Problem SABRopt(optdata, bc, p0);
	QuantLib::LevenbergMarquardt solver;
	QuantLib::EndCriteria::Type solvedCrit = solver.minimize(SABRopt, endCrit);

	QuantLib::Array p1 = SABRopt.currentValue();

	struct SABRparamSeries res;
	for (int i = 0; i < (*ExpSize); i++){
		res.p[i].alpha = (float)(p1[i]);
		res.p[i].rho = (float)(p1[(*ExpSize)]);
		res.p[i].nu = (float)(std::exp((*nu_a)*std::log(data[i].p.TTM) + (*nu_b)));
	}
	return res;
}


extern "C" __declspec(dllexport) double __stdcall BlackVolCalc(struct OptData* opt, struct MarketParam* p){
	BlackVol IV(opt->CP, opt->TTM, opt->Strike);
	double res = IV.BlackImpv(p->F, p->rf, opt->price);
	return res;
}

extern "C" __declspec(dllexport) struct SABRparamSeries __stdcall ParamExt(struct SABRparamSeries* p0, int* skipind, 
																		float* nu_a_, float* nu_b_, float* nu_skipind_,
																		float* alpha_first, float* exp_first,
																		float* S0, struct RateCurve* rc, float* q){
	
	// passed p0 should contain parameters whose expiries are longer than 30 days
	// alpha_first: -1 if first market expirys is longer than 5 days
	//				contain positive data when other case
	
	struct SABRparamSeries p1;

	// fixed t: just my purpose
	float t1[10] = {30, 60, 90, 180, 270, 360, 540, 720, 900, 1080};
	p1.Size = 10;

	std::vector<double> nu_list((*skipind)+1);
	std::vector<double> t_list((*skipind)+1);

	for (int i = 0; i < (*skipind); i++){
		nu_list[i] = p0->p[i].nu;
		t_list[i] = p0->t[i];
	}
	nu_list[(*skipind)] = (double)(*nu_skipind_);
	t_list[(*skipind)] = p0->t[(*skipind)];
	QuantLib::MonotonicCubicNaturalSpline nuInterp(t_list.begin(), t_list.end(), nu_list.begin());

	float alpha, rho, nu, F;
	float alpha0, alpha1, rho0, rho1, nu0, nu1;
	float var0, var1, var_diff;

	int t_ind = 0;
	p1.t[0] = t1[0];
	// first expiries (t = 30)
	// if shortest expiry is longer than 30 days
	if (p0->t[t_ind] > 30){
		if ((*alpha_first) > 0){
			var0 = (*exp_first) * (*alpha_first) * (*alpha_first);
			var1 = (p0->t[t_ind]) * (p0->p[t_ind].alpha) * (p0->p[t_ind].alpha);
			var_diff = (var1 - var0) / (p0->t[t_ind] - (*exp_first));
			alpha = std::sqrt((var0 + var_diff * (p1.t[0] - (*exp_first)))/p1.t[0]);
		}
		else{
			var0 = (p0->t[t_ind]) * (p0->p[t_ind].alpha) * (p0->p[t_ind].alpha);
			var1 = (p0->t[t_ind+1]) * (p0->p[t_ind+1].alpha) * (p0->p[t_ind+1].alpha);
			var_diff = (var1 - var0) / (p0->t[t_ind+1] - p0->t[t_ind]);
			alpha = std::sqrt((var0 + var_diff * (p1.t[0] - p0->t[t_ind]))/p1.t[0]);
		}
		rho = p0->p[t_ind].rho;
		nu = p0->p[t_ind].nu;
		F = (*S0) + (p0->p[t_ind].F - (*S0)) / (p0->t[t_ind]) * p1.t[0];
		//F = (*S0) * std::exp((rate_integr(rc, t1[0]) - (*q))*t1[0]/360.0f);
	}
	else{
		var0 = (p0->t[t_ind]) * (p0->p[t_ind].alpha) * (p0->p[t_ind].alpha);
		var1 = (p0->t[t_ind+1]) * (p0->p[t_ind+1].alpha) * (p0->p[t_ind+1].alpha);
		var_diff = (var1 - var0) / (p0->t[t_ind+1] - p0->t[t_ind]);
		alpha = std::sqrt((var0 + var_diff * (p1.t[0] - p0->t[t_ind]))/p1.t[0]);

		rho0 = p0->p[t_ind].rho;
		rho1 = p0->p[t_ind+1].rho;
		rho = rho0 + (rho1-rho0) / ((p0->t[t_ind+1]) - (p0->t[t_ind])) * ((p1.t[0]) - (p0->t[t_ind]));

		nu = (float)nuInterp(p1.t[0]);
		F = (*S0) + (p0->p[t_ind].F - (*S0)) / p0->t[t_ind] * p1.t[0];
		//F = (*S0) * std::exp((rate_integr(rc, t1[0]) - (*q))*t1[0]/360.0f);
	}

	p1.p[0].alpha = alpha;
	p1.p[0].rho = rho;
	p1.p[0].nu = nu;
	p1.p[0].F = F;
	p1.p[0].rf = rate_interp(rc, p1.t[0]);
	//p1.p[0].rf = rate_integr(rc, p1.t[0]);

	// rest of expiries

	for (int i = 1; i < 10; i++){
		p1.t[i] = t1[i];
		// find interval
		if (p1.t[i] > p0->t[t_ind]){
			while (p1.t[i] > p0->t[t_ind] && t_ind < p0->Size){
				t_ind++;
			}
		}
		
		// no-standard part
		if (t1[i] > 360){
			F = (*S0) * std::exp((rate_integr(rc, t1[i]) - (*q))*t1[i]/360.0f);
		}
		else{
			F = p0->p[t_ind-1].F + (p0->p[t_ind].F - p0->p[t_ind-1].F)/(p0->t[t_ind] - p0->t[t_ind-1]) * (p1.t[i] - p0->t[t_ind-1]);
		}

		// larger than market available expiry
		if (t_ind == p0->Size){
			var0 = (p0->t[t_ind-2]) * (p0->p[t_ind-2].alpha) * (p0->p[t_ind-2].alpha);
			var1 = (p0->t[t_ind-1]) * (p0->p[t_ind-1].alpha) * (p0->p[t_ind-1].alpha);
			var_diff = (var1 - var0) / (p0->t[t_ind-1] - p0->t[t_ind-2]);
			alpha = std::sqrt((var0 + var_diff * (p1.t[i] - p0->t[t_ind-2]))/p1.t[i]);

			rho = p0->p[t_ind-1].rho;
			nu = std::exp((*nu_a_)*std::log(p1.t[i]) + (*nu_b_));
		}
	
		// in-between market avilable expiry
		else{
			var0 = (p0->t[t_ind-1]) * (p0->p[t_ind-1].alpha) * (p0->p[t_ind-1].alpha);
			var1 = (p0->t[t_ind]) * (p0->p[t_ind].alpha) * (p0->p[t_ind].alpha);
			var_diff = (var1 - var0) / (p0->t[t_ind] - p0->t[t_ind-1]);
			alpha = std::sqrt((var0 + var_diff * (p1.t[i] - p0->t[t_ind-1]))/p1.t[i]);

			rho0 = p0->p[t_ind-1].rho;
			rho1 = p0->p[t_ind].rho;
			rho = rho0 + (rho1-rho0) / ((p0->t[t_ind]) - (p0->t[t_ind-1])) * ((p1.t[i]) - (p0->t[t_ind-1]));

			// border of raw data and trend-incorporating exponential estimation
			if((t_ind+1) == (*skipind)){ 
				// estimation part
				nu = std::exp((*nu_a_)*std::log(p1.t[i]) + (*nu_b_)) * (p1.t[i] - p0->t[t_ind-1]) / (p0->t[t_ind] - p0->t[t_ind-1]);
				// monotonic cubic interpolation part
				nu += (float)nuInterp(p1.t[i]) * (p0->t[t_ind]-p1.t[i]) / (p0->t[t_ind] - p0->t[t_ind-1]);
			}
			else if (t_ind+1 > (*skipind)){
				nu = std::exp((*nu_a_)*std::log(p1.t[i]) + (*nu_b_));
			}
			else{
				// monotonic cubic interpolation
				nu = (float)nuInterp(p1.t[i]);
			}
		}

		p1.p[i].alpha = alpha;
		p1.p[i].rho = rho;
		p1.p[i].nu = nu;
		p1.p[i].F = F;
		p1.p[i].rf = rate_interp(rc, p1.t[i]);
		//p1.p[i].rf = rate_integr(rc, p1.t[i]);
	}

	return p1;
}

extern "C" __declspec(dllexport) struct VolSurf __stdcall LocalVolCalcSABR(struct SABRparamSeries* p, float *S0, float* q){
	float alpha, alpha_p, alpha_m, rho, rho_p, rho_m, nu, nu_p, nu_m, F, F_p, F_m, t, t_p, t_m, r;
	float t_eps = 0.001, K_eps;
	int Ksize = 17, tsize = p->Size;
	
	struct VolSurf vs;
	for (int i = 0; i < Ksize; i++){
		vs.K[i] = (float)((*S0) * (0.4 + (float)(i) * 0.05));
	}

	for (int i = 0; i < p->Size; i++){
		alpha = p->p[i].alpha;
		rho = p->p[i].rho;
		nu = p->p[i].nu;
		F = p->p[i].F;
		t = p->t[i]/360.0f;
		r = p->p[i].rf;

		vs.t[i] = p->t[i];		
		if (i == p->Size - 1){
			alpha_p = alpha + (alpha - (p->p[i-1].alpha))/(t - p->t[i-1]/360.0f) * t_eps;
			rho_p = rho + (rho - (p->p[i-1].rho))/(t - p->t[i-1]/360.0f) * t_eps;
			nu_p = nu + (nu - (p->p[i-1].nu))/(t - p->t[i-1]/360.0f) * t_eps;
			F_p = F + (F - (p->p[i-1].F))/(t - p->t[i-1]/360.0f) * t_eps;
		}
		else{
			alpha_p = alpha + ((p->p[i+1].alpha) - alpha)/(p->t[i+1]/360.0f - t) * t_eps;
			rho_p = rho + ((p->p[i+1].rho) - rho)/(p->t[i+1]/360.0f - t) * t_eps;
			nu_p = nu + ((p->p[i+1].nu) - nu)/(p->t[i+1]/360.0f - t) * t_eps;
			F_p = F + ((p->p[i+1].F) - F)/(p->t[i+1]/360.0f - t) * t_eps;
		}
		t_p = t + t_eps;

		if (i == 0){
			alpha_m = alpha + ((p->p[i+1].alpha) - alpha)/(p->t[i+1]/360.0f - t) * -t_eps;
			rho_m = rho + ((p->p[i+1].rho) - rho)/(p->t[i+1]/360.0f - t) * -t_eps;
			nu_m = nu + ((p->p[i+1].nu) - nu)/(p->t[i+1]/360.0f - t) * -t_eps;
			F_m = F + ((p->p[i].F) - (*S0))/(p->t[i]/360.0f) * -t_eps;
		}
		else{
			alpha_m = alpha + (alpha - (p->p[i-1].alpha))/(t - p->t[i-1]/360.0f) * -t_eps;
			rho_m = rho + (rho - (p->p[i-1].rho))/(t - p->t[i-1]/360.0f) * -t_eps;
			nu_m = nu + (nu - (p->p[i-1].nu))/(t - p->t[i-1]/360.0f) * -t_eps;
			F_m = F + (F - (p->p[i-1].F))/(t - p->t[i-1]/360.0f) * -t_eps;
		}
		t_m = t - t_eps;

		float K, K_p, K_m;
		float impv, impv_tp, impv_tm, impv_Kp, impv_Km;
		float w, w_tp, w_tm, w_Kp, w_Km;
		float locv, locv_nom, locv_denom;
		float y, f, fs, fl, df, d2f, d1;
		float dwdy, d2wdy2, dwdt;
		float dvdK, d2vdK2, dvdt;

		for (int j = 0; j < Ksize; j++){			
			K = vs.K[j];
			y = std::log(F/K);

			if (std::abs(1/y) > 0.01)
				K_eps = std::abs(1.0/y) * 0.01;
			else
				K_eps = 0.001;
			
			y = nu/alpha * std::log(F/K);
			if (std::abs(y) <= 0.01){
				f = (float)(1 - 0.5 * rho * y + (-0.25 * rho * rho + 1.0/6.0) * y * y);
				df = (float)(-0.5 * rho + 2.0 * (-0.25 * rho * rho + 1.0/6.0) * y - 0.125 * (6.0 * rho * rho - 5.0) * rho * y * y);
				d2f = (float)(2.0 * (-0.25 * rho * rho + 1.0/6.0) - 0.25 * (6.0 * rho * rho - 5.0) * rho * y + 
							  12.0 * (-5.0/16.0 * pow(rho, 4) + pow(rho, 2)/3.0 - 17.0/360.0) * y * y);
			}
			else{
				fs = std::sqrt(1 - 2*rho*y + y*y);
				fl = std::log((fs-rho+y)/(1-rho));
				f = y/fl;
				df = (fl*fs-y)/(fl*fl*fs);
				d2f = (fl*(3.0f*rho*y-y*y-2.0f)+2.0f*fs*y)/pow(fl*fs,3.0f);					
			}
			impv = (float)(alpha * f * (1+(0.25*rho*nu*alpha+(2-3*rho*rho)/24.0*nu*nu)*t));
			w = impv * impv * t;
			dvdK = -nu * df / K * (1+(0.25*rho*nu*alpha+(2-3*rho*rho)/24.0*nu*nu)*t);
			d2vdK2 = (nu * df / K / K + nu * nu * d2f / alpha / K / K) * (1+(0.25*rho*nu*alpha+(2-3*rho*rho)/24.0*nu*nu)*t);
			d1 = (std::log((*S0)/K) + (r-(*q)-0.5*impv*impv)*t)/(impv*std::sqrt(t));

			K_p = K * std::exp(K_eps);
			y = nu/alpha * std::log(F/K_p);
			if (std::abs(y) <= 0.01){
				f = (float)(1 - 0.5 * rho * y + (-0.25 * rho * rho + 1.0/6.0) * y * y);
			}
			else{
				fs = std::sqrt(1 - 2*rho*y + y*y);
				fl = std::log((fs-rho+y)/(1-rho));
				f = y/fl;		
			}
			impv_Kp = (float)(alpha * f * (1.0+(0.25*rho*nu*alpha+(2.0-3.0*rho*rho)/24.0*nu*nu)*t));
			w_Kp = impv_Kp * impv_Kp * t;

			K_m = K * std::exp(-K_eps);
			y = nu/alpha * std::log(F/K_m);
			if (std::abs(y) <= 0.01){
				f = (float)(1 - 0.5 * rho * y + (-0.25 * rho * rho + 1.0/6.0) * y * y);
			}
			else{
				fs = std::sqrt(1 - 2*rho*y + y*y);
				fl = std::log((fs-rho+y)/(1-rho));
				f = y/fl;		
			}
			impv_Km = (float)(alpha * f * (1.0+(0.25*rho*nu*alpha+(2.0-3.0*rho*rho)/24.0*nu*nu)*t));
			w_Km = impv_Km * impv_Km * t;

			y = nu_p/alpha_p * std::log(F_p/K);
			if (std::abs(y) <= 0.01){
				f = (float)(1 - 0.5 * rho_p * y + (-0.25 * rho_p * rho_p + 1.0/6.0) * y * y);
			}
			else{
				fs = std::sqrt(1 - 2*rho_p*y + y*y);
				fl = std::log((fs-rho_p+y)/(1-rho_p));
				f = y/fl;		
			}
			impv_tp = (float)(alpha_p * f * (1.0+(0.25*rho_p*nu_p*alpha_p+(2.0-3.0*rho_p*rho_p)/24.0*nu_p*nu_p)*t_p));
			w_tp = impv_tp * impv_tp * t_p;
			
			y = nu_m/alpha_m * std::log(F_m/K);
			if (std::abs(y) <= 0.01){
				f = (float)(1 - 0.5 * rho_m * y + (-0.25 * rho_m * rho_m + 1.0/6.0) * y * y);
			}
			else{
				fs = std::sqrt(1 - 2*rho_m*y + y*y);
				fl = std::log((fs-rho_m+y)/(1-rho_m));
				f = y/fl;		
			}
			impv_tm = (float)(alpha_m * f * (1.0+(0.25*rho_m*nu_m*alpha_m+(2.0-3.0*rho_m*rho_m)/24.0*nu_m*nu_m)*t_m));
			w_tm = impv_tm * impv_tm * t_m;

			y = std::log(K/F);
			if (w_tp - w_tm > 0)
				dwdt = (w_tp - w_tm) / (2*t_eps);
			else if (w_tp - w > 0)
				dwdt = (w_tp - w) / t_eps;
			else if (w - w_tm > 0)
				dwdt = (w - w_tm) / t_eps;

			dvdt = (impv_tp * impv_tp - impv * impv) / (2 * t_eps * impv);
			
			dwdy = (w_Kp - w_Km) / (2*K_eps);
			d2wdy2 = (w_Kp - 2*w + w_Km) / (K_eps * K_eps);
			//locv_nom = (float)(dwdt);
			//locv_denom = (float)(1.0 - y/w*dwdy + 0.25*(-0.25 - 1.0/w + y*y/w/w)*dwdy*dwdy + 0.5*d2wdy2);
			locv_nom = impv*impv + 2*impv*t*dvdt + 2*(r-(*q))*impv*K*t*dvdK;
			locv_denom = pow(1.0f+d1*K*std::sqrt(t)*dvdK, 2.0f) + impv*K*K*t*(d2vdK2 - d1*std::sqrt(t)*dvdK*dvdK);
			locv = std::sqrt(locv_nom / locv_denom);

			vs.vol[j + i*Ksize] = locv;
		}
	}

	return vs;
}

extern "C" __declspec(dllexport) struct VolSurf __stdcall LocalVolCalcVar(struct VolSurf* vs, float *F0, float* S0){
	float F, F_p, F_m, t, t_p, t_m;
	float t_eps = 0.0001, K_eps;
	int Ksize = 17, tsize = 9;		// Hard-coded!!
	
	std::vector<float> Ks(Ksize);
	std::vector<float> ts(tsize);
	QuantLib::Matrix IVs(tsize, Ksize);

	for (int i = 0; i < tsize; i++){
		ts[i] = ((float)(vs->t[i]))/360.0f;
		for (int j = 0; j < Ksize; j++){
			Ks[j] = (float)((*S0) * (0.4 + (float)(j) * 0.05));
			IVs[i][j] = (QuantLib::Real)(vs->vol[i*Ksize + j]);
		}
	}

	QuantLib::BicubicSpline volspline(Ks.begin(), Ks.end(), ts.begin(), ts.end(), IVs);
	struct VolSurf vs_res;

	for (int i = 0; i < tsize; i++){
		F = F0[i];
		t = vs->t[i]/360.0f;

		if (i == tsize - 1){
			F_p = F + (F - F0[i-1])/(t - vs->t[i-1]/360.0f) * t_eps;
		}
		else{
			F_p = F + (F0[i+1] - F)/(vs->t[i+1]/360.0f - t) * t_eps;
		}
		t_p = t + t_eps;

		if (i == 0){
			F_m = F + (F - (*S0))/(vs->t[i]/360.0f) * -t_eps;
		}
		else{
			F_m = F + (F - F0[i-1])/(t - vs->t[i-1]/360.0f) * -t_eps;
		}
		t_m = t - t_eps;

		float K, K_p, K_m, y;
		float impv, impv_tp, impv_tm, impv_Kp, impv_Km;
		float w, w_tp, w_tm, w_Kp, w_Km;
		float locv, locv_nom, locv_denom;
		float dwdy, d2wdy2, dwdt;

		for (int j = 0; j < Ksize; j++){			
			K = Ks[j];
			y = std::log(K/F);
			if (std::abs(y) > 0.001)
				K_eps = y * 0.0001;
			else
				K_eps = 1e-6;
			
			impv = (float)(IVs[i][j]);
			w = impv * impv * t;

			K_p = K * std::exp(K_eps);
			impv_Kp = (float)(volspline(K_p, t, true));
			w_Kp = impv_Kp * impv_Kp * t;

			K_m = K * std::exp(-K_eps);
			impv_Km = (float)(volspline(K_m, t, true));
			w_Km = impv_Km * impv_Km * t;

			impv_tp = (float)(volspline(K, t_p, true));
			w_tp = impv_tp * impv_tp * t_p;
			
			impv_tm = (float)(volspline(K, t_m, true));
			w_tm = impv_tm * impv_tm * t_m;

			y = std::log(K/F);
			if (w_tp - w_tm > 0)
				dwdt = (w_tp - w_tm) / (2*t_eps);
			else if (w_tp - w > 0)
				dwdt = (w_tp - w) / t_eps;
			else if (w - w_tm > 0)
				dwdt = (w - w_tm) / t_eps;
			
			dwdy = (w_Kp - w_Km) / (2*K_eps);
			d2wdy2 = (w_Kp - 2*w + w_Km) / (K_eps * K_eps);
			locv_nom = (float)(dwdt);
			locv_denom = (float)(1.0 - y/w*dwdy + 0.25*(-0.25 - 1.0/w + y*y/w/w)*dwdy*dwdy + 0.5*d2wdy2);
			locv = std::sqrt(locv_nom / locv_denom);

			vs_res.vol[j + i*Ksize] = locv;
		}
	}

	return vs_res;
}



float rate_interp(struct RateCurve* rc, float t){
	int i = 0;
	float res;
	while (t > rc->t[i]){
		i++;
		if (i == rc->size)
			break;
	}
	if (i == 0){
		res = rc->r[0];
	}
	else if (i == rc->size){
		res = rc->r[(rc->size)-1];
	}
	else{
		res = rc->r[i-1] + (rc->r[i]-rc->r[i-1])/(rc->t[i]-rc->t[i-1])*(t-rc->t[i-1]);
	}
	return res;
}

float rate_integr(struct RateCurve* rc, float t){
	int i = 0;
	float res = 0, r_prior = rc->r[0], t_prior = 0;
	while (t > rc->t[i]){
		res += (r_prior + rc->r[i])/2.0f * (rc->t[i]-t_prior);
		r_prior = rc->r[i]; t_prior = rc->t[i];
		i++;
		if (i == rc->size)
			break;
	}
	if (i == rc->size){
		res += (r_prior + rc->r[i])/2.0f * (t-t_prior);
	}
	else{
		res += (rate_interp(rc, t)+r_prior)/2.0f * (t-t_prior);
	}
	return res / t;
}

//int main(){
//	QuantLib::Array t(3);
//	QuantLib::Array K(4);
//	for (int i = 0; i < 3; i++){
//		t[i] = (QuantLib::Real)(i+1);
//	}
//	for (int i = 0; i < 4; i++){
//		K[i] = (QuantLib::Real)(0.8 + i * 0.1);
//	}
//	QuantLib::Matrix v(3, 4);
//	for (int i = 0; i < 3; i++){
//		for (int j = 0; j < 4; j++){
//			v[i][j] = (QuantLib::Real)(i * i * 1.0 + j * j * 0.1);
//		}
//	}
//
//	for (int i = 0; i < 3; i++){
//		for (int j = 0; j < 4; j++){
//			std::cout << t[i] << " " << K[j] << " " << v[i][j] << std::endl;
//		}
//	}
//
//	QuantLib::BicubicSpline test(K.begin(), K.end(), t.begin(), t.end(), v);
//	std::cout << test.isInRange(1.05, 2) << std::endl;
//	if (test.isInRange(1.05, 2))
//		std::cout << test(1.05, 2) << std::endl;
//
//	return 0;
//}