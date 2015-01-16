#include "stdafx.h"

#include "BSProcessMC.h"

// MC_test.cpp : 콘솔 응용 프로그램에 대한 진입점을 정의합니다.
//

int main()
{
	time_t t0, t1;
	time(&t0);

	double d_t[] = {1.0/12.0, 1.0/6.0, 0.25, 0.5, 0.75, 1, 1.5, 2, 2.5, 3};
	double d_K[] = {60, 65, 70, 75, 80, 85, 90, 95, 100, 105, 110, 115, 120, 125, 130, 135, 140};
	double d_v[] = {0.464279503, 0.337724603, 0.323311974, 0.303414695, 0.407712061, 0.404356426, 0.425828078, 0.428642032, 0.279895579, 0.258287935, 
					0.49118271, 0.332024749, 0.319894195, 0.315473015, 0.39952287, 0.385189883, 0.369974691, 0.358038613, 0.237488881, 0.221990926, 
					0.527715085, 0.322719871, 0.305939245, 0.322057061, 0.368278349, 0.333163925, 0.294739983, 0.289372273, 0.203607241, 0.194544425, 
					0.577517724, 0.300239585, 0.279043341, 0.31177296, 0.311032813, 0.266535909, 0.233713317, 0.237521094, 0.179745056, 0.175295545, 
					0.635228258, 0.249579116, 0.262935024, 0.26312606, 0.25412267, 0.217738379, 0.195110315, 0.202736528, 0.164027958, 0.162337112, 
					0.578519191, 0.214129949, 0.241642594, 0.203798307, 0.212473559, 0.191599566, 0.172967025, 0.180639608, 0.154135734, 0.15389387, 
					0.367291203, 0.174528026, 0.178454393, 0.167089393, 0.184206129, 0.17972115, 0.161024098, 0.167290967, 0.148278457, 0.148656949, 
					0.167286606, 0.129317979, 0.141541141, 0.147629222, 0.165950629, 0.175022703, 0.155228411, 0.159956098, 0.14524848, 0.14572975, 
					0.103795111, 0.111566853, 0.129066232, 0.138276316, 0.15536632, 0.173655859, 0.153160977, 0.156803699, 0.14424134, 0.144497134, 
					0.115502172, 0.105933776, 0.127127207, 0.135959698, 0.15068711, 0.173859974, 0.153430452, 0.156622866, 0.144732503, 0.144539061, 
					0.182390922, 0.115823124, 0.132098843, 0.141029144, 0.15127563, 0.17496038, 0.15516687, 0.158536838, 0.146399863, 0.14562661, 
					0.35323843, 0.161187037, 0.149852376, 0.157176166, 0.158007645, 0.176910636, 0.157924948, 0.161878024, 0.148903727, 0.147455462, 
					0.500288143, 0.241844574, 0.20047703, 0.191018755, 0.173563841, 0.179984456, 0.161402356, 0.166143739, 0.152089721, 0.149929158, 
					0.506102065, 0.295347318, 0.260788263, 0.245213983, 0.202419067, 0.184904121, 0.165503621, 0.170919017, 0.155729728, 0.152870256, 
					0.490244228, 0.322114712, 0.296445441, 0.298450592, 0.247225399, 0.19295586, 0.170281874, 0.175950843, 0.15967593, 0.156149269, 
					0.476176673, 0.331465596, 0.321699707, 0.324609883, 0.299370127, 0.205876926, 0.175914559, 0.181121738, 0.163830301, 0.159668037, 
					0.464849529, 0.334575963, 0.338798318, 0.330204558, 0.340288187, 0.224827138, 0.182686569, 0.186442636, 0.168153704, 0.163359945};

	alglib::real_1d_array t;
	alglib::real_1d_array K;
	alglib::real_1d_array v;

	t.setcontent(10, d_t);
	K.setcontent(17, d_K);
	v.setcontent(170, d_v);

	boost::shared_ptr<InterestConst> rf1(new InterestConst(0.02));
	boost::shared_ptr<DividendConst> div1(new DividendConst(0.01));
	//boost::shared_ptr<VolatilityConst> vol1(new VolatilityConst(0.15));
	boost::shared_ptr<VolatilitySurf> vol1(new VolatilitySurf(t, K, v));
	double S1 = 100;
	Underlying Stock1(S1, rf1, div1, vol1);	

	boost::shared_ptr<InterestConst> rf2(new InterestConst(0.02));
	boost::shared_ptr<DividendConst> div2(new DividendConst(0.01));
	//boost::shared_ptr<VolatilityConst> vol2(new VolatilityConst(0.25));
	boost::shared_ptr<VolatilitySurf> vol2(new VolatilitySurf(t, K, v));
	double S2 = 100;
	Underlying Stock2(S2, rf2, div2, vol2);	

	std::vector<Underlying> Stock;
	Stock.push_back(Stock1);
	Stock.push_back(Stock2);

	std::vector<Schedule> ProdSched;
	Schedule sched_atom;
	
	double d_Strike1[] = {100, 100};
	std::vector<double> Strike1(d_Strike1, d_Strike1+2);
	boost::shared_ptr<Payoff> payoff1(new Digital(Strike1, 1, true, 0.5));
	sched_atom.setttm(0.5);
	sched_atom.setPayoff(payoff1);
	ProdSched.push_back(sched_atom);

	double d_Strike2[] = {100, 100};
	std::vector<double> Strike2(d_Strike2, d_Strike2+2);
	double d_Strike_Put[] = {80, 80};
	std::vector<double> Strike_Put(d_Strike_Put, d_Strike_Put+2);
	boost::shared_ptr<Payoff> payoff2(new DigiCall_Put(Strike2, 1, true, Strike_Put, 1, 1));
	sched_atom.setttm(1);
	sched_atom.setPayoff(payoff2);
	ProdSched.push_back(sched_atom);
	
	double d_BasePrice[] = {100, 100};
	std::vector<double> BasePrice(d_BasePrice, d_BasePrice+2);
	
	double d_Corr[] = {1, 0.5, 0.5, 1};
	alglib::real_2d_array Corr;
	Corr.setcontent(2, 2, d_Corr);

	boost::shared_ptr<Product> Prod(new Product(Stock, BasePrice, ProdSched, Corr));	

	double d_Price[] = {S1, S2};
	std::vector<double> Price(d_Price, d_Price+2);
	std::vector<double> MinPrice(d_Price, d_Price+2);
	std::vector<double> MaxPrice(d_Price, d_Price+2);

	BSProcess MC(Price, MinPrice, MaxPrice, Prod);
	double MCresult = 100.0 + MC.Calculate(100000);
	time(&t1);
	
	std::cout << MCresult << " / time: " << difftime(t1, t0) << "\n";

	return 0;
}