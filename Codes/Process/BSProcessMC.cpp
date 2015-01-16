#include <stdafx.h>

#include <Process/BSProcessMC.h>

BSProcess::BSProcess(std::vector<double> Price,
					 std::vector<double> MinPrice,
					 std::vector<double> MaxPrice,
					 boost::shared_ptr<Product> Prod){
	 Price_ = Price;
	 Price0_ = Price;
	 MinPrice_ = MinPrice;
	 MaxPrice_ = MaxPrice;
	 Prod_ = Prod;

	 Cholesky();
}

void BSProcess::MCinnov(double t0, double dt, alglib::hqrndstate state){
	alglib::hqrndrandomize(state);
	double rf, q, vol, lnS;

	std::vector<Underlying> Stock = Prod_->getStock();
	for (unsigned int i = 0; i < Stock.size(); i++){
		rf = Stock[i].getInterest(t0);
		q = Stock[i].getDividend(t0);
		vol = Stock[i].getVolatility(t0, Price_[i]);

		//std::cout << i << "/" << Price_[i] << "/" << vol << "\n";
		
		lnS = log(Price_[i]);
		lnS += (rf-q-vol*vol/2.0)*dt;
		for (unsigned int j = 0; j <= i; j++){
			lnS += vol*Corr_Cholesky[i][j]*alglib::hqrndnormal(state)*sqrt(dt);
		}
		Price_[i] = exp(lnS);
		if (MinPrice_[i] > Price_[i]){
			MinPrice_[i] = Price_[i];
		}
		if (MaxPrice_[i] < Price_[i]){
			MaxPrice_[i] = Price_[i];
		}
	}
}
	
void BSProcess::MC(double t0, double t1, alglib::hqrndstate state){
	//alglib::hqrndrandomize(state);
	double dt = 1.0/365.0;
	while(t0 < t1){
		BSProcess::MCinnov(t0, dt, state);
		t0 += dt;
	}
}

void BSProcess::Cholesky(){
	//Cholesky decomposition
	Corr_Cholesky = Prod_->getCorr();
	std::vector<Underlying> Stock = Prod_->getStock();
	alglib::spdmatrixcholesky(Corr_Cholesky, Prod_->getStock().size(), false);
}

double BSProcess::Calculate(int n){
	alglib::hqrndstate state;
	
	double t0, t1;
	double MC_result = 0;
	std::vector<Schedule> Sched = Prod_->getSchedule();
	unsigned int j;
	for (int i = 0; i < n; i++){
		Price_ = Price0_;
		MinPrice_ = Price0_;
		MaxPrice_ = Price0_;
		t0 = 0;
		j = 0;
		do{
			t1 = Sched[j].getttm();
			MC(t0, t1, state);
			if(Sched[j].isTerminate(Price_)){
				MC_result += Sched[j].CalculateCF(Price_);
			}
			t0 = t1;
			j++;
		}while(!(Sched[j-1].isTerminate(Price_)));
	}

	return MC_result / (double)n;
}