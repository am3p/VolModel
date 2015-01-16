#include "stdafx.h"

#include <Codes/Product.h>

Product::Product(std::vector<Underlying> Stock,
				 std::vector<double> BasePrice,
				 std::vector<Schedule> ProdSched,
				 alglib::real_2d_array Corr){
	Stock_ = Stock;
	BasePrice_ = BasePrice;
	ProdSched_ = ProdSched;
	Corr_ = Corr;
};

std::vector<Underlying> Product::getStock(){
	return Stock_;
}

std::vector<Schedule> Product::getSchedule(){
	return ProdSched_;
}

alglib::real_2d_array Product::getCorr(){
	return Corr_;
}