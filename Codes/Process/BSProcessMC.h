#include <alglibmisc.h>

#include <boost/shared_ptr.hpp>

#include <Product/Product.h>

class BSProcess{
public:	
	BSProcess(std::vector<double>, std::vector<double>, std::vector<double>, boost::shared_ptr<Product>);
	void MC(double, double, alglib::hqrndstate);
	double Calculate(int);
	void Cholesky();
protected:
	void MCinnov(double, double, alglib::hqrndstate);
private:
	std::vector<double> Price_;
	std::vector<double> Price0_;
	std::vector<double>	MinPrice_;
	std::vector<double> MaxPrice_;
	alglib::real_2d_array Corr_Cholesky;
	boost::shared_ptr<Product> Prod_;
};