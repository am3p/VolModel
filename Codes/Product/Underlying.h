#include <boost/shared_ptr.hpp>

#include <Process/BSParameter.h>

// Underlying Class
class Underlying
{
public:
	Underlying(const double, 
			   const boost::shared_ptr<InterestRate>,
			   const boost::shared_ptr<DividendRate>, 
			   const boost::shared_ptr<Volatility>);
	double getStock();
	double getInterest(double);
	double getDividend(double);
	double getVolatility(double, double);
	void chkInterestType();
	void chkDividendType();
	void chkVolatilityType();
private:
	double S_;
	boost::shared_ptr<InterestRate> rf_;
	boost::shared_ptr<DividendRate> div_;
	boost::shared_ptr<Volatility> vol_;
	bool isrfterm;
	bool isdivterm;
	bool isvolterm;
	bool isvolsurf;
};