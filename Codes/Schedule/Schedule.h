#include <vector>

#include <boost/shared_ptr.hpp>

#include <Product/Payoff.h>

class Schedule
{
public:
	Schedule() {};
	Schedule(boost::shared_ptr<Payoff>, double);

	double getttm();
	boost::shared_ptr<Payoff> getPayoff();

	void setttm(double);
	void setPayoff(boost::shared_ptr<Payoff>);

	void isBermudan();
	bool isTerminate(std::vector<double>);
	double CalculateCF(std::vector<double>);
private:
	boost::shared_ptr<Payoff> payoff_;
	bool isBermudan_;
	double ttm_;
};