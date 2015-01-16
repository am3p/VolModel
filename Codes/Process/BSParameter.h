#include <linalg.h>

#include <boost/shared_ptr.hpp>

// Interest Rate Class
class InterestRate
{
public:
	virtual double getInterest() = 0;
	virtual double getInterest(double) = 0;
private:
};

class InterestConst : public InterestRate
{
public:
	InterestConst(double);
	double getInterest();
	double getInterest(double) {return -1;};
	void setInterest(double);
private:
	double r_;
};

class InterestTerm : public InterestRate
{
public:
	InterestTerm(alglib::real_1d_array, alglib::real_1d_array);
	double getInterest() {return -1;};
	double getInterest(double);
	void setInterest(alglib::real_1d_array, alglib::real_1d_array);
private:
	alglib::real_1d_array t_;
	alglib::real_1d_array r_;
	alglib::spline1dinterpolant s;
};

// Dividend Rate Class
class DividendRate
{
public:
	virtual double getDividend() = 0;
	virtual double getDividend(double) = 0;
private:
};

class DividendConst : public DividendRate
{
public:
	DividendConst(double);
	double getDividend();
	double getDividend(double) {return -1;};
	void setDividend(double);
private:
	double q_;
};

class DividendTerm : public DividendRate
{
public:
	DividendTerm(alglib::real_1d_array, alglib::real_1d_array);
	double getDividend() {return -1;};
	double getDividend(double);
	void setDividend(alglib::real_1d_array, alglib::real_1d_array);
private:
	alglib::real_1d_array t_;
	alglib::real_1d_array q_;
	alglib::spline1dinterpolant s;
};

class DividendDiscrete : public DividendRate
{
public:
	DividendDiscrete(alglib::real_1d_array, alglib::real_1d_array);
private:
	alglib::real_1d_array t_;
	alglib::real_1d_array q_;
};

// Volatility Class
class Volatility
{
public:
	virtual double getVolatility() = 0;
	virtual double getVolatility(double) = 0;
	virtual double getVolatility(double, double) = 0;
private:
};

class VolatilityConst : public Volatility
{
public:
	VolatilityConst(double);
	double getVolatility();
	double getVolatility(double) {return -1;};
	double getVolatility(double, double) {return -1;};
	void setVolatility(double);
private:
	double v_;
};

class VolatilityTerm : public Volatility
{
public:
	VolatilityTerm(alglib::real_1d_array, alglib::real_1d_array);
	double getVolatility() {return -1;};
	double getVolatility(double);
	double getVolatility(double, double) {return -1;};
	void setVolatility(alglib::real_1d_array, alglib::real_1d_array);
private:
	alglib::real_1d_array t_;
	alglib::real_1d_array v_;
	alglib::spline1dinterpolant s;
};

class VolatilitySurf : public Volatility
{
public:
	VolatilitySurf(alglib::real_1d_array, alglib::real_1d_array, alglib::real_1d_array);
	double getVolatility() {return -1;};
	double getVolatility(double) {return -1;};
	double getVolatility(double, double);
	void setVolatility(alglib::real_1d_array, alglib::real_1d_array, alglib::real_1d_array);
private:
	alglib::real_1d_array t_;
	alglib::real_1d_array K_;
	alglib::real_1d_array v_;
	alglib::spline2dinterpolant s;
};