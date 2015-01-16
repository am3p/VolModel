#include <ql/Math/array.hpp>

class SABR{
public:
	SABR(QuantLib::Array param) : alpha(param[0]), rho(param[1]), nu(param[2]) {};
	double vol (double F, double K, double t);
private:
	double alpha;
	double rho;
	double nu;
};