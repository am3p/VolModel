#include <vector>

#include <linalg.h>

#include <Schedule/Schedule.h>
#include <Product/Underlying.h>

class Product
{
public:
	Product(std::vector<Underlying>, std::vector<double>, std::vector<Schedule>, alglib::real_2d_array);
	std::vector<Underlying> getStock();
	std::vector<Schedule> getSchedule();
	alglib::real_2d_array getCorr();
private:
	std::vector<Underlying> Stock_;
	std::vector<double> BasePrice_;
	std::vector<Schedule> ProdSched_;
	alglib::real_2d_array Corr_;
};