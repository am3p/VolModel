#include <vector>

class Interpolation1D{
public:
	Interpolation1D(std::vector<double> xs, std::vector<double> ys) : xs_(xs), ys_(ys) {};
	double Interpolate(double x);
	void data_x(std::vector<double> xs) {xs_ = xs;}
	void data_y(std::vector<double> ys) {ys_ = ys;}
	enum Method{
		Linear,
		Cubic,
	};
	void setmethod(Method method) {method_ = method;}
private:
	double LinearInterpolate(double x);
	double CubicInterpolate(double x);

	Method method_;
	std::vector<double> xs_;
	std::vector<double> ys_;
};

class Interpolation2D{
public:
	double Interpolate(double x, double y);
	void Input(std::vector<double> xs, std::vector<double> ys, std::vector<std::vector<double>> zs)
		{xs_ = xs; ys_ = ys; zs_ = zs;};
	enum Method{
		Linear,
		Cubic,
	};
private:
	double LinearInterpolate(double x, double y);

	Method method_;
	std::vector<double> xs_;
	std::vector<double> ys_;
	std::vector<std::vector<double>> zs_;
};