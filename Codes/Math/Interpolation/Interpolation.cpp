
#include <cmath>

#include <Math/Interpolation/Interpolation.h>


double Interpolation1D::Interpolate(double x){
	switch(method_){
		case(Interpolation1D::Linear):
			return Interpolation1D::LinearInterpolate(x);
			break;
		case(Interpolation1D::Cubic):
			return Interpolation1D::CubicInterpolate(x);
			break;
		default:
			return -1;
			break;
	}
}

double Interpolation1D::CubicInterpolate(double x){
	std::vector<double> u(xs_.size()-1), yt(xs_.size());
	std::vector<double>::iterator iter_xs;
	std::vector<double>::iterator iter_ys;
	std::vector<double>::iterator iter_yt;
	std::vector<double>::iterator iter_u;

	double sig, p, u_tmp;

	iter_xs = xs_.begin()+1;
	iter_ys = ys_.begin()+1;
	iter_yt = yt.begin()+1;
	iter_u = u.begin()+1;

	while (iter_xs < xs_.end()-1){
		sig = ((*iter_xs) - (*(iter_xs-1))) / ((*(iter_xs+1)) - (*(iter_xs-1)));
		p = sig * (*(iter_yt-1)) + 2;
		(*iter_yt) = ((sig-1.0)/p);
		u_tmp = ((*(iter_ys+1)) - (*iter_ys)) / ((*(iter_xs+1))-(*iter_xs)) - 
				((*iter_ys) - (*(iter_ys-1))) / ((*iter_xs)-(*(iter_xs-1))); 
		(*iter_u) = ((6*u_tmp) / (*(iter_xs+1)-*(iter_xs-1)) - sig*(*(iter_u-1)))/p;
		iter_xs++; iter_ys++; 
		iter_yt++; iter_u++;
	}
	yt[yt.size()-1] = 0;

	iter_u = u.end()-1;
	iter_yt = yt.end()-2;
	
	while (iter_yt != yt.begin()){
		(*iter_yt) = (*iter_yt) * (*(iter_yt+1)) + (*iter_u);
		iter_yt--;
		iter_u--;
	}

	double h, b, a, result;
	iter_xs = xs_.begin();
	iter_ys = ys_.begin();
	iter_yt = yt.begin();
	while(*iter_xs < x && iter_xs < xs_.end()-1){
		iter_xs++;
		iter_ys++;
		iter_yt++;
	}
	h = (*iter_xs) - (*(iter_xs-1));
	a = ((*iter_xs) - x)/h;
	b = (x - (*(iter_xs-1)))/h;
	
	result = a*(*(iter_ys-1)) + b*(*iter_ys) + ((std::pow(a,3)-a)*(*(iter_yt-1)) + (std::pow(b,3)-b)*(*iter_yt))*pow(h,2)/6;
	return result;
}

double Interpolation1D::LinearInterpolate(double x){
	double h, a, b, result;
	std::vector<double>::iterator iter_xs;
	std::vector<double>::iterator iter_ys;
	iter_xs = xs_.begin();
	iter_ys = ys_.begin();

	while(*iter_xs < x && iter_xs < xs_.end()){
		iter_xs++;
		iter_ys++;
	}

	h = (*iter_xs) - (*(iter_xs-1));
	a = ((*iter_xs) - x)/h;
	b = (x - (*(iter_xs-1)))/h;
	result = a*(*(iter_ys-1)) + b*(*iter_ys);

	return result;
}

double Interpolation2D::Interpolate(double x, double y){
	if(method_ == Interpolation2D::Linear){
		
	}
	else if(method_ == Interpolation2D::Cubic){

	}
	return 0;
}

double Interpolation2D::LinearInterpolate(double x, double y){
	double hx, hy, a, b, p, q, result;
	std::vector<double>::iterator iter_xs;
	std::vector<double>::iterator iter_ys;
	std::vector<std::vector<double>>::iterator iter_zs_row;
	std::vector<double>::iterator iter_zs_col;
	iter_xs = xs_.begin();
	iter_ys = ys_.begin();

	iter_zs_row = zs_.begin();
	while(*iter_xs < x && iter_xs < xs_.end()){
		iter_xs++;
		iter_zs_row++;
	}

	iter_zs_col = iter_zs_row->begin();
	while(*iter_ys < y && iter_ys < ys_.end()){
		iter_ys++;
		iter_zs_col++;
	}

	hx = (*iter_xs) - (*(iter_xs-1));
	hy = (*iter_ys) - (*(iter_ys-1));
	a = ((*iter_xs) - x)/hx;
	b = (x - (*(iter_xs-1)))/hx;
	p = ((*iter_ys) - y)/hy;
	q = (y - (*(iter_ys-1)))/hy;

	double A, B, C, D;
	//A = (iter_zs_row->iter_zs_col);
	//B = (iter_zs_row->(iter_zs_col-1));
	//C = ((iter_zs_row-1)->iter_zs_col);
	//D = ((iter_zs_row-1)->(iter_zs_col-1));

	return result;
}