
class BlackScholesVol{
public:
	BlackScholesVol(int CP_, double TTM_, double Strike_) : CP(CP_), TTM(TTM_), Strike(Strike_) {};
	double BlackScholesPrice(double S, double rf, double q, double Impv);
	double BlackScholesImpv(double S, double rf, double q, double price);
protected:
	double BlackScholesImpvProb(double S, double rf, double q, double price, double Impv);
private:
	int CP;
	double TTM;
	double Strike;
};