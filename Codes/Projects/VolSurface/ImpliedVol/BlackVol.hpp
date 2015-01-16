
class BlackVol{
public:
	BlackVol(int CP_, double TTM_, double Strike_) : CP(CP_), TTM(TTM_), Strike(Strike_) {};
	double BlackPrice(double F, double rf, double Impv);
	double BlackImpv(double F, double rf, double price);
protected:
	double BlackImpvProb(double F, double rf, double price, double Impv);
private:
	int CP;
	double TTM;
	double Strike;
};