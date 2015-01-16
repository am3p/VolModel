#include <vector>

// Payoff Class
class Payoff{
public:
	virtual void CalcPayoff(std::vector<double>) = 0;
	virtual bool BermudanCond(std::vector<double>) = 0;
	virtual bool TerminalCond(std::vector<double>, std::vector<double>) = 0;
	std::vector<double> getCF();
protected:
	std::vector<double> strike_;
	double sign_;
	double particip_;
	bool CP_; // true: Call, false: Put
	std::vector<double> CF;
};

class Bermudan : public Payoff{
public:
protected:
};

class Vanilla : public Bermudan
{
public:
	Vanilla(std::vector<double>, double, double, bool);
	void CalcPayoff(std::vector<double>);
	bool BermudanCond(std::vector<double>);
	bool TerminalCond(std::vector<double>, std::vector<double>) {return false;};
private:
};

class Digital : public Bermudan
{
public:
	Digital(std::vector<double>, double, bool, double);
	void CalcPayoff(std::vector<double>);
	bool BermudanCond(std::vector<double>);
	bool TerminalCond(std::vector<double>, std::vector<double>) {return false;};
private:
	double Coupon_;
};


// Final Payoff Class
class Terminal : public Payoff{
public:
protected:
	bool DownBTouched;
	bool UpBTouched;
};

class DigiCall_Put : public Terminal{
public:
	DigiCall_Put(std::vector<double>, double, bool, std::vector<double>, double, double);
	bool BermudanCond(std::vector<double>) {return false;};
	bool TerminalCond(std::vector<double>, std::vector<double>);
	void CalcPayoff(std::vector<double>);
private:
	std::vector<double> strike_put_;
	double Coupon_;
	double Coupon_Dummy_;
};
