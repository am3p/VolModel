struct MarketParam{
    float ttm;
    float F;
    float rf;
};

struct OptData{
	int CP;
	float TTM;
	float Strike;
	float price;
	float impv;
};

struct SABRparam{
	float TTM;
	float alpha;
	float rho;
	float nu;
	float F;
	float rf;
};

struct SABRparamSeries{
	struct SABRparam p[20];
	float t[20];
	int Size;
};

struct SABRData{
	int size;
	struct SABRparam p;
	struct OptData opt[500];
};

struct RateCurve{
	float t[50];
	float r[50];
	int size;
};

struct VolSurf{
	float t[20];
	float K[20];
	float vol[400];
};