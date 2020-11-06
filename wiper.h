#include <vector>

#define DIM 3

class Wiper{
public:
	Wiper();
	void setTransform();
private:
	double gamma_0;
	double beta_0;
	double id_norm[DIM];
	double norm[DIM];
	int rows;
	int cols;
	int tot;
	double max_h;
	double width;
	double height;
	double dx;
	double dy;
	double y;
	std::vector<std::vector<double>> id_top_points;
	std::vector<std::vector<double>> top_points;
};
