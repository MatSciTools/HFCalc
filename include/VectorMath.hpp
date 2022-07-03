#include <vector>
#include <Eigen/Dense>

class VectorMath {
public:
    VectorMath(){};
    static std::vector <long double> subtractVectors(std::vector <long double> a, 
    std::vector <long double> b);
    static std::vector <long double> calMeanVector(std::vector <long double> a, 
    long double alpha, std::vector <long double> b, long double beta);
    static long double calNorm(std::vector <long double> a);
    static long double calDistance(std::vector <long double> a, 
    std::vector <long double> b);
    static long double calExpectationValue(std::vector <long double> C, Eigen::MatrixXd S, int dim);
};