#include <iostream>
#include <math.h>
#include <vector>
#include "VectorMath.hpp"

std::vector <long double> VectorMath::subtractVectors(std::vector <long double> a, 
std::vector <long double> b){
    std::vector <long double> diff; long double temp; int i;
    if (a.size() != b.size()){
        std::cerr << "The two vectors are of different dimensions!" << a.size() << "  " 
        << b.size() << std::endl;
    }
    for (i=0;i<a.size();i++){
        diff.push_back(a[i] - b[i]);
    }
    return diff;
}

std::vector <long double> VectorMath::calMeanVector(std::vector <long double> a, 
long double alpha, std::vector <long double> b, long double beta) {
    std::vector <long double> P; int i;
    for (i=0;i<a.size();i++){
        P.push_back((alpha*a[i] + beta*b[i])/(alpha + beta));
    }
    return P;
}
long double VectorMath::calNorm(std::vector <long double> a){
    long double norm = 0.0; int i;
    for(i=0;i<a.size();i++){
        norm = norm + pow(a[i], 2);
    }
    return sqrt(norm);
}

long double VectorMath::calDistance(std::vector <long double> a, 
std::vector <long double> b){
    std::vector <long double> diff = VectorMath::subtractVectors(a, b);
    long double norm = VectorMath::calNorm(diff);
    return norm;
}

long double VectorMath::calExpectationValue(std::vector <long double> C, Eigen::MatrixXd S, int dim){
    long double norm = 0.0;
    int i, j;
    for (i=0;i<dim;i++){
        for (j=0;j<dim;j++){
            norm = norm + C[i]*S(i, j)*C[j];
        }
    }
    return norm;
}