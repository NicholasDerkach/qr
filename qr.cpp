#include <stdio.h>
#include <iostream>
#include <vector>
#include <ctime>
#include <cmath>
#include <iomanip>
#include <omp.h>
#include <cstdlib>

using std::setw;
using std::setprecision;

using namespace std;



std::vector<double> operator-(std::vector<double> vec_a, std::vector<double> vec_b) {
    std::vector<double> result(vec_a.size());
    for (int i = 0; i < vec_a.size(); i++) {
        result[i] = vec_a[i] - vec_b[i];
    }
    return result;
}

double operator*(std::vector<double> vec_a, std::vector<double> vec_b) {
    double result = 0;
    int i;
    for (i = 0; i < vec_a.size(); i++) {
        result += vec_a[i] * vec_b[i];
    }
    return result;
}

std::vector<double> operator*(std::vector<double> vec, double scalar) {
    std::vector<double> result(vec.size());
    int i;
    for (i = 0; i < vec.size(); i++) {
        result[i] = vec[i] * scalar;
    }
    return result;
}


std::vector<double> operator|(std::vector<double> which, std::vector<double> onto) {
    return onto * ((which * onto)/(onto*onto));
}



void printVector(std::vector<double> vec) {
    for (int i = 0; i < vec.size(); i++) {
        std::cout << setprecision(2) << setw(10) << vec[i];
    }
}


void printMatrix(std::vector<std::vector<double>> matrix) {
    for (int i = 0; i < matrix.size(); i++) {
        printVector(matrix[i]);
        std::cout << "\n";
    }
    std::cout << "\n";

}

std::vector<std::vector<double>> generateMartix(int scale) {
    srand(time(NULL));
    std::vector<std::vector<double>> matrix(scale, std::vector<double>(scale));
    for ( int i = 0 ; i < scale ; i++ )
         for ( int j = 0 ; j < scale ; j++ )
            matrix[i][j] = rand() % 99;

    return matrix;
}

std::vector<std::vector<double>> ortogonalization(std::vector<std::vector<double>> input) {
    std::vector<std::vector<double>> result(input.size(), std::vector<double>(input.size()));
    int i, j;
    std::vector<double> intermidiateVector;
    result[0] = input[0];
#pragma omp parallel for shared(result) private(intermidiateVector, j)
    for ( i = 1 ; i < input.size() ; i++ ) {
        intermidiateVector = input[0];
        for ( j = 1 ; j <= i; j++ ) {
            intermidiateVector = intermidiateVector - (input[j] | intermidiateVector);
        }
        result[i] = intermidiateVector;
    }
   
    return result;
}

std::vector<std::vector<double>> leftTransponation(std::vector<std::vector<double>> input) {
    std::vector<std::vector<double>> result(input.size(), std::vector<double>(input.size()));
    for (int i = input.size() - 1; i >= 0 ; i--) {
        for (int j = 0; j < input.size() ; j++) {
            result[input.size() - 1 - i][j] = input[j][i];
        }
    }
   
    return result;
}


std::vector<std::vector<double>> rightTransponation(std::vector<std::vector<double>> input) {
    std::vector<std::vector<double>> result(input.size(), std::vector<double>(input.size()));
    for ( int i = 0 ; i < input.size() ; i++ ) {
        for (int j = input.size() - 1; j >= 0 ; j--) {
            result[i][input.size() - 1 - j] = input[j][i];
        }
    }
   
    return result;
}



int main() {
    int numberOfThreads, poryadok;
    std::cout << "Number of threads = ";
    std::cin >> numberOfThreads;
    std::cout << "Poryadok = ";
    std::cin >> poryadok;
    std::vector<std::vector<double>> matrix = generateMartix(poryadok);
    omp_set_num_threads(numberOfThreads);
    double start = omp_get_wtime();
    std::vector<std::vector<double>> ortogonalizedMatrix = ortogonalization(matrix);
    double end = omp_get_wtime();
    double time = end - start;
    std::cout << "\n";
    std::cout << "OmpTime:" << setprecision(5) << time << endl;
    std::cout << "Done" << endl;
}
