#pragma once
#include <vector>
/*
    Notebook.hpp.
    Author: Omer Rabin.
*/
using namespace std;
#include <ostream>
#include <iostream>

namespace zich {
    class Matrix {
        vector<vector<double>> *matrix; 
    public:
        Matrix(){
        }
        Matrix(vector<double> v, int dim_x, int dim_y) {
            if(dim_y*dim_x!=v.size()){
                 throw runtime_error("bad input for constructor");
             }
            matrix = new vector<vector<double>>((unsigned long)dim_x);
            for (int i = 0 ; i < dim_x ; i++ )
            {
                matrix->at((unsigned long)i).resize((unsigned long)dim_y);
            }
            int index = 0;
            for(int i=0 ; i < dim_x; i++){
                for(int j = 0; j< dim_y; j++){
                    matrix->at((unsigned long)i).at((unsigned long)j) = v[(unsigned long)index];
                    index ++;
                }
            }
        }
        
        ~Matrix(){
        } // destructor
       
    Matrix operator+(Matrix const &other);
    Matrix operator+();
    Matrix operator++(int x);
    Matrix operator++(); 
    Matrix operator-(Matrix const &other);
    Matrix operator-();
    Matrix operator--(int x);
    Matrix operator--(); 
    Matrix operator+=(Matrix const &other);
    Matrix operator-=(Matrix const &other);
    bool operator>(Matrix &other);
    bool operator<(Matrix &other);
    bool operator>=(Matrix &other);
    bool operator<=(Matrix &other);
    bool operator==(const Matrix &other)const;
    bool operator!=(const Matrix &other)const;
    Matrix operator*(Matrix const &other);
    friend Matrix operator*(const double x, const Matrix& mat);
    Matrix operator*(const double scalar);
    Matrix operator*=(Matrix const &other);
    Matrix operator*=(double x);
    friend std::ostream& operator<< (ostream& out, const Matrix& m);
    friend std::istream& operator>> (istream& in, Matrix& m);

    };


}
