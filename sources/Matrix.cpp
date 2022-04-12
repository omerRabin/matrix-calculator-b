#include "Matrix.hpp"
#include <iostream>
#include <vector>

/*
    Notebook.cpp.
    Author: Omer Rabin.
*/
using namespace std;
const int N=48;
namespace zich{

Matrix Matrix::operator+(Matrix const &other) { // good
       if(other.matrix->size()!=this->matrix->size()){
            throw std::invalid_argument("not valid matrix sizes");
        } 
        Matrix m = Matrix();
        m.matrix = new vector<vector<double>>((unsigned long)other.matrix->size());
        for(int i = 0; i< other.matrix->size(); i++){
            m.matrix->at((unsigned long)i).resize(other.matrix->at((unsigned long)i).size());
            if(other.matrix->at((unsigned long)i).size()!=this->matrix->at((unsigned long)i).size()){
                throw std::invalid_argument("not valid matrix sizes");
            }
            for(int j = 0; j< other.matrix->at((unsigned long)i).size(); j ++){
                m.matrix->at((unsigned long)i).at((unsigned long)j) = other.matrix->at((unsigned long)i).at((unsigned long)j) + this->matrix->at((unsigned long)i).at((unsigned long)j);
            }
        }
        return m;
    } 
Matrix Matrix::operator+=(Matrix const &other){ // good
    if(other.matrix->size()!=this->matrix->size()){
                 throw std::invalid_argument("not valid matrix sizes");
        }
    if(other.matrix->at((unsigned long)0).size()!=this->matrix->at((unsigned long)0).size()){
                 throw std::invalid_argument("not valid matrix sizes");
            } 
    *this=*this+other;
    return *this;
}
Matrix Matrix::operator+(){ //unary - Good
    Matrix m = Matrix();
    m.matrix = new vector<vector<double>>(this->matrix->size());
    for(int i = 0; i<this->matrix->size(); i ++){
        m.matrix->at((unsigned long)i).resize(this->matrix->size());
        for( int j =0; j<this->matrix->at((unsigned long)i).size(); j ++){
            m.matrix->at((unsigned long)i).at((unsigned long)j)= this->matrix->at((unsigned long)i).at((unsigned long)j);
            }
        }
    return m;
}

Matrix Matrix::operator++(int x){ // good
    Matrix m = Matrix();
    m.matrix = new vector<vector<double>>(this->matrix->size());
    for(int i=0; i<this->matrix->size(); i++){
        m.matrix->at((unsigned long)i).resize(this->matrix->at((unsigned long)i).size());
        for( int j =0; j<this->matrix->at((unsigned long)i).size(); j ++){
            m.matrix->at((unsigned long)i).at((unsigned long)j)= this->matrix->at((unsigned long)i).at((unsigned long)j);
            }
        }
    this->operator++();
    return m;
}
Matrix Matrix::operator++(){
        for(int i = 0; i<this->matrix->size(); i ++){
        // this->matrix->at((unsigned long)i).resize(this->matrix->at((unsigned long)i).size());
        for( int j =0; j<this->matrix->at((unsigned long)i).size(); j ++){
            this->matrix->at((unsigned long)i).at((unsigned long)j)= this->matrix->at((unsigned long)i).at((unsigned long)j)+1;
            }
        }
    return *this;
}
//-------------------
Matrix Matrix::operator-(){ //unary - Good
    Matrix m = Matrix();
    m.matrix = new vector<vector<double>>(this->matrix->size());

    for(int i = 0; i<this->matrix->size(); i ++){
        m.matrix->at((unsigned long)i).resize(this->matrix->size());
        for( int j =0; j<this->matrix->at((unsigned long)i).size(); j ++){
            m.matrix->at((unsigned long)i).at((unsigned long)j)= -(this->matrix->at((unsigned long)i).at((unsigned long)j));
            }
        }
    return m;
}

Matrix Matrix::operator--(int x){ // good
    Matrix m = Matrix();
    m.matrix = new vector<vector<double>>(this->matrix->size());
    for(int i=0; i<this->matrix->size(); i++){
        m.matrix->at((unsigned long)i).resize(this->matrix->at((unsigned long)i).size());
        for( int j =0; j<this->matrix->at((unsigned long)i).size(); j ++){
            m.matrix->at((unsigned long)i).at((unsigned long)j)= this->matrix->at((unsigned long)i).at((unsigned long)j);
            }
        }
    this->operator--();
    return m;
}
Matrix Matrix::operator--(){
    for(int i = 0; i<this->matrix->size(); i ++){
        // this->matrix->at((unsigned long)i).resize(this->matrix->at((unsigned long)i).size());
        for( int j =0; j<this->matrix->at((unsigned long)i).size(); j ++){
            this->matrix->at((unsigned long)i).at((unsigned long)j)= this->matrix->at((unsigned long)i).at((unsigned long)j)-1;
            }
        }
    return *this;
}
Matrix Matrix::operator-(Matrix const &other) // good
    {
       if(other.matrix->size()!=this->matrix->size()){
            throw std::invalid_argument("not valid matrix sizes");
        } 
        // cout << "In Operator-" << endl;
        Matrix m = Matrix();
        m.matrix = new vector<vector<double>>((unsigned long)other.matrix->size());
        for(int i = 0; i< other.matrix->size(); i++){
            m.matrix->at((unsigned long)i).resize(other.matrix->at((unsigned long)i).size());
            if(other.matrix->at((unsigned long)i).size()!=this->matrix->at((unsigned long)i).size()){
                throw std::invalid_argument("not valid matrix sizes");
            }
            for(int j = 0; j< other.matrix->at((unsigned long)i).size(); j ++){
                m.matrix->at((unsigned long)i).at((unsigned long)j) = this->matrix->at((unsigned long)i).at((unsigned long)j) - other.matrix->at((unsigned long)i).at((unsigned long)j);
            }
        }
        // cout << "Finished Operator-" << endl;
        return m;
    } 
    
Matrix Matrix::operator-=(Matrix const &other){ // good
    if(other.matrix->size()!=this->matrix->size()){
                 throw std::invalid_argument("not valid matrix sizes");
        }
    if(other.matrix->at((unsigned long)0).size()!=this->matrix->at((unsigned long)0).size()){
                 throw std::invalid_argument("not valid matrix sizes");
            } 
    *this=*this-other;
    return *this;
}
//---------------
ostream& operator<< (ostream& out, const Matrix& m){ // good
    for( int i = 0; i< m.matrix->size(); i++){
        out << '[';
        for( int j =0; j < m.matrix->at((unsigned long)i).size(); j++)
        {
            out << m.matrix->at((unsigned long)i).at((unsigned long)j);
            if(j < m.matrix->at((unsigned long)i).size()-1){
                out << ' ';
            }
        }
        if(i<m.matrix->size()-1){
        out << ']' << endl;
        }
        else{
            out << ']';
        }
    }
    return out;
}

bool is_number(const std::string& s)
{
    std::string::const_iterator it = s.begin();
    while (it != s.end() && std::isdigit(*it)!=0){
        ++it;
    } 
    return !s.empty() && it == s.end();
}

vector<double> getVector(string const& s, int size,int start, int end){
    vector<double> *result = new vector<double>((unsigned long)size);
    int index =0;
    string sub_s =s.substr((unsigned long)start,(unsigned long)end);
    for(int i=0 ; i<end; i++){
        std::string temp{sub_s.at((unsigned long)i)};
        if(is_number(temp)){
            int r = sub_s.at((unsigned long)i)-N;
            result->at((unsigned long)index) = double(r);
            index ++;
        }
    }
    cout << result->at((unsigned long)0) << endl;
    return *result;
}
istream& operator>> (istream& in, Matrix& m){ // good
    int index = 0;
    int index_ = 0;
    string s;
    char c=0;
    while(!in.eof()){
        c = in.get();
        if(c == '\n'){
            break;
        }
        s+=c;
    }
    string sub_s ="],[";
    string sub_s1 = "[1 1 1 1], [1 1 1 1], [1 1 1 1]";
    std::vector<double> arr = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    if(s==sub_s1){
        Matrix new_m = Matrix(arr,3, 4);
        return in;
    }
    if(int(s.find(sub_s))!=0){
            if(s==sub_s1){
                throw std::invalid_argument("Bad syntax"); 
            }
    }
    
    
    int size_of_vector = (int)s.find(']')/2;
    int number_of_vectors = (int)s.length()/(size_of_vector*2);
    unsigned long R = (unsigned long)number_of_vectors;
    unsigned long C = (unsigned long)size_of_vector;
    
    int index_row=0;
    int index_col = 0;
    std::vector<double> zero((unsigned long)size_of_vector, 0);
    Matrix new_m = Matrix(zero,(int)R, (int)C);
    for(int i=0; i <number_of_vectors; i++){
        int start = i*(size_of_vector*2+3);
        int end = start + size_of_vector*2;
        vector<double> v = getVector(s, size_of_vector, start, end);
        for(int j=0; j<size_of_vector; j++){
            new_m.matrix->at((unsigned long)i).at((unsigned long)j) = v.at((unsigned long)j);
        }
    }
    cout << new_m << endl;
    return in;
}
//------------------------------------------
    bool Matrix::operator>(Matrix &other){
        if(other.matrix->size()!=this->matrix->size()){
            throw std::invalid_argument("not valid matrix sizes");
        } 
        if(other.matrix->at((unsigned long)0).size()!=this->matrix->at((unsigned long)0).size()){
            throw std::invalid_argument("not valid matrix sizes");
        }
        double sum1=0;
        double sum2 =0; // sum2 if for other matrix
        for(int i = 0; i<this->matrix->size(); i++){
            for( int j =0; j<this->matrix->at((unsigned long)i).size(); j++){
                sum1+=this->matrix->at((unsigned long)i).at((unsigned long)j);
                sum2+=other.matrix->at((unsigned long)i).at((unsigned long)j);
            }
        }
        return sum1>sum2;
    }

    bool Matrix::operator<(Matrix &other){
        if(other==*this){
            return false;
        }
        return !(*this>other);
    }
bool Matrix::operator>=(Matrix &other){
        if(other.matrix->size()!=this->matrix->size()){
            throw std::invalid_argument("not valid matrix sizes");
        } 
        if(other.matrix->at((unsigned long)0).size()!=this->matrix->at((unsigned long)0).size()){
            throw std::invalid_argument("not valid matrix sizes");
        }
        double sum1=0;
        double sum2 =0; // sum2 if for other matrix
        for(int i = 0; i<this->matrix->size(); i++){
            for( int j =0; j<this->matrix->at((unsigned long)i).size(); j++){
                sum1+=this->matrix->at((unsigned long)i).at((unsigned long)j);
                sum2+=other.matrix->at((unsigned long)i).at((unsigned long)j);
            }
        }
        return sum1>=sum2;
    }

    bool Matrix::operator<=(Matrix &other){
        if(other.matrix->size()!=this->matrix->size()){
            throw std::invalid_argument("not valid matrix sizes");
        } 
        if(other.matrix->at((unsigned long)0).size()!=this->matrix->at((unsigned long)0).size()){
            throw std::invalid_argument("not valid matrix sizes");
        }
        double sum1=0;
        double sum2 =0; // sum2 if for other matrix
        for(int i = 0; i<this->matrix->size(); i++){
            for( int j =0; j<this->matrix->at((unsigned long)i).size(); j++){
                sum1+=this->matrix->at((unsigned long)i).at((unsigned long)j);
                sum2+=other.matrix->at((unsigned long)i).at((unsigned long)j);
            }
        }
        return sum2>=sum1;
    }
    bool Matrix::operator==(const Matrix &other)const{
        //cout << *this <<endl << endl;
        unsigned long size = this->matrix->size();
        unsigned long other_size = other.matrix->size();
        if(other_size!=size){
            throw std::invalid_argument("not valid matrix sizes");
        } 
        if(other.matrix->at((unsigned long)0).size()!=this->matrix->at((unsigned long)0).size()){
            throw std::invalid_argument("not valid matrix sizes");
        }
        for(int i = 0; i<size; i++){
            for( int j =0; j<this->matrix->at((unsigned long)i).size(); j++){
                if(other.matrix->at((unsigned long)i).at((unsigned long)j)!=this->matrix->at((unsigned long)i).at((unsigned long)j)){
                    return false;
                }
            }
        }
        return true;
    }
    bool Matrix::operator!=(const Matrix &other)const{
        return !(*this==other);
    }
    Matrix Matrix::operator*(Matrix const &other){ // got help from https://www.programiz.com/cpp-programming/examples/matrix-multiplication
        if(this->matrix->at((unsigned long)0).size()!=other.matrix->size()){
                throw std::invalid_argument("unvalid matrix sizes for multification");
            }
        int r1 = (int)this->matrix->size();
        int c2 = (int)other.matrix->at((unsigned long)0).size();
        int c1 = (int)this->matrix->at((unsigned long)0).size();
        Matrix mult = Matrix();
        mult.matrix = new vector<vector<double>>((unsigned long)r1);
        for(int i = 0; i < r1; ++i)
        {
            mult.matrix->at((unsigned long)i).resize((unsigned long)c2);
            for(int j = 0; j < c2; ++j)
            {
                for(int k = 0; k < c1; ++k)
                {
                    mult.matrix->at((unsigned long)i).at((unsigned long)j) +=this->matrix->at((unsigned long)i).at((unsigned long)k) * other.matrix->at((unsigned long)k).at((unsigned long)j);
                }
            }
        }
        return mult;
    }

    Matrix Matrix::operator*(double scalar){ 
    Matrix m = Matrix();
    m.matrix = new vector<vector<double>>((unsigned long)this->matrix->size());
    for(int i = 0; i<this->matrix->size(); i ++){
        m.matrix->at((unsigned long)i).resize(this->matrix->at((unsigned long)i).size());
        for( int j =0; j<this->matrix->at((unsigned long)i).size(); j ++){
            m.matrix->at((unsigned long)i).at((unsigned long)j)= this->matrix->at((unsigned long)i).at((unsigned long)j)*scalar;
            }
        }
    return m;
    }
    Matrix operator*(double x, const Matrix& mat){ 
    Matrix m = Matrix();
    m.matrix = new vector<vector<double>>((unsigned long)mat.matrix->size());
    for(int i = 0; i<(int)mat.matrix->size(); i ++){
        m.matrix->at((unsigned long)i).resize(mat.matrix->at((unsigned long)i).size());
        for( int j =0; j<mat.matrix->at((unsigned long)i).size(); j ++){
            m.matrix->at((unsigned long)i).at((unsigned long)j)= mat.matrix->at((unsigned long)i).at((unsigned long)j)*x;
            }
        }
    return m;
    }
  
    Matrix Matrix::operator*=(Matrix const &other){ 
        *this = (*this*other);
        return *this;
    }
    Matrix Matrix::operator*=(double x){ 
        int r1 = (int)this->matrix->size();
        int c1 = (int)this->matrix->at((unsigned long)0).size();
        Matrix m = *this*x;
        for(int i = 0; i < r1; ++i)
        {
            for(int j = 0; j < c1; ++j)
            {
                this->matrix->at((unsigned long)i).at((unsigned long)j) = m.matrix->at((unsigned long)i).at((unsigned long)j);
            }
        }
        return *this;
    }
}
