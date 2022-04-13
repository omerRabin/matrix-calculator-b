#include "Matrix.hpp"
#include <iostream>
#include <vector>
#include <cstring>  

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
int getSecondOccurance(string s,char ch){
    int counter= 0;
    for(int i=0; i<s.length();i++){
        if(s.at((unsigned long)i)==ch){
            counter++;
            if(counter==2){
                return i;
            }
        }
    }
    return -1;
}
bool checkEscapeForVector(string s, int vectorSize, int indexStart){
    int counter = 0;
    for(int i=indexStart ; i<s.length();i+=2){
        if(s.at((unsigned long)i)!=' '){
            return false;
        } 
        counter++;
        if ( counter == vectorSize -1 ){
        return true;
        }
    }
    return false;
}
bool checkEscape(string const &s, int vectorSize){
    for ( int i =2 ; i < s.length(); i+=2*vectorSize-1  + 4 ){
        if(!checkEscapeForVector(s, vectorSize, i)){
            return false;
        }
    }
    return true;
}
bool checkSogar(string s, char sogar, int range, int firstIndex){
    for(int i = firstIndex ; i< s.length(); i+=range){
        if(s.at((unsigned long)i)!= sogar){
            return false;
        }
    }
    return true;
}
bool checkComma(string s, int range, int firstIndex){
    for(int i = firstIndex ; i< s.length(); i+=range){
        if(s.at((unsigned long)i)!= ','){
            return false;
        }
    }
    return true;
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
//"[1 1 1 1], [1 1 1 1], [1 1 1 1]\n"
    int size_of_vector = (int)s.find(']')/2;
    if (getSecondOccurance(s,'[') ==-1 ||getSecondOccurance(s,']')==-1){
        throw std::invalid_argument("Bad syntax of input1"); 
    } 
    int rangeOpenSogar = getSecondOccurance(s,'[')- (int)s.find('[');
    int rangeCloseSogar = getSecondOccurance(s, ']') - (int)s.find(']');
    int rangeComma = getSecondOccurance(s, ',') - (int)s.find(',');
    if(!checkEscape(s, size_of_vector)){
        throw std::invalid_argument("Bad syntax of input2"); 
    }
    if(!checkSogar(s,'[',rangeOpenSogar, 0)){
        throw std::invalid_argument("Bad syntax of input3"); 
    }
    if(!checkSogar(s,']',rangeCloseSogar, s.find(']'))){
        throw std::invalid_argument("Bad syntax of input4"); 
    }
    if(!checkComma(s,rangeComma,s.find(','))){
        throw std::invalid_argument("Bad syntax of input5"); 
    }
    int number_of_vectors = (int)(s.length())/(size_of_vector*2);
    if(number_of_vectors%2 ==0){
        number_of_vectors -=1;
    }
    unsigned long R = (unsigned long)number_of_vectors;
    unsigned long C = (unsigned long)size_of_vector;
    int index_row=0;
    int index_col = 0;
    std::vector<double> zero(C*R, 0);
    int index_vector = 0;
    for(int i=0; i <s.length(); i++){
        char ch_temp = s.at((unsigned long)i);
        std::string str_temp{ch_temp};
        if(is_number(str_temp)){
            double double_temp = (double)ch_temp;
            zero.at((unsigned long)index_vector) = double_temp-N;
            index_vector++;
        }
    }
    Matrix new_m = Matrix(zero,(int)R, (int)C);
    m = new_m;
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
