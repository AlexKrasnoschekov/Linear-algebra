#ifndef LINAL_HEADER_H_INCLUDED
#define LINAL_HEADER_H_INCLUDED

#include <iostream>
#include <cmath>
template <typename T>
class Point{
    T x_;
    T y_;
public:
    Point() {}
    ~Point() {
        x_ = 0;
        y_ = 0;
    }
    Point(const T& x, const T& y) {
        x_ = x;
        y_ = y;
    }
    Point(const Point& other) {
        x_ = other.x_;
        y_ = other.y_;
    }
    const Point& operator=(const Point<T>& other) {
        x_ = other.x_;
        y_ = other.y_;
        return *this;
    }
    const Point& operator=(Point<T>& other) {
        x_ = other.x_;
        y_ = other.y_;
        return *this;
    }
    T GetX() {
        return x_;
    }
    T GetY() {
        return y_;
    }
};
template <typename T>
class Line {
    T A_;
    T B_;
    T C_;
public:
    Line() {}
    Line<T>(Point<T>& a, Point<T>& b) {
        if (a.GetX() != b.GetX() && a.GetY() != b.GetY()) {
            A_ = 1;
            B_ = (a.GetX() - b.GetX()) / (b.GetY() - a.GetY());
            C_ = - a.GetX() - (a.GetX() - b.GetX()) / (b.GetY() - a.GetY()) * a.GetY();
        } else if (a.GetX() == b.GetX() && a.GetY() != b.GetY()) {
            A_ = 1;
            B_ = 0;
            C_ = - a.GetX();
        } else if (a.GetX() != b.GetX() && a.GetY() == b.GetY()) {
            A_ = 0;
            B_ = 1;
            C_ = - a.GetY();
        }
    }
    Line<T>(const T& a, const T& b, const T& c) {
        A_ = a;
        B_ = b;
        C_ = c;
    }
    Line(Line& other) {
        A_ = other.A_;
        B_ = other.B_;
        C_ = other.C_;
    }
    const T GetA() {
        return A_;
    }
    const T GetB() {
        return B_;
    }
    const T GetC() {
        return C_;
    }
    T GetIntersectionX(Line<T>& other) {
        if (GetA() != 0 && GetB() != 0 && other.GetA() != 0 && other.GetB() != 0) {
            return  - (GetC() + GetB() * GetIntersectionY(other)) / GetA();
        } else if (GetA() == 0 && GetB() != 0 && other.GetA() != 0 && other.GetB() != 0) {
            return  - other.GetC() / other.GetA() - other.GetB() * GetIntersectionY(other) / other.GetA();
        } else if (GetA() != 0 && GetB() == 0 && other.GetA() != 0 && other.GetB() != 0) {
            return - GetC() / GetA();
        } else if (GetA() != 0 && GetB() != 0 && other.GetA() == 0 && other.GetB() != 0) {
            return - GetC() / GetA() - GetB() * GetIntersectionY(other) / GetA();
        } else if (GetA() != 0 && GetB() != 0 && other.GetA() != 0 && other.GetB() == 0) {
            return - other.GetC() / other.GetA();
        }
    }
    T GetIntersectionY(Line<T>& other) {
        if (GetA() != 0 && GetB() != 0 && other.GetA() != 0 && other.GetB() != 0) {
            return (GetC() * other.GetA() / GetA() - other.GetC()) / (other.GetB() - other.GetA() * GetB()/ GetA());
        } else if (GetA() == 0 && GetB() != 0 && other.GetA() != 0 && other.GetB() != 0) {
            return - GetC() / GetB();
        } else if (GetA() != 0 && GetB() == 0 && other.GetA() != 0 && other.GetB() != 0) {
            return - other.GetC() / other.GetB() - other.GetA() * GetIntersectionX(other) / other.GetB();
        } else if (GetA() != 0 && GetB() != 0 && other.GetA() == 0 && other.GetB() != 0) {
            return - other.GetC() / other.GetB();
        } else if (GetA() != 0 && GetB() != 0 && other.GetA() != 0 && other.GetB() == 0) {
            return - GetC() / GetB() - GetA() * GetIntersectionX(other) / GetB();
        }
    }
    T GetDirectingX() {
        return B_;
    }
    T GetDirectingY() {
        return - A_;
    }
    T Interval(const Line<T>& other) {
        if (GetB() != 0 && GetA() != 0) {
            return fabs((GetC() / GetB()) - (other.GetC()) / other.GetB()) / sqrt(1 + (GetA() * GetA()) / (GetB() * GetB()));
        } else if (GetA() != 0) {
            return fabs(GetC() / GetA() - other.GetC() / other.GetA());
        } else if (GetB() != 0) {
            return fabs(GetC() / GetB() - other.GetC() / other.GetB());
        }
    }
    ~Line<T>() {
        A_ = 0;
        B_ = 0;
        C_ = 0;
    }
    bool CheckInclusion(Point<T>& point) {
        if (GetA() * point.GetX() + GetB() * point.GetY() + GetC() == 0) {
            return true;
        } else {
            return false;
        }
    }
    bool CheckParallels(Line<T>& other) {
        if ((GetA() / GetB()) == (other.GetA() / other.GetB()) && GetC() * GetB() / GetA() != other.GetC() * other.GetB() / other.GetA()) {
            return true;
        } else {
            return false;
        }
    }
    bool operator==(Line<T>& other) {
         if (GetA() == other.GetA() && GetB() == other.GetB() && GetC() == other.GetC()) {
            return true;
         } else {
            return false;
         }
    }
    bool operator!=(Line<T>& other) {
         if (GetA() != other.GetA() || GetB() != other.GetB() || GetC() != other.GetC()) {
            return true;
         } else {
            return false;
         }
    }
};
template <typename T>
class Vector {
    Point<T> start_;
    Point<T> end_;
public:
    T StartX() {
        return start_.GetX();
    }
    T EndX() {
        return end_.GetX();
    }
    T StartY() {
        return start_.GetY();
    }
    T EndY() {
        return end_.GetY();
    }
    Vector();
    Vector(const Point<T>& a, const Point<T>& b) {
        start_ = a;
        end_ = b;
    }
    Vector(Point<T>& a, Point<T>& b) {
        start_ = a;
        end_ = b;
    }
    Vector(Vector& other) {
        start_ = other.start_;
        end_ = other.end_;
    }
    Vector(const Vector& other) {
        start_ = other.start_;
        end_ = other.end_;
    }
    ~Vector() {
        start_.~Point();
        end_.~Point();
    }
    T GetXCoordinate() {
        return end_.GetX() - start_.GetX();
    }
    T GetYCoordinate() {
        return end_.GetY() - start_.GetY();
    }
    T GetLength();
    const Vector<T>& operator=(const Vector<T>& other) {
        start_ = other.start_;
        end_ = other.end_;
        return *this;
    }
    const Vector<T>& operator=(Vector<T>& other) {
        start_ = other.start_;
        end_ = other.end_;
        return *this;
    }
    const Vector<T>& operator+=(Vector<T>& other) {
        T coord_x = end_.GetX() + other.GetXCoordinate();
        T coord_y = end_.GetY() + other.GetYCoordinate();
        end_ = Point<T>(coord_x, coord_y);
        return * this;
    }
    const Vector<T>& operator+(Vector<T>& second) {
        *this += second;
        return *this;
    }
    bool CheckCoordinates(Point<T>& point) {
        if (StartX() <= EndX() && StartY() <= EndY()) {
            if (StartX() <= point.GetX() && EndX() >= point.GetX() && StartY() <= point.GetY() && EndY() >= point.GetY()) {
                return true;
            } else {
                return false;
            }
        } else if (StartX() <= EndX() && StartY() >= EndY()) {
            if (StartX() <= point.GetX() && EndX() >= point.GetX() && StartY() >= point.GetY() && EndY() <= point.GetY()) {
                return true;
            } else {
                return false;
            }
        } else if (StartX() >= EndX() && StartY() <= EndY()) {
            if (StartX() >= point.GetX() && EndX() <= point.GetX() && StartY() <= point.GetY() && EndY() >= point.GetY()) {
                return true;
            } else {
                return false;
            }
        } else if (StartX() >= EndX() && StartY() >= EndY()) {
            if (StartX() >= point.GetX() && EndX() <= point.GetX() && StartY() >= point.GetY() && EndY() <= point.GetY()) {
                return true;
            } else {
                return false;
            }
        } else {
            return false;
        }
    }
    bool CheckInclusion(Point<T>& point) {
        Line<T> line(start_, end_);
        if (line.CheckInclusion(point) && CheckCoordinates(point)) {
            return true;
        } else {
            return false;
        }
    }
    bool CheckIntersection(Vector<T>& other) {
        Line<T> line1(start_, end_);
        Line<T> line2(other.start_, other.end_);
        if (line1.CheckParallels(line2) || line2.CheckParallels(line1)) {
            return false;
        } else {
            if (VectorProduct(*this, other) == 0) {
                if (CheckCoordinates(other.start_) || CheckCoordinates(other.end_) || other.CheckCoordinates(start_) || other.CheckCoordinates(end_)) {
                    return true;
                } else {
                    return false;
                }
            } else {
                Vector<T> v_1(start_, other.start_);
                Vector<T> v_2(start_, other.end_);
                Vector<T> v_3(other.start_, start_);
                Vector<T> v_4(other.start_, end_);
                if ((VectorProduct(*this, v_1) >= 0 && VectorProduct(*this, v_2) <= 0 || VectorProduct(*this, v_1) <= 0 && VectorProduct(*this, v_2) >= 0) &&
                    (VectorProduct(other, v_3) >= 0 && VectorProduct(other, v_4) <= 0 || VectorProduct(other, v_3) <= 0 && VectorProduct(other, v_4) >= 0)) {
                        return true;
                    } else {
                        return false;
                    }
            }
        }
    }
};

template <typename T>
T Vector<T>::GetLength() {
    return sqrt(GetXCoordinate() * GetXCoordinate() + GetYCoordinate() * GetYCoordinate());
}
template <typename T>
T ScalarProduct(Vector<T>& first_vector, Vector<T>& second_vector) {
    return first_vector.GetXCoordinate() * second_vector.GetXCoordinate() + first_vector.GetYCoordinate() * second_vector.GetYCoordinate();
}
template <typename T>
T VectorProduct(Vector<T>& first_vector, Vector<T>& second_vector) {
    return (first_vector.GetXCoordinate() * second_vector.GetYCoordinate() - first_vector.GetYCoordinate() * second_vector.GetXCoordinate());
}
class MatrixAllocationError : std::exception {
};
class MatrixWrongSizeError : std::exception {
};
class MatrixIndexError : std::exception{
};

class MatrixIsDegenerateError : std::exception {
};
//=============== Matrix class ===============//
template <typename T>
class Matrix {
protected:
    int rowsCnt_;
    int colsCnt_;
    T** array_;
public:
    friend Matrix<T> operator+(const Matrix<T>& left, const Matrix<T>& right){
        if (left.rowsCnt_ != right.rowsCnt_ || left.colsCnt_ != right.colsCnt_){
            throw MatrixWrongSizeError();
        }
        Matrix<T> result(left.rowsCnt_, left.colsCnt_);
        for(int i = 0; i < left.rowsCnt_; ++i){
            for(int j = 0; j < left.colsCnt_; ++j){
                result.array_[i][j] = left.array_[i][j] + right.array_[i][j];
            }
        }
        return result;
    }
    friend Matrix<T> operator-(const Matrix<T>& left, const Matrix<T>& right) {
        if (left.rowsCnt_ != right.rowsCnt_ || left.colsCnt_ != right.colsCnt_) {
            throw MatrixWrongSizeError();
        }
        Matrix<T> result(left.rowsCnt_, left.colsCnt_);
        for(int i = 0; i < left.rowsCnt_; ++i) {
            for(int j = 0; j < left.colsCnt_; ++j){
                result.array_[i][j] = left.array_[i][j] - right.array_[i][j];
            }
        }
        return result;
    }
    friend Matrix<T> operator*(const Matrix<T>& left, const T digit) {
        Matrix<T> result(left.rowsCnt_, left.colsCnt_);
        for(int i = 0; i < left.rowsCnt_; ++i) {
            for(int j = 0; j < left.colsCnt_; ++j) {
                result.array_[i][j] = left.array_[i][j] * digit;
            }
        }
        return result;
    }
    friend Matrix<T> operator*(const Matrix<T>& left, const Matrix<T>& right) {
        if (left.colsCnt_ != right.rowsCnt_) {
            throw MatrixWrongSizeError();
        }
        Matrix<T> result(left.rowsCnt_, right.colsCnt_);
        for(int i = 0; i < result.rowsCnt_; ++i){
            for(int j = 0; j < result.colsCnt_; ++j) {
                for(int k = 0; k < left.colsCnt_; ++k) {
                    result.array_[i][j] += left.array_[i][k] * right.array_[k][j];
                }
            }
        }
        return result;
    }
    friend Matrix<T> operator*(const T digit, const Matrix<T>& right){
        Matrix<T> result(right.rowsCnt_, right.colsCnt_);
        for(int i = 0; i < right.rowsCnt_; ++i){
            for(int j = 0; j < right.colsCnt_; ++j){
                result.array_[i][j] = right.array_[i][j] * digit;
            }
        }
        return result;
    }
    friend std::istream&operator>>(std::istream& in, Matrix<T>& matrix) {
        for(int i = 0; i < matrix.rowsCnt_; ++i) {
            for(int j = 0; j < matrix.colsCnt_; ++j) {
                std::cin >> matrix.array_[i][j];
            }
        }
        return in;
    }
    friend std::ostream&operator<<(std::ostream& out, const Matrix<T>& matrix) {
        for(int i = 0; i < matrix.rowsCnt_; ++i) {
            for(int j = 0; j < matrix.colsCnt_; ++j) {
                std::cout << matrix.array_[i][j] << ' ';
                if (j == matrix.colsCnt_ - 1){
                    std::cout << std::endl;
                }
            }
        }
        return out;
    }
    Matrix<T>(const Matrix<T>& matrix) {
        rowsCnt_ = matrix.rowsCnt_;
        colsCnt_ = matrix.colsCnt_;
        array_ = new T*[rowsCnt_];
        for(int i = 0; i < rowsCnt_; ++i) {
            array_[i] = new T[colsCnt_];
        }
        for(int i = 0; i < rowsCnt_; ++i) {
            for(int j = 0; j < colsCnt_; ++j) {
                array_[i][j] = matrix.array_[i][j];
            }
        }
    }
    Matrix<T>(const int rows, const int columns) {
        rowsCnt_ = rows;
        colsCnt_ = columns;
        array_ = new T*[rowsCnt_];
        for(int i = 0; i < rowsCnt_; ++i) {
            array_[i] = new T[colsCnt_];
        }
        for(int i = 0; i < rowsCnt_; ++i) {
            for(int j = 0; j < colsCnt_; ++j){
                array_[i][j] = 0;
            }
        }
    }
    ~Matrix<T>() {
        for(int i = 0; i < rowsCnt_; ++i) {
            delete[] array_[i];
        }
        delete[] array_;
    }
    int getRowsNumber() const {
        return  rowsCnt_;
    }
    int getColumnsNumber() const {
        return colsCnt_;
    }
    Matrix<T>&operator=(const Matrix<T>& right) {
        if (&right == this) {
            return *this;
        }
        if(right.rowsCnt_ != getRowsNumber() || right.colsCnt_ != getColumnsNumber()) {
            for(int i = 0; i < getRowsNumber(); ++i) {
                delete[] array_[i];
            }
            delete[] array_;
            rowsCnt_ = right.getRowsNumber();
            colsCnt_ = right.getColumnsNumber();
            array_ = new T*[getRowsNumber()];
            for(int i = 0; i < getRowsNumber(); ++i) {
                array_[i] = new T[getColumnsNumber()];
            }
        }
        for(int i = 0; i < getRowsNumber(); ++i) {
            for(int j = 0; j < getColumnsNumber(); ++j) {
                array_[i][j] = right.array_[i][j];
            }
        }
        return *this;
    }
    Matrix<T> getTransposed() const {
        Matrix<T> result(getColumnsNumber(), getRowsNumber());
        for (int i = 0; i < getRowsNumber(); ++i) {
            for (int j = 0; j < getColumnsNumber(); ++j) {
                result.array_[j][i] = array_[i][j];
            }
        }
        return result;
    }
    Matrix<T>& transpose() {
        *this = this->getTransposed();
        return *this;
    }
    Matrix<T>&operator+=(const Matrix<T>& right) {
        Matrix<T> temp = *this;
        *this = temp + right;
        return *this;
    }
    Matrix<T>&operator-=(const Matrix<T>& right) {
        Matrix<T> temp = *this;
        *this = temp - right;
        return *this;
    }
    Matrix<T>&operator*=(const T& n) {
        for(int i = 0; i < getRowsNumber(); ++i) {
            for(int j = 0; j <getColumnsNumber(); ++j) {
                array_[i][j] *= n;
            }
        }
        return *this;
    }
    Matrix<T>&operator*=(const Matrix<T>& right) {
        Matrix<T> temp = *this;
        *this = temp * right;
        return *this;
    }
    T& operator()(const int& i, const int& j) {
        if(i >= getRowsNumber() || j >= getColumnsNumber()) {
            throw MatrixIndexError();
        }
        if (i < 0 || j < 0) {
            throw MatrixIndexError();
        }
        return array_[i][j];
    }
    T operator()(const int& i, const int& j) const {
        if(i >= getRowsNumber() || j >= getColumnsNumber()) {
            throw MatrixIndexError();
        }
        if (i < 0 || j < 0) {
            throw MatrixIndexError();
        }
        return array_[i][j];
    }
};
//=============== SquareMatrix class ===============//
template <typename T>
class SquareMatrix : public Matrix<T> {
public:
    int getSize() const {
        return this->rowsCnt_;
    }
    SquareMatrix<T>(const SquareMatrix<T>& matrix) {
        this->rowsCnt_ = matrix.rowsCnt_;
        this->colsCnt_ = this->rowsCnt_;
        for (int i = 0; i < this->rowsCnt_; ++i) {
            for (int j = 0; j < this->colsCnt_; ++j) {
                this->array_[i][j] = matrix.array_[i][j];
            }
        }
    }
    explicit SquareMatrix<T>(const int size) : Matrix<T>(size, size) {}
    ~SquareMatrix<T>() {
        for(int i = 0; i < this->rowsCnt_; ++i) {
            delete[] this->array_[i];
        }
        delete[] this->array_;
    }
    SquareMatrix<T> operator=(const SquareMatrix<T>& right) {
        static_cast<Matrix<T>>(*this) = static_cast<Matrix<T>>(right);
    }
    friend SquareMatrix<T> operator+(const SquareMatrix<T>& left, const SquareMatrix<T>& right) {
        return static_cast<SquareMatrix<T>>(static_cast<Matrix<T>>(left) + static_cast<Matrix<T>>(right));
    }
    friend SquareMatrix<T> operator-(const SquareMatrix<T>& left, const SquareMatrix<T>& right) {
        return static_cast<SquareMatrix<T>>(static_cast<Matrix<T>>(left) - static_cast<Matrix<T>>(right));
    }
    friend SquareMatrix<T> operator*(const SquareMatrix<T>& right, const T digit) {
        return static_cast<SquareMatrix<T>>(static_cast<Matrix<T>>(right) * digit);
    }
    friend SquareMatrix<T> operator*(const T digit, const SquareMatrix<T>& right) {
        return static_cast<SquareMatrix<T>>(digit * static_cast<Matrix<T>>(right));
    }
    friend SquareMatrix<T> operator*(const SquareMatrix<T>& left, const SquareMatrix<T>& right) {
        return static_cast<SquareMatrix<T>>(static_cast<Matrix<T>>(left) * static_cast<Matrix<T>>(right));
    }
    friend std::ostream&operator<<(std::ostream& out, SquareMatrix<T>& matrix) {
        std::cout << static_cast<Matrix<T>>(matrix);
    }
    friend std::istream&operator>>(std::istream& in, SquareMatrix<T>& matrix) {
        std::cin >> static_cast<Matrix<T>>(matrix);
    }
    T operator()(const int& i, const int& j) const {
        return static_cast<Matrix<T>>(*this)(i, j);
    }
    SquareMatrix<T>&operator+=(const SquareMatrix<T>& other) {
        return (*this) + other;
    }
    SquareMatrix<T>&operator-=(const SquareMatrix<T>& other) {
        return (*this) - other;
    }
    SquareMatrix<T>&operator*=(const SquareMatrix<T>& other) {
        return (*this) * other;
    }
    SquareMatrix<T>&operator*=(const T digit) {
        return  digit * (*this);
    }
    SquareMatrix<T> getTransposed() const {
        return static_cast<SquareMatrix<T>>(static_cast<Matrix<T>>(*this).getTransposed());
    }
    const T getDeterminant() const {
        SquareMatrix<T> matrix = *this;
        T determinant = 0;
        return determinant;
    }
    const SquareMatrix<T> GetInverse() const {
        if (this->Determinant() == 0) {
            throw MatrixIsDegenerateError();
        }
        int size_ = getSize();
        SquareMatrix<T> E(size_);
        SquareMatrix<T> answer = *this;
        for (int i = 0; i < size_; ++i) {
            for (int j = 0; j < size_; ++j) {
                if (i == j) {
                    E(i, j) += 1;
                }
            }
        }
        for(int i = 0; i < size_; ++i) {
            T tmp = answer.array_[i][i];
            for(int j = size_; j > 0; --j) {
                E(i, j - 1) /= tmp;
                answer(i, j - 1) /= tmp;
            }
            for(int j = 0; j < size_ ; ++j) {
                if (j != i) {
                    tmp = answer(j, i);
                    for (int k = size_; k > 0; --k) {
                        E(j, k - 1) -= E(i, k - 1) * tmp;
                        answer(j, k - 1) -= answer(i, k - 1) * tmp;
                    }
                }
            }
        }
        for (int i = 0; i < size_; ++i) {
            for (int j = 0; j < size_; ++j) {
                answer(i, j) = E(i, j);
            }
        }
    return answer;
    }
};

#endif // LINAL_HEADER_H_INCLUDED
