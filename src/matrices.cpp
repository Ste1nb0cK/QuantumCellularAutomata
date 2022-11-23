#include <cmath>
#include <complex> //Standard c++ library for complex numbers.
#include <iomanip>
#include <iostream>
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>

typedef std::complex<double> complex;
typedef Eigen::VectorXcd Vector;
typedef Eigen::MatrixXcd Matrix;

const int L = 9; // space size. For Symmetric V we should use only odd L.
const double theta = M_PI / 4;
const complex p(std::cos(theta), 0); // Transition amplitudes
const complex q(0, std::sin(theta));
const int Q = 2; // Number of directions
const int N = Q*L; //Dimension of the vectors we will be using

double Potential(double x);

class Evolution{

private:
  Vector Psi;
  Matrix M;
  Matrix C;
  Matrix V;

public:
  Evolution(void);
  void PrintMatrix(void);
  Matrix GetMatrix(void);
  
  



};

Evolution::Evolution(void){

 M.resize(N,N);
 C.resize(N,N);
 V.resize(N,N);
 int i, j;

 //Fill the colision matrix
 for (j = 0; j<N ; j++){
    //divide into odd and even case
    if(j%2 ==0){
      for(i=0; i<N; i++){
        if(i==j){ C(i,j)=p;}
        else if(i-j == 1){C(i,j)=q;}
        else {C(i,j)=0;}
      }
    }
    if(j%2!=0){
      for(i=0; i<N; i++){
        if(i==j){ C(i,j)=p;}
        else if(j-i == 1){C(i,j)=q;}
        else {C(i,j)=0;}
      }
    }
  }
 //Fill the advection matrix
 
 M = Matrix::Zero(N,N);
  for (j = 0; j <  N; j++) {
    if (j % 2 == 0) {
      M((j + 2 +  N) % (N),j) = (1, 1);
    } else if (j % 2 == 1) {
      M((j + 2 * L - 2 + N) % (N),j) = (1, 1);
    }
  }
  //Fill the potential matrix
 double epsilon = 1;
  V = Matrix::Zero(N,N);
  for (i=0; i<L; i++){
    double V_x = Potential(i);
    complex aux(std::cos(epsilon*V_x), -std::sin(epsilon*V_x));
    V(2*i+1,2*i+1) = V(2*i,2*i) = aux;

  }


}

void Evolution::PrintMatrix(void){

  std::cout<<V;


}


Matrix Evolution::GetMatrix(void){

  return M*V*C;

}


int main(){

 std::cout << std::fixed<< std::setprecision(2);

  Evolution QHO;
  // QHO.PrintMatrix();

  //Declare the evolution matrix and compute its eigenvalues and eigen vectors.
  
  Matrix U = QHO.GetMatrix();
  Eigen::ComplexEigenSolver<Matrix> Solver(U, false);
  Solver.compute(U);
  // std::cout<<Solver.eigenvectors();


  Vector Psi2 = Solver.eigenvectors().col(4);

  Vector Psi;
  Psi.resize(L);

  //Creating the wave vector in the position basis. 

  for(int i = 0; i<L; i++){
    Psi(i) = Psi2(2*i)+Psi2(2*i+1);

  }

  for (int i =0; i<L; i++){
    std::cout<<i<<" "<<std::norm(Psi(i))<<std::endl;
   
    }
  

  
  



 return 0;


}




double Potential(double x){
  return 0.5*std::pow((x-(L/2)), 2);
}
