#include <cmath>
#include <complex> //Standard c++ library for complex numbers.
#include <iomanip>
#include <iostream>
#include <Eigen/Dense>

//alias for complex numbers
typedef std::complex<double> complex;

const int L = 300; // space size
const double theta = M_PI / 4;
complex p(std::cos(theta), 0); // Transition amplitudes
complex q(0, std::sin(theta));
const int Q = 2; // Number of directions
const int N = Q*L; //Dimension of the vectors we will be using

// typedef ('aliases')
typedef Eigen::VectorXcd Vector;
typedef Eigen::MatrixXcd Matrix;


class QLB {
private:
  Vector Psi ; // Wave vectors are created as private attributes.
  Vector Psi_new;
  Matrix M;
  Matrix C;

public:
  QLB(void); //constructor. Initialize state as zero
  void Start(void);     // This imposes the initial conditions.
  void Get_Psi(void);   // Returns the Wave vector. Just for test.
  void Get_Psi_new();   // Returns the Wave vector. Just for test.
  complex Rho(int ix);  // Returns the probability density at each cell.
  void Collision(void); // Creates the Collsion operator and modifies Psi_new.
  void Advection(void); // Creates the Advection operator and modifies Psi.
  void Print_Rho(void); // Prints the real part fo the wave function, later, the
                        // function will pirnt the probability density.
  // void PrintM(void); //for testing
};

// void QLB::PrintM(void){
//   std::cout << M << std::endl;
// }
QLB::QLB(void){
  int i,j;
  //Declare the size of the matrices, otherwise this wont work.
  Psi.resize(N);
  Psi_new.resize(N);
  M.resize(N,N);
  C.resize(N,N);
  for (i=0; i<N; i++){
    Psi(i) = (0,0);
    Psi_new(i) = (0,0);
  }
  // initialize the Collision Catrix
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
 // initialize the Advection
  M = Matrix::Zero(N,N);
  for (int j = 0; j <  N; j++) {
    if (j % 2 == 0) {
      M((j + 2 +  N) % (N),j) = (1, 1);
    } else if (j % 2 == 1) {
      M((j + 2 * L - 2 + N) % (N),j) = (1, 1);
    }
}
}
void QLB::Start(void) {
  // A right traveling plane wave is created as  initial condition.
  complex z;
  double k = (2 * M_PI / L);
  for (int ix = 0; ix < N; ix++) {
    if (ix % 2 == 1) {
      z = (std::cos(k * ix), -1 * std::sin(k * ix));
      Psi(ix) = z;
    }
  }
}

void QLB::Get_Psi(void) {
  for (int ix = 0; ix < N; ix++)
    std::cout << Psi(ix,0) << std::endl;
}

void QLB::Get_Psi_new(void) {
  for (int ix = 0; ix < N; ix++)
    std::cout << Psi_new(ix,0) << std::endl;
}
complex QLB::Rho(int ix) { return Psi(ix,0) + Psi(ix + 1,0); }

void QLB::Collision(void) {
  Psi_new = C * Psi;
}

void QLB::Advection(void) {
  Psi = M * Psi_new;
}

void QLB::Print_Rho(void) {
  for (int ix = 0; ix < N; ix += 2)
    std::cout << ix / 2 << " " << std::real(Rho(ix)) << std::endl;
}

int main() {
  std::cout << std::fixed
            << std::setprecision(
                   3); // This is to choose the precision of complex numbers.
  QLB free_particle;
  free_particle.Start();
  for (int t = 0; t < 200; t++) {
    free_particle.Collision();
    free_particle.Advection();
  }

  free_particle.Print_Rho();
  return 0;
}
