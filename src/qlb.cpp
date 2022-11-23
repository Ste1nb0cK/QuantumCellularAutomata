#include <cmath>
#include <complex> //Standard c++ library for complex numbers.
#include <iomanip>
#include <iostream>
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>
//Trying to implement the Harmonic Oscillator
//---------------------Typdefs--------------------------------------------------
typedef std::complex<double> complex;
typedef Eigen::VectorXcd Vector;
typedef Eigen::MatrixXcd Matrix;
//----------------------Simulation Conditions-----------------------------------
const int L = 3; // space size. For Symmetric V we should use only odd L.
const double theta = M_PI / 4;
const complex p(std::cos(theta), 0); // Transition amplitudes
const complex q(0, std::sin(theta));
const int Q = 2; // Number of directions
const int N = Q*L; //Dimension of the vectors we will be using
//-----------------------System Conditions--------------------------------------
//The system is asumed to be in natural units, for the SHO case:
//For this case choice of natural units for the system makes the potential look:
// V(x) = x^2
// the mass, the frecuency of the oscillator  and \hbar\omega are equal to unity.
// The lattice units are simply normalizations of this ones by N or tmax. This
// is relevant to apply the Schrödinger Operator. As the center of our
// simulation is L/2 we take that as reference i.e. the potential is
// V(x) = (x-L/2)^2
double Potential(double x);
class QLB {
  private:
  Vector Psi ; // Wave vectors are created as private attributes.
  Vector Psi_new;
  Matrix M; //TODO: Implement this using sparse and diagonal matrices.
  Matrix C;
  Matrix V;
  Matrix C2;

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
  void PrintM(void); //for testing
  void Evolution(void); //Apply a complete evolution step
  Matrix Get_Evolution(void);
};

Matrix QLB::Get_Evolution(void){




  
 
  return (V*(M*C2)*(M*C)).pow(1000);
}



void QLB::PrintM(void){



  std::cout<<V;
}


QLB::QLB(void){
  int i,j;
  //Declare the size of the matrices, otherwise this wont work.
  Psi.resize(N);
  Psi_new.resize(N);
  M.resize(N,N);
  C.resize(N,N);
  V.resize(N,N);
  //Set vectors to zero
  Psi = Vector::Zero(N);
  Psi_new = Vector::Zero(N);
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
  for (j = 0; j <  N; j++) {
    if (j % 2 == 0) {
      M((j + 2 +  N) % (N),j) = (1, 1);
    } else if (j % 2 == 1) {
      M((j + 2 * L - 2 + N) % (N),j) = (1, 1);
    }
  }

  double epsilon = 1;
  V = Matrix::Zero(N,N);
  for (i=0; i<L; i++){
    double V_x = Potential(i);
    complex aux(std::cos(epsilon*V_x), -std::sin(epsilon*V_x));
    V(2*i+1,2*i+1) = V(2*i,2*i) = aux;

  }

  C2.resize(N,N);
    for(int j =0; j<N; j++){
    for(int i =0; i<N; i++){
      C2(j,i)=C((i-1+N)%N, j);

    }
    }
}
void QLB::Start(void) {
  //TODO: Pass the initial conditions more generically
  // A right traveling plane wave is created as  initial condition.
  double k = (2 * M_PI / L);
  double mu = L/2;
  double sigma2 = L*L/(100);
  for (int ix = 0; ix < Q * L; ix++) {
    if (ix % 2 == 1) {
      complex z(std::cos(k * ix), -1 * std::sin(k * ix));
      Psi(ix) = z*std::exp(-std::pow(ix-mu, 2)/(2*sigma2));
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

void QLB::Evolution(void){
  Psi = V*M*C*Psi;
} //Apply a complete evolution step
void QLB::Print_Rho(void) {
  for (int ix = 0; ix < N; ix += 2)
    std::cout << ix / 2 << " " << std::norm(Rho(ix)) << std::endl;
}
int main() {
  std::cout << std::fixed
            << std::setprecision(
				 6); // This is to choose the precision of complex numbers.
  QLB particle;

  particle.Start();

  Matrix U = particle.Get_Evolution();
  Eigen::ComplexEigenSolver<Matrix> Solver(U, false);
  Solver.compute(U);

  std::cout<<Solver.eigenvectors();
  Vector Psi2 = Solver.eigenvectors().col(0
					  );
  //std::cout<<Solver.eigenvalues()[0]*Psi2 - U*Psi2<<std::endl;
 



  Vector Psi;
  Psi.resize(L);
  


  

   for(int i = 0; i<L; i++){
    Psi(i) = Psi2(2*i)+Psi2(2*i+1);

  }
 


   /*for (int i =0; i<L; i++){
    std::cout<<i<<" "<<std::norm(Psi(i))<<std::endl;
   
    }*/
   


     // particle.PrintM();
 
 


  return 0;

  }
  







double Potential(double x){
  return 0.5*std::pow((x-(L/2)), 2);

  

  
}


