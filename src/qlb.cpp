#include <cmath>
#include <complex> //Standard c++ library for complex numbers.
#include <iomanip>
#include <iostream>
#include <Eigen/Dense>

//alias for complex numbers
typedef std::complex<double> complex;

const int L = 500; // space size
const double theta = M_PI / 4;
complex p(std::cos(theta), 0); // Transition amplitudes
complex q(0, std::sin(theta));
const int Q = 2; // Number of directions
const int N = Q*L; //Dimension of the vectors we will be using

// typedef ('aliases')
typedef Eigen::VectorXcd Vector;
typedef Eigen::MatrixXcd Matrix;
double Potential(double x);


class QLB {
private:
  Vector Psi ; // Wave vectors are created as private attributes.
  Vector Phi;
  Matrix M;
  Matrix C;
  Matrix V;
  Matrix C2;
  

public:
  QLB(void); //constructor. Initialize state as zero
  void Start(void);     // This imposes the initial conditions.
  void Get_Psi(void);   // Returns the Wave vector. Just for test.
  complex Rho(int ix);  // Returns the probability density at each cell.
  void Print_Rho(void); // Prints the real part fo the wave function, la                       
  void Evolution(void);
  void DFT(void);
  void Print_Rho_Moment(void);
  double Variance(Vector V);
  double Uncertainty();
  
};

QLB::QLB(void){
  int i,j;
  //Declare the size of the matrices, otherwise this wont work.
  Psi.resize(N);
  M.resize(N,N);
  C.resize(N,N);
  Phi.resize(L);
  for (i=0; i<N; i++){
    Psi(i) = (0,0);
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
 
  //----------------------------Right Traveling Planewave---------------------//
   
  /*  double k = (2 * M_PI / L);
   for (int ix = 0; ix < N; ix++) {
     if(ix%2==1){
      complex z(std::cos(k * ix), -1 * std::sin(k * ix));
       Psi(ix) = z;}
     

     else{ Psi(ix)=0;}
     } */
   
  //----------------------------Gaussian--------------------------------------//

  double k = ( 3*2*M_PI / L);
  // double k = 0;
  double mu = L/2;
  double sigma2 = (L*L)/(100);
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


complex QLB::Rho(int ix) { return Psi(ix,0) + Psi(ix + 1,0); }


void QLB::Print_Rho(void) {
  for (int ix = 0; ix < N; ix += 2)
    std::cout << ix/2 << " " << std::norm(Rho(ix)) << std::endl;
    //Add two blank lines for animating in gnuplot
    std::cout << "\n"<< "\n";
}

void QLB::Evolution(void){

  Psi =(V*M*C)*Psi;

}


void QLB::DFT(void){

  complex lenght(1/L,0);
  Vector AUX;
  AUX.resize(L);
  for(int ix = 0; ix<L; ix++){

    AUX(ix) = Psi(2*ix) + Psi(2*ix +1);
  }
  
 for(int k=0; k<L; k++){
   complex sum;
   sum = 0;
   for(int n=0; n<L; n++){
     complex phase(cos((2*M_PI/L)*k*n), -sin((2*M_PI/L)*k*n)); 
     sum += AUX(n)*phase; 

   }
   Phi[k]=sum;
  }
}

double QLB::Variance(Vector V){

  //Calcular N
  double N = 0;
  for(int ix = 0; ix<L; ix++){

    N+= std::norm(V(ix));

  }
  double prom = 0;
   //Calcular prom
    for(int ix = 0; ix<L; ix++){

   prom += std::norm(V(ix))*ix;

  }
    
  prom/=N;

  //Calcular Sigma2
  double Sigma2 = 0;
  for(int ix = 0; ix<L; ix++)
   Sigma2+=pow(ix-prom,2.0)*std::norm(V(ix));


 Sigma2/=N;
 return prom;
}

double QLB::Uncertainty(void){

  Vector AUX;
  AUX.resize(L);
  for(int ix = 0; ix<L; ix++){

    AUX(ix) = Psi(2*ix) + Psi(2*ix +1);
  }
  
  // return std::sqrt(Variance(AUX))*std::sqrt(Variance(Phi));
  return Variance(AUX);
  

}
void QLB::Print_Rho_Moment(void){
 for (int ix = 0; ix < L; ix ++)
   std::cout << ix << " " << std::norm(Phi((ix -30 + L)%L)) << std::endl;
    //Add two blank lines for animating in gnuplot
    std::cout << "\n"<< "\n";

}

int main() {
  std::cout << std::fixed
            << std::setprecision(
                   16); // This is to choose the precision of complex numbers.
  QLB free_particle;
  free_particle.Start();
  //free_particle.Get_Psi();


   for (int t = 0; t < 100; t+=1) {
  
    free_particle.Evolution();
    // free_particle.DFT();
    //free_particle.Print_Rho_Moment();
    // free_particle.Print_Rho();
    std::cout<<t<<" "<<free_particle.Uncertainty()<<std::endl;
    //free_particle.Get_Psi();
    }

  return 0;
}

double Potential(double x){


  /* if (x<L/2) return 0;
     else{return 10;}*/

  /* if (x == 0 or x == L) {return 1e10;}
     else {return 0;}*/



  // return 30*0.5*std::pow(x-(L/4), 2);

  return 0;
  

}
