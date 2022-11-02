#include<iostream>
#include <cmath>
#include<complex> //Standart c++ library for complex numbers.
#include <iomanip>

const int L = 4; //space size
const double theta = M_PI/4;
std::complex<double> p(std::cos(theta), 0); //Transition amplitudes
std::complex<double> q(0, std::sin(theta));
const int Q = 2; //Number of directions


class QLB{
private:
  std::complex<double> Psi[Q*L]; //Wave vectors are created as private attributes
  std::complex<double> Psi_new[Q*L];

 

public:
  void Start(void); //This impose the initial conditions.
  void Get_Psi(void); //Returns the Wave vector.
  double rho(int ix); //Returns the probability density at each cell.
  void Collision(void); //Creates the Collsion operator and modifies Psi_new.
  void Advection(void); // Creates the Advection operator and modifies Psi.
  void Show_matrix(std::complex<double> M); //Print a given matrix.
  


};


void QLB::Start(void){
  std::complex<double> z(1,1);
  for(int ix=0; ix<Q*L; ix++)
    Psi[ix] = z;

}

void QLB::Get_Psi(void){
 for(int ix=0; ix<Q*L; ix++)
   std::cout<<Psi[ix]<<std::endl;

}

double QLB::rho(int ix){

  if (ix%2==0){return std::norm(Psi[ix]+Psi[ix+1]);}
  else{return std::norm(Psi[ix]+Psi[ix-1]);}


}

void QLB::Collision(void){
 std::complex<double> C[Q*L][Q*L];
 for (int i=0; i<Q*L; i++)
   for(int j=0; j<Q*L; j++){
     if(i==j and j%2==0) {C[i][j]=p; C[i+1][j]=q;}
      else if (i==j and j%2==1) {C[i][j]=p; C[i-1][j]=q;}
     }

 for (int i=0; i<Q*L; i++)
   for(int j=0; j<Q*L; j++){
     Psi_new[i] += C[i][j]*Psi[j];
   }
  
}

void QLB::Advection(void){
 std::complex<double> M[Q*L][Q*L];
 for(int j=0; j<Q*L; j++){
   if(j%2==0){M[(j+2+Q*L)%(Q*L)][j]=(1,1);}
   else if(j%2==1){M[(j+2*L-2+Q*L)%(Q*L)][j]=(1,1);}
 }

 for (int i=0; i<Q*L; i++)
   for(int j=0; j<Q*L; j++){
     Psi[i] += M[i][j]*Psi_new[j];
   }
 

}

void QLB::Show_matrix(std::complex<double> M){

  for(int j=0; j<Q*L; j++){
   for(int i=0; i<Q*L; i++)
     std::cout<<M[i][j]<<" ";
     std::cout<<std::endl; }
     

}



int main(){
  std::cout << std::fixed << std::setprecision(1); //This is to choose the precision of complex numbers.
    




     
 
 

  return 0;
}
