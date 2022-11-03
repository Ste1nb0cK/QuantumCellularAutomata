#include<iostream>
#include <cmath>
#include<complex> //Standart c++ library for complex numbers.
#include <iomanip>

const int L = 300; //space size
const double theta = M_PI/4;
std::complex<double> p(std::cos(theta), 0); //Transition amplitudes
std::complex<double> q(0, std::sin(theta));
const int Q = 2; //Number of directions


class QLB{
private:
  std::complex<double> Psi[Q*L]; //Wave vectors are created as private attributes.
  std::complex<double> Psi_new[Q*L];

 

public:
  void Start(void); //This impose the initial conditions.
  void Get_Psi(void); //Returns the Wave vector. Just for test.
  void Get_Psi_new();//Returns the Wave vector. Just for test.
  std::complex<double> Rho(int ix); //Returns the probability density at each cell.
  void Collision(void); //Creates the Collsion operator and modifies Psi_new.
  void Advection(void); // Creates the Advection operator and modifies Psi.
  void Print_Rho(void); //Prints the real part fo the wave function, later, the function will pirnt the probability density. 
  


};


void QLB::Start(void){
  //A right traveling plane wave is created as a initial condition.
  double k = (2*M_PI/L); 
  for(int ix=0; ix<Q*L; ix++){
    if (ix%2==1){
    std::complex<double> z(std::cos(k*ix),-1*std::sin(k*ix));
    Psi[ix]=z;
  }
    

 }
}

void QLB::Get_Psi(void){
 for(int ix=0; ix<Q*L; ix++)
   std::cout<<Psi[ix]<<std::endl;

}

void QLB::Get_Psi_new(void){
 for(int ix=0; ix<Q*L; ix++)
   std::cout<<Psi_new[ix]<<std::endl;

}
std::complex<double> QLB::Rho(int ix){

  return Psi[ix] + Psi[ix+1];

}

void QLB::Collision(void){
 std::complex<double> C[Q*L][Q*L];
 for (int i=0; i<Q*L; i++)
   for(int j=0; j<Q*L; j++){
     if(i==j and j%2==0) {C[i][j]=p; C[i+1][j]=q;}
      else if (i==j and j%2==1) {C[i][j]=p; C[i-1][j]=q;}
     }

 for(int i=0; i<Q*L; i++){Psi_new[i]=(0,0);}

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

 

  for(int i=0; i<Q*L; i++){Psi[i]=(0,0);}

  for (int i=0; i<Q*L; i++)
   for(int j=0; j<Q*L; j++){
     Psi[i] += M[i][j]*Psi_new[j];
   }
 

}

void QLB::Print_Rho(void){

  for(int ix=0;ix<L*Q;ix+=2)
    std::cout<<ix/2<<" "<<std::real(Rho(ix))<<std::endl;


}


  

     



int main(){
  std::cout << std::fixed << std::setprecision(3); //This is to choose the precision of complex numbers.


  QLB free_particle;

  free_particle.Start();
  

  for(int t=0; t<200; t++){

    free_particle.Collision();
    free_particle.Advection();

  }

  free_particle.Print_Rho();


  
 
 
  //Piece of code for print matrices.
  /* for(int j=0; j<Q*L; j++){
   for(int i=0; i<Q*L; i++)
     std::cout<<M[i][j]<<" ";
     std::cout<<std::endl; }*/
     

     

  


 
 

  

     
 
 

  return 0;

}
