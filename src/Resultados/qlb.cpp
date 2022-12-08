#include <cmath>
#include <complex> //Standard c++ library for complex numbers.
#include <iomanip>
#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>
//Trying to implement the Harmonic Oscillator
//---------------------Typdefs--------------------------------------------------
typedef std::complex<double> complex;
typedef Eigen::VectorXcd Vector;
typedef Eigen::MatrixXcd Matrix;
typedef std::fstream stream;
//----------------------Simulation Conditions-----------------------------------
const int L = 81*4; // space size. For Symmetric V we should use only odd L.
const double sigma2 = L*L/900.0;
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

//--------------------FUNCIONES FUERA DE LA CLASE-------------------
double Potential(double x);
void plotInitilizeMap3D(stream & Map, const char * File);
void plotInitilizeAnimation2D(stream & Animation, const char * File);
void plotCloseMap3D(stream & Map);
//------------------------------  CLASE --------------------
class QLB {
  private:
  Vector Psi ; // Wave vectors are created as private attributes.
  Vector Psi_new;
  Vector Psi_PositionBase;
  Matrix M; //TODO: Implement this using sparse and diagonal matrices.
  Matrix C;
  Matrix V;
  Matrix C2;

public:
  QLB(void); //constructor. Initialize state as zero
  void Start(void);     // This imposes the initial conditions.
  void Get_Psi(void);   // Returns the Wave vector. Just for test.
  void Get_Psi_new(void);   // Returns the Wave vector. Just for test.
  complex Rho(int ix);  // Returns the probability density at each cell.
  void Collision(void); // Creates the Collsion operator and modifies Psi_new.
  void Advection(void); // Creates the Advection operator and modifies Psi.
  void Print_Rho(void); // Prints the real part fo the wave function, later, the
                        // function will pirnt the probability density.
  double Varianza(void);
  void PrintM(void); //for testing
  void plotPerfil(stream & Map, int t);
  void plotFrame2D(stream & Animation, int t);
  void Evolution(void); //Apply a complete evolution step
  void Update_Psi_PositionBase(void);
  Matrix Get_Evolution(void);
};


int main() {
  int tmax = 700;
  std::cout << std::fixed
            << std::setprecision(
				 6); // This is to choose the precision of complex numbers.
  QLB particle;

  //The streams that send data are opened
   stream Animation;
  stream Map;
  const char * AnimationFile = "Animation.p";
  const char * MapFile = "Map.p";
  Animation.open(AnimationFile, std::ios_base::out);
  Map.open(MapFile, std::ios_base::out);
  Animation<<std::setprecision(8);
  Map<<std::setprecision(11);
  
  plotInitilizeMap3D(Map, MapFile);
  plotInitilizeAnimation2D(Animation, AnimationFile);
  
  //Run
  particle.Start();
  particle.Update_Psi_PositionBase();
  
  particle.plotFrame2D(Animation,0);
  particle.plotPerfil(Map,0);
  
  for(int t = 0; t<tmax; t++){
    particle.Evolution();
    if( t % 25 == 0){
      particle.plotFrame2D(Animation, t + 1);
      particle.plotPerfil(Map, t+1);
    }
  }
  
  plotCloseMap3D(Map);
  
  Map.close();
  Animation.close();
  return 0;

}

//-------------------FUNCIONES FUERA DE CLASE------------------
double Potential(double x){
  double Oscilador = 0.5*std::pow((x-(L/2)), 2);
  return 0.0;
}


void plotInitilizeMap3D(stream & Map, const char * File){
  Map<<"reset session"<<std::endl;		     
  Map<<"set xlabel 'x [Cells]' "<<std::endl;	     
  Map<<"set ylabel 'Time [clicks]' "<<std::endl;	     
  Map<<"set zlabel 'Re[{/Symbol Y}(x)]' "<<std::endl;
  Map<<"unset key"<<std::endl;
  Map<<"set pm3d"<<std::endl;
  Map<<"##set terminal jpeg enhanced"<<std::endl;
  Map<<"##set out '"<<File<<".jpg' "<<std::endl;
  Map<<"set terminal qt"<<std::endl;
  Map<<"unset surface"<<std::endl;
}
void plotInitilizeAnimation2D(stream & Animation, const char * File){
  Animation<<"reset session"<<std::endl;		     
  Animation<<"set xlabel 'x [Cells]' "<<std::endl;	     
  Animation<<"set ylabel 'Re[{/Symbol Y}(x)]' "<<std::endl;
  Animation<<"set key box "<<std::endl;
  Animation<<"set terminal gif animate"<<std::endl;
  Animation<<"set out '" <<File<<".gif'  "<<std::endl;
  Animation<<"set yrange [-0.01:"<< 1.0/std::pow((81.0*81.0/900.0*2*M_PI), 0.5) <<"]"<<std::endl;
}

void plotCloseMap3D(stream & Map){
  Map<<"e"<<std::endl;
}

//-------------------FUNCIONES DE CLASE-------------------
QLB::QLB(void){
  int i,j;
  //Declare the size of the matrices, otherwise this wont work.
  Psi.resize(N);
  Psi_new.resize(N);
  Psi_PositionBase.resize(L);
  M.resize(N,N);
  C.resize(N,N);
  V.resize(N,N);
  //Set vectors to zero
  Psi = Vector::Zero(N);
  Psi_new = Vector::Zero(N);
  Psi_PositionBase = Vector::Zero(L);
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

Matrix QLB::Get_Evolution(void){
  return (V*(M*C2)*(M*C)).pow(1000);
}
void QLB::Update_Psi_PositionBase(void){
  for(int i = 0; i < L; i++){
    Psi_PositionBase(i,0) = Psi_new(2*i,0) + Psi_new(2*i +1,0);
  }
}

void QLB::PrintM(void){
  std::cout<<V;
}

void QLB::Start(void){
  double sigma2 = 81.0*81.0/900.0;
  double k = ( 2* M_PI / L);
  double mu = 81.0/2.0;
  for (int ix = 0; ix < Q * L; ix++) {
    if (ix % 2 == 1) {
      Psi(ix) =  std::exp(-std::pow(ix/4.0/2.0 - mu, 2)/(4*sigma2))/std::pow(sigma2*2*M_PI, 0.25);
    }
  }
}
      //complex z(std::cos(k * ix)/2.0, -1 * std::sin(k * ix)/2.0);
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

double QLB::Varianza(void){
  double Prom = 0.0;
  double N = 0.0;
  double Var = 0.0;

  for(int ix = 0; ix < N; ix += 2){
    N += std::norm(Rho(ix));
  }
  std::cout<<N<<"\t";
  
  for(int ix = 0; ix < N; ix += 2){
    Prom +=  std::norm(Rho(ix))*(ix/2.0);
  }
  Prom /= N;
  
  std::cout<<Prom<<"\t";
  for(int ix = 0; ix < N; ix += 2){
    Var += std::norm(Rho(ix))*(ix/2.0 - Prom)*(ix/2.0 - Prom);
  }
  Var /= N;
  
  return Var;
}

void QLB::plotFrame2D(stream & Animation, int t){
  Animation<<"plot '-' t 't = "<<t<<"' w l"<<std::endl;
  for(int x = 0; x < N; x += 2){
    Animation<< x/2 << "\t" << std::norm(Rho(x))<<std::endl;
  }
  Animation<<"e"<<std::endl;
}

void QLB::plotPerfil(stream & Map, int t){

  Map<<"splot '-' w l "<<std::endl;
  for(int x = 0; x < N; x += 2 ){
    Map<< x/2 <<"\t"<< t <<"\t"<< std::norm(Rho(x))<<std::endl;
  }
  Map<<"\n";
}
