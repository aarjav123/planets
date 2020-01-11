// Planet_leapfrog_submit.cc
//  Simulates planet (inverse square law)
//    using leapfrog method.
//  example usage:
//       mkdir output
//       ./Planet > output/1
// recommended gnuplot commands:
//  set size ratio -1
//  plot   'output/1' u 2:3 w linesp lt 3 pt 4
//  replot 'output/1' u 2:3:4:5 w vector lt 5

#include <iostream>
#include <cmath>
using namespace std;

#define D  2  // number of dimensions
struct particle {
  double x[D] ; // (x,y) coordinates
  double p[D] ; // momentum
  double F[D] ; // force
  double im   ; // inverse mass
  double GMm  ; // gravitational parameter of this particle
  double v[D] ; // velocity
  double T    ; // kinetic energy
  double V    ; // potential energy
  double r    ; // absolute distance from origin
} ; // Note the definition of a structure ends with a semicolon

struct control {
  int verbose ; // program verbosity
  int printing ; // period with which to print
};

void Force( particle &a )
  // sets the force vector and the potential energy
{
  double R = 0.0 ; 
  for ( int i = 0 ; i < D ; i++ ) {
    R += a.x[i] * a.x[i] ;
  }
  double r = sqrt(R) ;
  a.V = - a.GMm / r ;
  a.r = r ; 
  for ( int i = 0 ; i < D ; i++ ) {
    a.F[i] = - a.GMm * a.x[i] / (r*R) ;  // inverse sq
  }
}

void PositionStep ( particle &a , double dt )
{ // the "&" indicates 'this is a pass by reference'
  for ( int i = 0 ; i < D ; i++ ) { 
    a.x[i] += dt * a.p[i] * a.im ; 
  } 
}

void MomentumStep ( particle &a , double dt )
{ // the "&" indicates 'this is a pass by reference. 
  for ( int i = 0 ; i < D ; i++ ) { 
    a.p[i] += dt * a.F[i] ;
  }
}

void v2p( particle &a )
  // propagate a change in velocity into the momentum vector
{
  for ( int i = 0 ; i < D ; i++ ) {
    a.p[i] =  a.v[i] / a.im ;
  }
}

void p2v ( particle &a ) {
  for ( int i = 0 ; i < D ; i++ ) { 
    a.v[i]  = a.p[i] * a.im ;
  }
}

void pv2T ( particle &a ) {
  a.T=0.0;
  for ( int i = 0 ; i < D ; i++ ) { 
    a.T    += 0.5*a.v[i] * a.p[i] ; 
  }
}

void showState ( particle &a )
{ 
  for ( int i = 0 ; i < D ; i++ ) { 
    cout << "\t"<<a.x[i];
  }
  for ( int i = 0 ; i < D ; i++ ) { 
    cout << "\t"<<a.v[i];
  }
  cout << "\t" << a.T << "\t" << a.V << "\t" << a.T+a.V << endl;
  cout << "#" ;
  for ( int i = 0 ; i < D ; i++ ) { 
    cout << "\tx["<<i<<"]\t";
  }
  for ( int i = 0 ; i < D ; i++ ) { 
    cout << "\tv["<<i<<"]\t";
  }
  cout << "\tT\t\tV\t\tT+V" ; 
  cout << endl;
}

void LeapfrogDynamics( particle &a , 
		       double dt , // step size
		       double &t , // current time
		       int N ,  // number of steps
		       control &c  // parameters controlling output, etc
		       )
{
  for ( int i=0 ; i < N ; i ++ ) {
    if( c.printing && !(i%c.printing) ) {
      // this iteration we will print out the state
      p2v(a);   pv2T(a) ;     Force( a ) ;
      cout << t ; showState( a ) ;
    }
    // each iteration we move the position a half-step
    PositionStep( a , dt*0.5 ) ;

    t += dt*0.5 ; //INCREMENT t BY HALF dt AS SUGGESTED IN THE NOTES

    Force( a ) ; // (computes the force at this position)
    // then we move the momentum full-step
    MomentumStep( a , dt )     ;
    // then we move the position a half-step
    PositionStep( a , dt*0.5 ) ;
    t += dt*0.5 ; //IINCREMENT t BY HALF dt (NOT dt AS DONE ORIGINALLY), FOLLOWING THE NOTES.
  }
}

void showVelocity ( particle &a )
{ // the "&" indicates 'this is a pass by reference. 
  for ( int i = 0 ; i < D ; i++ ) { 
    cout << "v["<<i<<"] = "<<a.v[i]<<endl; 
  }
}


int main()
{
  particle   a ;
  control    c ;
  c.verbose=1;
  c.printing=5 ; // how many iterations to go
  //               between state-printings

  a.im   = 1.0 ;
  a.v[0] = 0.0;
  a.v[1] = 1.0/3.0;
  a.x[0] = 9.0;
  a.x[1] = 0;

  a.GMm = 1.0;
  v2p(a) ; // make sure the momentum is set up

  double dt = 0.001 ; // time per step
  double t = 0.0 ; // start time
  int N = 500000 ;  // number of steps to make

  LeapfrogDynamics( a , dt , t , N , c ) ; 
  
  return 0;
}
