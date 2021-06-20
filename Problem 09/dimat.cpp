/*SE solver with the 'IMPLICIT METHOS' - Diego Andreoni 2021*/

#include <iostream>
#include <Eigen/Eigenvalues>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <fstream>
#include <cmath>
#include <sstream>
#include <string>
#include <iomanip>



using namespace std;
using namespace Eigen;

int main(){
	//Paramters of the problem
	int N=3000; //Number of points of the mesh
	double rmin=0.4, rmax=5; //Interval
	double h; //Distance between the poin of the mesh

	h=(rmax-rmin)/N;
	double f=1/(h*h);
	/*The problem we are approaching -0.5Au+BVu=EBu, multiplying on the LHS
	by inverse of B we get: -0.5(Inv(B)*A+V=EIdentity. We have a nice eigenvalues problem*/

	//Define Matrix A
	MatrixXd A(N,N);
	A.setZero();
	A.diagonal(1).setConstant(1);
	A.diagonal(-1).setConstant(1);
	A.diagonal(0).setConstant(-2);
	A=f*A;

	//Define Matrix B
	MatrixXd B(N,N);
	B.setZero();
	B.diagonal(1).setConstant(1);
	B.diagonal(0).setConstant(10);
	B.diagonal(-1).setConstant(1);
	B=0.0833333333*B;

	//Invertiamo la matrice B
	MatrixXd InvB(N,N);
	InvB=B.inverse();

	//Define Matrix of the potential
	MatrixXd V(N,N);
	V.setZero();
	double p;
	for(int i=0; i<N; i++){
		p=(1/(rmin+h*i))*(7.39*exp(-3.11*(rmin+h*i))-3.22*exp(-1.555*(rmin+i*h)));	//Potenziale di Coulomb
		//p=-1/(rmin+i*h);
		V(i,i)=p;
	}
	//Matrix of the Hamiltonian
	MatrixXd H(N,N);
	H.setZero();
	H=-0.5*InvB*A+V;

//------>RESULTS:
/*Nota, dato che eigen mette a caso gli autovalori, fino a che non faccio il Bubble sort per ordinarli, devo selezionare gli stato a mano*/
	//Stampiamo autovalori della matrice
	SelfAdjointEigenSolver <MatrixXd> es;
	es.compute(H);

	for(int i=0; i<N; i++){
		if(es.eigenvalues()[i]<0){
		cout<<es.eigenvalues()[i]<<endl;
		}
	}

}
