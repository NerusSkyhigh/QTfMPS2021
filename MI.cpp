//Implicit method

#include <iostream>
#include <Eigen/Eigenvalues>
#include <Eigen/Core>

using namespace std;
using namespace Eigen;

int main(){
	int N=20; //Numero di punti
	double rmin=-1, rmax=1;
	double h;
	double U[N][N];

	h=(rmax-rmin)/N;

	//Definiamo la matrice
//	for(int i=0; i<N; i++){
//		for(int j=0; j<N; j++){
//			if(i==j){
//			U[i][j]=-10;
//			}
//			else if(i==j+1){
//			U[i][j]=1;
//			}
//			else if(i==j-1){
//			U[i][j]=1;
//			}
//		}
//	}

	MatrixXd A(N,N);
	A.setZero();
	A.diagonal(1).setConstant(-1/(h*h));
	A.diagonal(-1).setConstant(-1/(h*h));
	A.diagonal(0).setConstant((2/h)-70);
	cout<<A<<endl<<endl;

	EigenSolver <MatrixXd> es;
	es.compute(A);

	cout<<es.eigenvalues()<<endl;
	//cout<<es.eigenvectors().col(1)<<endl;

	A(1,1)=100000000;

	cout<<A(1,1);
}

