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
	int N=1000; //Number of points of the mesh
	double rmin=-5, rmax=5; //Interval
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
	for(int i=0; i<N; i++){
		double p=0;
		V(i,i)=p;
	}
	//Matrix of the Hamiltonian
	MatrixXd H(N,N);
	H.setZero();
	H=-0.5*InvB*A+V;

//------>RESULTS:
/*Nota, dato che eigen mette a caso gli autovalori, fino a che non faccio il Bubble sort per ordinarli, devo selezionare gli stato a mano*/
	//Stampiamo autovalori della matrice
	SelfAdjointEigenSolver <MatrixXd> es(H);

//----->BUBBLE SORT VETTORE AUTOVALORI (che la libreria Eigen non restituisce in ordine)
	double eigenval[N];

	//Generiamo un vettore che per elementi gli autovalori senza ordine
	ofstream E;
	E<<es.eigenvalues();

	for(int i=0; i<N; i++){
//		eigenval[i]=es.eigenvalues().real()[i];
		eigenval[i]=es.eigenvalues()[i];
//		cout<<eigenval[i]<<" "<<es.eigenvalues().real()[i]<<endl;
	}

//	cout<<"Autovlori e autovettori prima di essere riordinati con la libreria canonica"<<endl;
//	cout<<es.eigenvectors().real()<<endl<<endl;
//	cout<<es.eigenvalues().real()<<endl;
	//Ordiniamolo con l'algoritmo BubbleSort
	//Nota che per ogni autovalore riordinato bisogna anche riordinare la colonna che contiene l'autovettore
	double temp=0;
	int j=1;
	double Mat_eigenvec[N][N];
	double temp_mat;

	//Generiamo la matrice che ha per ogni colonna un autovettore
	for(int i=0; i<N; i++){
		for(int j=0; j<N; j++){
//		Mat_eigenvec[j][i]=es.eigenvectors().col(i)[j].real();
		Mat_eigenvec[j][i]=es.eigenvectors().col(i)[j];
		}
	}

	//Printiamo la matrice per vedere se coincide con gli autovettori stampati
//	cout<<"La matrice Mat_eigenvec"<<endl<<endl;
//	for(int i=0; i<N; i++){
//		for(int j=0; j<N; j++){
//		cout<<Mat_eigenvec[i][j]<<" ";
//		}
//		cout<<endl;
//	}

	while(j==1){
		j=0;
		for(int i=0; i<N-1; i++){
			if(eigenval[i]>eigenval[i+1]){
			//Scambia l'autovalore
			temp=eigenval[i];
			eigenval[i]=eigenval[i+1];
			eigenval[i+1]=temp;
			//Scambia le colonne della matrice (autovettore) quando viene scambiato il relativo autovettore
			for(int k=0; k<N-1; k++){
				temp_mat=Mat_eigenvec[k][i];
				Mat_eigenvec[k][i]=Mat_eigenvec[k][i+1];
				Mat_eigenvec[k][i+1]=temp_mat;
			}
			j=1;
			}
		}
	}
	//Facciamo un check per vedere se è ordinato
//	cout<<"Il vettore riordinato e la matrice riordinata sono:"<<endl<<endl;
//	for(int i=0; i<N; i++){
//		cout<<i<<"\t"<<eigenval[i]<<endl;
//	}

//	for(int i=0; i<N; i++){
//		for(int j=0; j<N; j++){
//		cout<<Mat_eigenvec[i][j]<<"  ";
//		}
//		cout<<endl;
//	}

//----->NORMALIZZA LE FUNZIONI D'ONDA
	//Normalizza la funzione d'onda 
	double S=0; //Costante di normalizzazione
	double v[N];

	int No;	//Numero di autovalor che si vuole stampare
	cout<<"Quante autofunzioni devo generare?"<<endl;
	cin>>No;

	string beg;
	beg="File_";
	string s;
	string tot;

for(int n=0; n<No; n++){
		//Parametri nome file
		s = n;
		tot=beg+s;
		//Generiamo un vettore contenente gli autovettori
		//Ogni colonna della matrice Mat_eigenvec contiene gli autovettori ordinati per autovalore in ordine crescente
		for(int i=0; i<N; i++){
			//v[i]=es.eigenvector().col(1)[i].real(); //In questo caso è la funzione d'onda ground state
			v[i]=Mat_eigenvec[i][n]; //In questo caso è il ground state
		}

		//Somma di Riemann
		for(int i=0; i<N; i++){
			S=S+h*v[i]*v[i];
		}
		cout<<"S: "<<S<<endl;
		//Dividiamo per la costante di normalizzazione
		for(int i=0; i<N; i++){
			v[i]=v[i]/sqrt(S);
		}

		//Check per la normalizzazione
		double newS=0;
		for(int i=0; i<N; i++){
			newS=newS+h*v[i]*v[i];
		}
		cout<<endl;
		cout<<"newS: "<<newS; //Se viene restituito 1 la normalizzazione è stata eseguita correttamente

		//.txt contenente la funzione d'onda normalizzata
		ofstream normgs(tot);
		normgs<<"#Autovalore"<<n<<endl<<endl;
		for(int i=0; i<N; i++){
			normgs<<rmin+i*h<<"\t"<<v[i]+eigenval[n]<<endl;
		}
		normgs.close();
	}
	return 0;
}

