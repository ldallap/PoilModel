
/*include essential libraries*/
#include<cstdio>
#include<cstdlib>
#include<cstring>
#include<cmath>
#include<iostream>
/**/

using namespace std;

/*NETWORK PARAMETERS*/
#define Lnetwork 10 /*Network size L: base of a rectangle*/
#define Hnetwork 5 /*Network height H: height of a rectangle*/
 
#define Lnetwork_height 3
#define Lnetwork_base 5

int Nneigh = int(Lnetwork_base*Lnetwork_height); /*Lnetwork_base * Lnetwork_height*/
int Nnetwork = int(Lnetwork*Hnetwork); /*Area: number of units*/

#define densExcSyn 0.5
#define densInhSyn 0.9
#define densExcNeurons 0.8
#define densInhNeurons 0.2



long seed = 123456;
int main(int argc, char *argv[]){


	int *ninneighbors;
	int **inneighbors; 
	int *nature; 
	int **outpotentialneighbors;

	int *relativeneighbors;
	relativeneighbors=new int[Nneigh];



	int Nneigh2=0;
	double ran2(long *idum);
	int xcoord(int i, int L);
	int ycoord(int i, int L);
	int icoord(int x, int y, int L); // icoord(x,y, Length side)


	// for(int i=0;i<Lnetwork_height;i++){
	// 	int deltax=int(-Lnetwork_base/2);
	// 	int deltay=int(-Lnetwork_height/2+i);
	// 	int indfirstinline=icoord(deltax,deltay,50);
	// 	for(int j=0;j<Lnetwork_base;j++){
	// 		relativeneighbors[i*Lnetwork_base+j]=indfirstinline+j;
	// 		// cout << i*Lnetwork_base+j << endl;
	// 	}
	// }


	// // for(int i=0;i<Nnetwork;i++){
	// 	for(int j=0;j<Nneigh;j++){
	// 		int i = 10;
	// 		int candidate = i+relativeneighbors[j];
	// 		cout << j << "\t" << candidate << endl;
	// 	}
	// // }

	// int auxidnei = 58;
	// int auxidneuron = 10;
	// cout << abs(xcoord(auxidneuron,Lnetwork)-xcoord(auxidnei,Lnetwork)) << "\t <= \t" << Lnetwork_base/2 << endl;
	// cout << abs(ycoord(auxidneuron,Lnetwork)-ycoord(auxidnei,Lnetwork)) << "\t <= \t" << Lnetwork_height/2 << endl;


	// return 1;

	long seed;
	nature = new int[Nnetwork];
	for(int k=0;k<Nnetwork;k++) nature[k]= (ran2(&seed) < densExcNeurons);
	/**/
	/**/
	outpotentialneighbors = new int*[Nnetwork];
	for(int k=0;k<Nnetwork;k++) outpotentialneighbors[k]=new int[Nneigh];
	/**/
	// int *relativeneighbors;
	// relativeneighbors=new int[Nneigh];

	// int *relativeneighbors2;
	// relativeneighbors2=new int[Nneigh];

	int *noutpotentialneighbors;
	noutpotentialneighbors=new int[Nnetwork];

	for(int i=0;i<Lnetwork_height;i++){
		int deltax=int(-Lnetwork_base/2);
		int deltay=int(-Lnetwork_height/2+i);
		int indfirstinline=icoord(deltax,deltay,Lnetwork);
		for(int j=0;j<Lnetwork_base;j++){
			relativeneighbors[i*Lnetwork_base+j]=indfirstinline+j;
			// cout << i*Lnetwork_base+j << endl;
		}
	}


	/*Calculates matrix of outgoing potential neighbors [Nnetwork][Nneighborhood]:*/
	for(int i=0;i<Nnetwork;i++){
		noutpotentialneighbors[i]=0;

		for(int j=0; j<Nneigh; j++){
			int candidate = 0;
			candidate = i+relativeneighbors[j];
			

			if((candidate >= 0) && (candidate<Nnetwork) && (candidate != i)){
				if((abs(xcoord(i,Lnetwork)-xcoord(candidate,Lnetwork))<=Lnetwork_base/2) && (abs(ycoord(i,Lnetwork)-ycoord(candidate,Lnetwork))<=Lnetwork_height/2)){
					outpotentialneighbors[i][noutpotentialneighbors[i]]=candidate;
					noutpotentialneighbors[i]++;
				}
			}
		}
	}


	int auxidd = 0;

	for(int i = 0; i < noutpotentialneighbors[auxidd]; i++)
		std::cout << outpotentialneighbors[auxidd][i] << std::endl;

	std::cout << "Test ended!" << std::endl;
	return 1;

	/**/
	ninneighbors = new int[Nnetwork];
	for(int k=0;k<Nnetwork;k++) ninneighbors[k]=0;
	/**/
	inneighbors = new int*[Nnetwork];
	for(int k=0;k<Nnetwork;k++) inneighbors[k]=new int[Nneigh];
	/**/
	/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
	/*%%%%%%%%%%%%%%%%-----Choose randomly among outgoing neighbors with UNIFORM probability densXSyn (X=Exc,Inh)-----%%%%%%%%%%%%%%%%%%%%%*/
	/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
	for(int i=0;i<Nnetwork;i++){
		for(int j=0;j<noutpotentialneighbors[i];j++){
			if(ran2(&seed)< nature[i]*densExcSyn + (1-nature[i])*densInhSyn){
				int selected_outgoing=outpotentialneighbors[i][j];
				inneighbors[selected_outgoing][ninneighbors[selected_outgoing]]=i;
				ninneighbors[selected_outgoing]++;
			}
		}
	}
	/**/



	// 325 exc
	// 323 inh
	int auxid = 220;

	for(int i=0;i<49;i++){
		if(outpotentialneighbors[auxid][i]!=0)std::cout << i << "\t" << outpotentialneighbors[auxid][i] << std::endl;
	}






 return 0;
}

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/*Organize the Network structure Coordinate functions:*/
int ycoord(int i, int L) { return (i / L); }
int xcoord(int i, int L) { return (i % L); }
int icoord(int x, int y, int L) { return L*y+x ;}
/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/*DO NOT EVEN COME HERE*/
/*You do not need to understand what is below - It comes from a book! Numerical C/C++ Recipes*/
/*Random Number:*/
#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)
double ran2(long *idum){
 int j;
 long k;
 static long idum2 = 123456789, iy = 0, iv[NTAB];
 double temp;
 if(*idum <= 0)
  {
   if(-(*idum) < 1)   *idum = 1;
   else   *idum = -(*idum);
   idum2 = (*idum);
    for(j = NTAB+7;j >= 0;j--)
     {
      k = (*idum)/IQ1;
      *idum = IA1*(*idum-k*IQ1)-k*IR1;
      if(*idum < 0)   *idum += IM1;
      if(j < NTAB)   iv[j] = *idum;
     }
   iy = iv[0];
  }

 k = (*idum)/IQ1;
 *idum = IA1*(*idum-k*IQ1)-k*IR1;
 if(*idum < 0)   *idum += IM1;

 k = idum2/IQ2;
 idum2 = IA2*(idum2-k*IQ2)-k*IR2;
 if(idum2 < 0)   idum2 += IM2;

 j = iy/NDIV;
 iy = iv[j]-idum2;
 iv[j] = *idum;
 if(iy < 1)   iy += IMM1;

 if((temp = AM*iy) > RNMX)   return RNMX;
 else   return temp;
}


















