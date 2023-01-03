/*
  CROS model described in:

	Leonardo Dalla Porta and Mauro Copelli (2018).
  Modeling neuronal avalanches and long-range temporal
	correlations at the emergence of collective oscillations:
	continuously varying exponents mimic M/EEG results.
  DOI: https://doi.org/10.1101/423921

  This reference should be cited whenever this code is used.

  COMPILE: g++ PoilModel_byLeonardoDPD_versionUpdate11Out2018.cpp
	-O3 -march=native -o poil

  RUN: ./poil SideSizeNetwork SideSizeNeighbor densExcSyn densInhSyn TimeSimulation SeedRandomNumber
  EX:  ./poil 50 7 0.3 0.4 10000 94520934

  Writen by Leonardo Dalla Porta and Mauro Copelli 2016   (updated 2018)
	Contact: leonardodallaporta@gmail.com
*/

/**************************************************************************/


/*include essential libraries*/
#include<cstdio>
#include<cstdlib>
#include<cstring>
#include<cmath>
#include<iostream>
/**/



int main(int argc, char *argv[]){

	FILE *outf;
	char Outfile[500];

	int Lnetwork; /*Lateral size of Network square*/
	int Lneighborhood; /*Lateral size of Neighboors square*/
	int tmax; /*Total number of steps - simulations time*/
	long seed; /*Random number seed*/
	double densExcSyn; /*Parameter rE - see Reference [1]*/
	double densInhSyn; /*Parameter rI - see Reference [1]*/

	int Nnetwork; /*Total size square network - LxL - see Reference [1]*/
	int Nneighborhood; /*Total size square neighborhood - lxl - see Reference [1]*/

	/*Variables to be used on the network construction*/
  int xcoord(int i, int L);
	int ycoord(int i, int L);
	int icoord(int x, int y, int L);
	int candidate,selected_outgoing;
	int deltax,deltay,indfirstinline;
	/**/

	/*Variables to set the neighborhood*/
	int **inneighbors,**outpotentialneighbors;
	int *relativeneighbors;
	int *noutpotentialneighbors;
	int *ninneighbors; /*Will store the Neighboors*/
	/**/

	int *nature;/*Neuron Nature 1-Exc 0-Inh*/
	int *spike;/*0=non-spike 1=spike*/

	double *current; /*Neuron Current - Variable I - see Reference [1]*/
	double *pspike; /*Neuron Prob Spike - Variable P - see Reference [1]*/

	int nspikes; /*Network activity A(t) per time step - See Reference [1]*/

	double ran2(long *idum); /*Call the randomnumber - describe on the page's end*/

	/*Functions to check the Network structure*/
	void printsites(int Lnetwork);
	void printsites_displaced(int Lnetwork);
	/**/


	/*&Constants to be used - Extract exaclty from Reference [2] - see also Reference [1]&*/
	/*Conditions for equations: see Reference [1] or [2]*/
	const double tauI = 0.009;
	double tauP [2];
	tauP [1] = 0.006;/*Constant TauP(excitatory) = 6ms = 0.006s*/
	tauP [0] = 0.012;/*Constant TauP(inhibitory) = 12ms = 0.012s*/
	const double dt = 0.001; /*time step*/
	const double Izero = 0.0;/*Izero = Initial Current*/
	double P[2];
	P[0] = 0.0;/*Poisson Probability for Inh-Cells*/
	P[1] = 0.000001;/*Poisson Probability for Exc-Cells*/
	double W[2][2];/*Synpase Weight - see Reference [1] or [2]*/
	W[0][0] = -2.0; /*Inh-Inh*/
	W[0][1] = 0.011;/*Inh-Exc*/
	W[1][1] = 0.02;/*Exc-Exc*/
	W[1][0] = -2.0;/*Exc-Inh*/
	const double dtovertauI= dt/tauI; /*See Reference [1]*/
	const double oneminusdtovertauI= 1.- dt/tauI;/*See Reference [1]*/
	double dtovertauP [2];/*See Reference [1]*/
	dtovertauP [1] = dt/tauP[1];/*See Reference [1]*/
	dtovertauP [0] = dt/tauP[0];/*See Reference [1]*/
	double oneminusdtovertauP [2];/*See Reference [1]*/
	oneminusdtovertauP [1] = 1.-dt/tauP[1];/*See Reference [1]*/
	oneminusdtovertauP [0] = 1.-dt/tauP[0];/*See Reference [1]*/
	const double Presetexc=-2.;/*Probability Reset for Inh - See Reference [1]*/
	const double Presetinh=-20;/*Probability Reset for Exc - See Reference [1]*/
	const double initialdensity = 0.01;/*# of neurons that start firing*/
	const double densExcNeurons=0.75; /*% of exc neurons*/
	const double densInhNeurons=1.-densExcNeurons;/*% of inh neurons*/


	/*************************************************/
	/*Extract the INPUT information*/
	/* Check command-line arguments */
  if(argc!=7){
		printf("ERROR\n");
		printf("Usage: poilmodel Lnetwork Lneighborhood densExcSyn densInhSyn tmax seed \n");
		printf("Ex: ./executable 50 7 0.3 0.4 10000 94520934\n");
    return 1;
  }
  else{
    Lnetwork = atoi(argv[1]);
    Lneighborhood = atoi(argv[2]);
    densExcSyn = (double) atof(argv[3]);
    densInhSyn = (double) atof(argv[4]);
    tmax = atoi(argv[5]);
    seed = (long) atoi(argv[6]);
  }
	/*************************************************/
	/*************************************************/
	for(int z=0;z<10000;z++) double nothing=ran2(&seed); /*warm up the randomnumber*/
	/*************************************************/
	Nnetwork = Lnetwork*Lnetwork;
	Nneighborhood = Lneighborhood*Lneighborhood;
	/*Name of the outfile containing the # of Spikes per dt*/
	sprintf(Outfile,"Lnetw%i_Lneig%i_dExcS%0.2f_dInhS%0.2f_seed%ld_tmax%g.ts",Lnetwork,Lneighborhood,densExcSyn,densInhSyn,seed,(double)tmax);
	/*************************************************/

	/*Check any inconsistency*/
	/*************************************************/
	/*Lneighborhood must be uneven (odd):*/
	if((Lneighborhood%2)==0){
		printf("Lneighborhood cannot be even!\n");
		return 1;
	}
  /*Lnetwork must be greater than or equal (>) Lneighborhood:*/
	if(Lnetwork < Lneighborhood){
		printf("Lnetwork must be larger than Lneighborhood!\n");
		return 1;
	}
	/*************************************************/

	/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
	/*Allocate memory and check for possible failure*/
	// Creates vector for relative position of neighbors [Nneighborhood]:
	relativeneighbors=new int[Nneighborhood];
	if(relativeneighbors == NULL){
		printf("ERROR: Can't allocate memory for the pointer relativeneighbors!");
		return 1;
	}
	/**/
	/*Allocate memory and check for possible failure*/
	// Creates matrix of outgoing potential neighbors [Nnetwork][Nneighborhood] :
	outpotentialneighbors = new int*[Nnetwork];
	if(outpotentialneighbors==NULL){
		printf("ERROR: Can't allocate memory for outpotentialneighbors");
		return 1;
	}
	for(int k=0;k<Nnetwork;k++){
		outpotentialneighbors[k]=new int[Nneighborhood];
		if(outpotentialneighbors[k]==NULL){
			printf("ERROR: Can't allocate 2nd position memory for the pointer outpotentialneighbors");
			return 1;
		}
	}
	/******************************/
	/*Allocate memory and check for possible failure*/
	// Creates vector of relative indexes for neighbors [Nneighborhood]:
	for(int i=0;i<Lneighborhood;i++){
		deltax=int(-Lneighborhood/2);
		deltay=int(-Lneighborhood/2+i);
		indfirstinline=icoord(deltax,deltay,Lnetwork);
		for(int j=0;j<Lneighborhood;j++){
			relativeneighbors[i*Lneighborhood+j]=indfirstinline+j;
			/*Uncomment next line to check any possible inconsistency on the network*/
			//printf("%d \t %d\n",i*Lneighborhood+j,indfirstinline+j);
		}
  }
	/******************************/
	/******************************/
	// Creates vector for number of outgoing potential neighbors [Nnetwork]:
	noutpotentialneighbors=new int[Nnetwork];
	if(noutpotentialneighbors==NULL){
		printf("ERROR: Can't allocate memory for noutpotentialneighbors");
		return 1;
	}
	/******************************/
	/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
	/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
	/*Calculates matrix of outgoing potential neighbors [Nnetwork][Nneighborhood]:*/
	for(int i=0;i<Nnetwork;i++){
		noutpotentialneighbors[i]=0;
		for(int j=0;j<Nneighborhood;j++){
			candidate= i+relativeneighbors[j];
			if((candidate >= 0) && (candidate<Nnetwork) && (candidate != i)){
				if((abs(xcoord(i,Lnetwork)-xcoord(candidate,Lnetwork))<=Lneighborhood/2) && (abs(ycoord(i,Lnetwork)-ycoord(candidate,Lnetwork))<=Lneighborhood/2)){
					outpotentialneighbors[i][noutpotentialneighbors[i]]=candidate;
					noutpotentialneighbors[i]++;
				}
			}
		}
	}
	/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%--------- MATRIX OF INCOMING NEIGHBORS: -----------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
	// Creates vector for number of incoming neighbors [Nnetwork]:
	ninneighbors = new int[Nnetwork];
	if(outpotentialneighbors==NULL){
		printf("Error: ERROR: Can't allocate memory for outpotentialneighbors.");
		return 1;
	}
/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
	// Creates matrix of incoming neighbors [Nnetwork][Nneighborhood] :
	inneighbors = new int*[Nnetwork];
	if(inneighbors==NULL){
		printf("ERROR: Can't allocate memory for outpotentialneighbors");
		return 1;
	}
	for(int k=0;k<Nnetwork;k++){
		inneighbors[k]=new int[Nneighborhood];
		if(inneighbors[k]==NULL){
			printf("ERROR: Can't allocate 2nd position memory for the pointer outpotentialneighbors");
			return 1;
		}
	}
/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
	// Creates vector for nature (excitatory or inhibitory) of the neuron [Nnetwork]:
	nature = new int[Nnetwork];
	if(nature==NULL){
		printf("ERROR: Can't allocate memory for outpotentialneighbors");
		return 1;
	}
	// Nature of neuros: +1 for excitatory neurons and 0 for inhibitory neurons:
	for(int k=0;k<Nnetwork;k++) nature[k]= (ran2(&seed) < densExcNeurons);
/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
	// Starting ninneighbors from 0:
	for(int k=0;k<Nnetwork;k++) ninneighbors[k]=0;
/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/*%%%%%%%%%%%%%%%%-----Choose randomly among outgoing neighbors with UNIFORM probability densXSyn (X=Exc,Inh)-----%%%%%%%%%%%%%%%%%%%%%*/
/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
	// for(int i=0;i<Nnetwork;i++){
	// 	for(int j=0;j<noutpotentialneighbors[i];j++){
	// 		if(ran2(&seed)< nature[i]*densExcSyn + (1-nature[i])*densInhSyn){
	// 			selected_outgoing=outpotentialneighbors[i][j];
	// 			inneighbors[selected_outgoing][ninneighbors[selected_outgoing]]=i;
	// 			ninneighbors[selected_outgoing]++;
	// 		}
	// 	}
	// }
/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/*%%%%%%%%%%%-----Choose randomly among outgoing neighbors with EXPONENTIALLY DECAYING probability densXSyn (X=Exc,Inh)-----%%%%%%%%%%%*/
/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
	for(int i=0;i<Nnetwork;i++){
		double sumnorm=0.;
		double Cnorm=0.;
		for(int j=0;j<noutpotentialneighbors[i];j++){
			sumnorm+=exp(-sqrt(pow(xcoord(i,Lnetwork)-xcoord(outpotentialneighbors[i][j],Lnetwork),2)+pow(ycoord(i,Lnetwork)-ycoord(outpotentialneighbors[i][j],Lnetwork),2)));
		}
		Cnorm = noutpotentialneighbors[i]*(nature[i]*densExcSyn + (1-nature[i])*densInhSyn)/sumnorm;
		for(int j=0;j<noutpotentialneighbors[i];j++){
			if(ran2(&seed)<(Cnorm*exp(-sqrt(pow(xcoord(i,Lnetwork)-xcoord(outpotentialneighbors[i][j],Lnetwork),2)+\
				pow(ycoord(i,Lnetwork)-ycoord(outpotentialneighbors[i][j],Lnetwork),2))) ) ){
				selected_outgoing=outpotentialneighbors[i][j];
				inneighbors[selected_outgoing][ninneighbors[selected_outgoing]]=i;
				ninneighbors[selected_outgoing]++;
			}
		}
	}
/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
	// Print sites and incoming neighbors:

	// printsites(Lnetwork);
	// for(int i=0;i<Nnetwork;i++){
	// 	printf("Neuron %d has %d incoming neighbors:\t",i,ninneighbors[i]);
	// 	for(int j=0;j<ninneighbors[i];j++){
	// 		if(nature[inneighbors[i][j]]){
	// 			printf("+");
	// 		}
	// 		else printf("-");
	// 		printf("%d\t",inneighbors[i][j]);
	// 	}
	// 	printf("\n");
	// }

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%------Creates the system of equations------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
	// Creates vector CURRENT:
	current = new double[Nnetwork];
	if(current==NULL){
		printf("ERROR: Can't allocate memory for outpotentialneighbors");
		return 1;
	}
/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  // Creates vector PSPIKE (Probability of Spike):
	pspike = new double[Nnetwork];
	if(pspike==NULL){
		printf("ERROR: Can't allocate memory for outpotentialneighbors");
		return 1;
	}
/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
  // Creates vector SPIKE (Sj):
	spike = new int[Nnetwork];
	if(spike==NULL){
		printf("ERROR: Can't allocate memory for outpotentialneighbors");
		return 1;
	}
/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
	/*Initializate the network with INITIALDENSITY neurons firing*/
	nspikes = 0;
	for(int i=0; i<Nnetwork; i++){
		current[i]=0.;			//Starting vector CURRENT with random number between 0 and 1
		pspike[i]=0.;				//Starting vector PSPIKE with random number between 0 and 1
		if(ran2(&seed)<initialdensity) spike[i]=1;
		else spike[i]=0;
	}
	/**/


	/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
	/*BLOCK FOR ITERATIONS OVER TIME*/
	/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
	/*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/
	/*HERE STARTS THE LOOPING OVER TIME*/
	outf = fopen(Outfile,"w");
	for(int t=0; t<tmax; t++){
		/*Firts Update the current of all neurons based on pre-synpatic spikes*/
		for(int i=0;i<Nnetwork; i++){
			for(int j=0;j<ninneighbors[i];j++){
				current[i] = current[i] + W[nature[i]][nature[inneighbors[i][j]]] * spike[inneighbors[i][j]];
			}
		}
		/**/
		/*reset vector spikes*/
	  memset((void*)spike, 0, Nnetwork*sizeof(int));
		/**/

	  nspikes= 0;
	  for(int i=0; i<Nnetwork; i++){
			current[i] = current[i] * oneminusdtovertauI; /*Equation for current*/
			pspike[i] = pspike[i] + current[i]; /*Equation for Prob of Spike*/
			pspike[i] = pspike[i] * oneminusdtovertauP[nature[i]] + P[nature[i]] * dtovertauP[nature[i]]; /*Eqaution for Pro of Spike*/

			if(ran2(&seed)<pspike[i]){ /*Theta Heavside*/
				if(nature[i]==1) pspike[i]= -2.0; /*for Exc-cells*/
				else pspike[i]= -20.0; /*for Inh-cells*/

				spike[i]=1;
				nspikes++;
			}
		}
		/*Print the # of spikes for each step of time - time is discrete*/
		fprintf(outf,"%i\n", nspikes);
	}
	fclose(outf);
	/*END LOOPING OVER TIME*/
	/*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/
	/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
	/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

	/*CHECK THE NETWORK FUNCTION*/
	// Print sites and potential outgoing neighbors:
	// printsites(Lnetwork);
	// for(int i=0;i<Nnetwork;i++){
	// 	printf("Neuron %d has %d potential neighbors:\t",i,noutpotentialneighbors[i]);
	// 	for(int j=0;j<noutpotentialneighbors[i];j++){
	// 		printf("\t%d",outpotentialneighbors[i][j]);
	// 	}
	// 	printf("\n");
	// }
	/***************************/


	/*De-allocate memory and destruct*/
	for(int i=0; i<Nnetwork; i++){
		delete [] inneighbors[i];
		delete [] outpotentialneighbors[i];
	}
	delete [] inneighbors;
	delete [] outpotentialneighbors;
	delete [] noutpotentialneighbors;
	delete [] ninneighbors;
  delete [] relativeneighbors;
  delete [] current;
  delete [] pspike;
  delete [] spike;
	delete [] nature;
}
/*end of main*/
/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%------ END OF MAIN ------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/*Organize the Network structure Coordinate functions:*/
int ycoord(int i, int L) { return (i / L); }
int xcoord(int i, int L) { return (i % L); }
int icoord(int x, int y, int L) { return L*y+x ;}
/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/*Function to check Network structure*/
/*Print sites:*/
// void printsites(int Lnetwork){
// 	for(int i=0;i<Lnetwork;i++){
// 		for(int j=0;j<Lnetwork;j++)	printf("%d\t",i*Lnetwork+j);
// 	}
// 	printf("\n");
// }
/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/*Function to check Network structure*/
/*Print sites:*/
// void printsites_displaced(int Lnetwork){
// 	for(int i=0;i<Lnetwork;i++){
// 		for(int j=0;j<Lnetwork;j++)
// 			printf("%d\t",i*Lnetwork+j - Lnetwork*Lnetwork/2);
// 	}
// 	printf("\n");
// }
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
/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
