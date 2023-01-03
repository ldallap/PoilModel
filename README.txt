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