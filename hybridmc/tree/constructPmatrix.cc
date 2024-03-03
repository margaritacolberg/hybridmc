// Copyright (c) 2018-2022 Jeremy Schofield
// SPDX-License-Identifier: BSD-3-Clause

#include <fstream>
#include <iostream>
#include <string>
#include <math.h>
#include <stdlib.h>
#include <algorithm>
#include <vector>
#include <iomanip>
#include <bitset>
#include <gsl/gsl_math.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>


//   Note, convention here is dot{c} = c*Kmatrix, so row sums of Kmatrix should be zero for prob. conservation

using Eigen::MatrixXd;

struct Configuration
{
    int index;
    double energy;
    double entropy;
    double freeEnergy;
    double prob;
    std::string binaryString;

};

bool sortByEnergy(const Configuration &lhs, const Configuration &rhs) { 
  if (lhs.energy == rhs.energy) return lhs.entropy > rhs.entropy;
  else return lhs.energy > rhs.energy; 
}

void removeRow(Eigen::MatrixXd& matrix, unsigned int rowToRemove)
{
    unsigned int numRows = matrix.rows()-1;
    unsigned int numCols = matrix.cols();

    std::cout << rowToRemove << " row to remove " << numRows << " num rows " << std::endl;
    if( rowToRemove < numRows )
        matrix.block(rowToRemove,0,numRows-rowToRemove,numCols) = matrix.block(rowToRemove+1,0,numRows-rowToRemove,numCols);

    matrix.conservativeResize(numRows,numCols);
}

void removeColumn(Eigen::MatrixXd& matrix, unsigned int colToRemove)
{
    unsigned int numRows = matrix.rows();
    unsigned int numCols = matrix.cols()-1;

    if( colToRemove < numCols )
        matrix.block(0,colToRemove,numRows,numCols-colToRemove) = matrix.block(0,colToRemove+1,numRows,numCols-colToRemove);

    matrix.conservativeResize(numRows,numCols);
}


int main(int argc, char** argv){
  
  if (argc < 4)
    {
      std::cerr << "Usage:  constructPmatrix <state file> <rate file> <beta value> [<fpt sink>]" << std::endl;
      exit(EXIT_FAILURE);
    }

  std::ifstream inFile(argv[1]);   
  if (inFile.bad()){
    std::cerr << "Unable to open file " << argv[1] << std::endl;
    exit(EXIT_FAILURE);
  }


  std::ifstream inFile2(argv[2]);
  if (inFile2.bad()){
    std::cerr << "Unable to open file " << argv[2] << std::endl;
    exit(EXIT_FAILURE);
  }


  int numStates = 0;
  std::string line;
  while (std::getline(inFile, line)) ++numStates;

  int sink = numStates-1;
  if (argc == 5) sink = atoi(argv[4]);

  double beta = atof(argv[3]);
  std::cerr << "  Beta value is " << beta << " for file " << argv[1] << " with " << numStates << " states:  sink for fpt is " << sink << std::endl;

  std::vector<Configuration> state(numStates);

  double *energy = new double[numStates];
  double *entropy = new double[numStates];

  inFile.clear();
  inFile.seekg(0, std::ios::beg);



  double Fmin = 100000000.0;
  double *freeEnergy = new double[numStates];
  double **kMatrix = new double*[numStates];
  double **pMatrix = new double*[numStates];

  double **tauMatrix = new double*[numStates];
  double kMax = 0.0;

  for (int i=0;i<numStates;i++){
    kMatrix[i] = new double[numStates];
    pMatrix[i] = new double[numStates];
    tauMatrix[i] = new double[numStates];
    for (int j=0;j<numStates;j++){
      kMatrix[i][j] = 0.0;
      pMatrix[i][j] = 0.0;
      tauMatrix[i][j] = 0.0;
    }

  }

  int minState = -1;

  for (int i=0; i < numStates;i++){
    int index;
    std::string stateBin;
    inFile >> index >> energy[i] >> entropy[i] >> stateBin;

    state[i].energy = energy[i];
    state[i].entropy = entropy[i];
    state[i].freeEnergy = 1000000.0;
    state[i].index = index;
    state[i].binaryString = stateBin;
 
    if (entropy[i] != -1.0) {
      freeEnergy[i] = beta*energy[i] - entropy[i];
      state[i].freeEnergy = freeEnergy[i];
    } else {
      freeEnergy[i] = 100000.0;
    }
    if (freeEnergy[i] < Fmin) {
      Fmin = freeEnergy[i];
      minState = i;
    }
    // std::cerr << "   state number " << index << "  e = " << energy[i] << "   S = " << entropy[i] << "   binary code is " << stateBin << std::endl;

  }

  while (!inFile2.eof()){
	int index1, index2, np;
        double dS, tauVal;
		//inFile2 >> index1 >> index2 >> tauVal >> dS >> np;
        inFile2 >> index1 >> index2 >> tauVal >> dS;
		//inFile2 >> index1 >> index2 >> tauVal >> np;
        tauMatrix[index1][index2] = tauVal;
  }


  double Z = 0.0;
  for (int i=0;i<numStates;i++){
    double deltaF = freeEnergy[i] - Fmin;
    if (deltaF < 20.0) Z += exp(-deltaF);
  }


  for (int i=0;i<numStates;i++){
    double prob_i = exp( Fmin - freeEnergy[i] )/Z;
    state[i].prob = prob_i;
  }

  //std::sort(state.begin(), state.end(), sortByEnergy);

  // record remapping in a file
  std::ofstream mappingFile("statemap.dat");
  for (int i=0;i<numStates;i++){
      mappingFile << i << " " << state[i].index << " " << state[i].energy 
                  << " " << state[i].entropy << std::endl;
  }

  std::ofstream probFile("Probability.dat");
  std::ofstream energyFile("Energy.dat");
  std::ofstream kfile("kProb.mtx");
  std::ofstream tsFile("ts.data");

  int elements = 0;
  double maxTS = -INFINITY;

  for (int i=0;i<numStates;i++) {
      probFile << state[i].prob << std::endl;
     
      // std::cerr << "   state number " << i << "  e = " << state[i].energy 
      //           << "   S = " << state[i].entropy << "  index = " << state[i].index << "   prob = " << state[i].prob
      //           << "   binary code is " << state[i].binaryString << std::endl;

      double ei = state[i].energy;
      if (state[i].entropy == -1) continue;
      energyFile << state[i].energy << std::endl;
      for (int j=0;j<i;j++){
		  double ej = state[j].energy;
		  //if (state[j].entropy == -1) continue;

		  //if ( (ej-ei) != 1.0) continue;

		  int ind1 = state[i].index;
		  int ind2 = state[j].index; // should be lower than ind1 since fewer bonds
		  std::bitset<26> state1(ind1);
		  std::bitset<26> state2(ind2);
		  std::bitset<26> connectedState = state1^state2;

		  //std::cerr << "  connected bit state for indices " << ind1 <<"-"<< ind2 << " is " << connectedState << std::endl;	
		  
		  if (connectedState.count() != 1) continue;
		  

		  double t_a = tauMatrix[i][j]; // index2 < index1, so tauM is average time to barrier from unbound state 
		  double t_b = tauMatrix[j][i]; // tauM[j][i] is average time to barrier from bound state (small) 

		  if (ind2 > ind1){
			  t_a = tauMatrix[j][i];
			  t_b = tauMatrix[i][j];
		  }

		  double i_kf = t_a + exp(state[j].freeEnergy - state[i].freeEnergy) * t_b;
		  double i_kr = exp(state[i].freeEnergy - state[j].freeEnergy) * t_a + t_b;

                  //std::cerr << "   t_a = " << t_a << "  and     tb=" << t_b << std::endl;

		  // state j has fewer bonds, so k_f for M_{i->j} since bond breaking
		  if (i_kf == 0.0){
			  std::cerr << " Error:  rate matrix " << i << "-" << j << "  bit states " << state1 << "-" << state2 << " has t_1 = " << t_a
						<< " and t_b = " << t_b << std::endl;
			  std::cerr << "  bonding integers are " << ind1 << " - " << ind2 << std::endl;
			  std::cerr << "  State " << i << ":  e = " << state[i].energy << "   S = " << state[i].entropy << "   bonds " << state[i].binaryString
						<< std::endl;
			  std::cerr << "  State " << j << ":  e = " << state[j].energy << "   S = " << state[j].entropy << "   bonds " << state[j].binaryString
						<< std::endl;
			  exit(1);
		  }
		  kMatrix[i][j] = 1.0/i_kf; //  note that j<i, so j state has few bonds, k_f is bond-breaking from state i to state j
		  kMatrix[j][i] = kMatrix[i][j]*state[i].prob/state[j].prob;
		  elements += 2;


		  if (state[i].energy < state[j].energy){
			  float zero = 0.0;
			  float ts = state[i].freeEnergy - log( kMatrix[i][j] )/beta;

			  if (maxTS < ts) maxTS = ts;
			  
			  tsFile << std::setw(10) << std::fixed << std::setprecision(3) << ts << " " 
					 << std::setw(10) << zero << " 1 " << i+1 << " "<< j+1;
			  tsFile << " " << std::setw(10) << zero << " " << std::setw(10) << zero << " " 
					 << std::setw(10) << zero << std::endl;

			  // now other direction
			  //tsFile << std::setw(10) << std::fixed << std::setprecision(3) << ts << " "
			  //       << std::setw(10) << zero << " 1 " << j+1 << " "<< i+1;
			  //tsFile << " " << std::setw(10) << zero << " " << std::setw(10) << zero << " "
			  //      << std::setw(10) << zero << std::endl;

		  }

		  // std::cerr << "  K[" << i << "][" << j << "] = " << kMatrix[i][j] << "  state index " << ind1 << "-" << ind2 << " t_a = " << t_a << "   t_b = " << t_b 
		  //           << "   pj/pi = " << state[j].prob/state[i].prob <<  "  Sj - Si = "
		  // 	  << state[j].entropy - state[i].entropy << "   k_f/k_r = " << i_kr/i_kf << std::endl; 	
      }
  }

 
  //

  //  compute the mean first-passage time from K matrix assuming an absorbing last state, starting from state 0
  std::cerr << " Built transition matrix array and ts.dat file." << std::endl;
 
  MatrixXd M(numStates, numStates);
  for (int i=0;i<numStates;i++)
      for (int j=0;j<numStates;j++) M(i,j) = kMatrix[j][i]; // transpose, so Matrix equation is dot{c} = M*c; // column sums zero for conserved prob.

  //
  //	remove sink elements
  std::cerr << "About to remove elements..." << std::endl;
  removeRow(M, sink);
  removeColumn(M, sink);

  std::cerr << "  Removed rows and columns." << std::endl;

 
  // check that columns sum to zero
  for (int i=0;i<numStates-1;i++){
	  double sum = 0.0;
	  for (int j=0;j<numStates-1;j++) sum += M(j,i);
	  M(i,i) = -sum;
          std::cerr << "  Sink matrix row " << i << " has Mii = " << sum << std::endl;
  }

  //  put in exit channels:  prob. not conserved in M dynamics
  for (int i=0;i<numStates-1;i++){
	  M(i,i) -= kMatrix[i][sink];// transition from i to final state
  }

  std::cout << "  sink matrix is:" << std::endl;
  int b = 10;
  //std::cout << M.block(numStates-1-b, numStates-1-b, b, b) << std::endl;
  
  double* data = new double[(numStates-1)*(numStates-1)];
  int index = 0;
  for (int i=0;i<numStates-1;i++){
    for (int j=0;j<numStates-1;j++) {
		data[index] = M(i, j); // pack array in format for gsl
		index++;
    }

  }



  gsl_matrix_view m = gsl_matrix_view_array (data, numStates-1, numStates-1);
  gsl_vector_complex *eval = gsl_vector_complex_alloc (numStates-1);
  gsl_matrix_complex *evec = gsl_matrix_complex_alloc (numStates-1, numStates-1);

  gsl_eigen_nonsymmv_workspace * w = gsl_eigen_nonsymmv_alloc (numStates-1);

  gsl_eigen_nonsymmv (&m.matrix, eval, evec, w);
  gsl_eigen_nonsymmv_free (w);
  gsl_eigen_nonsymmv_sort (eval, evec,
                           GSL_EIGEN_SORT_ABS_ASC);
   for (int i = 0; i < numStates-1; i++){
      gsl_complex eval_i = gsl_vector_complex_get (eval, i);
      //gsl_vector_complex_view evec_i = gsl_matrix_complex_column (evec, i);
      double eigenval_i = GSL_REAL(eval_i);
      std::cout << " Eigenvalue " << i << " is "  << eigenval_i << std::endl;
   }


  
  //gsl_matrix *Minv = gsl_matrix_calloc (numStates-1,numStates-1);
  //gsl_permutation *p = gsl_permutation_calloc(numStates-1);
  //int s;
  //gsl_linalg_LU_decomp(M, p, &s);
  //gsl_linalg_LU_invert(M, p, Minv);

  MatrixXd Minv = M.inverse();
  
  double tav = 0.0;
  for (int i=0;i<numStates-1;i++) {
	  tav -= Minv(i,0); // mfpt starting from state 0
  }
  std::cout << "#  mean first passage time = " << tav << "  for beta = " << beta << std::endl;

  //gsl_matrix_free (M);
  //gsl_matrix_free (Minv);

  

  double rowSumMax = 0.0;
  int maxRow = 0;
  for (int i=0;i<numStates;i++){
    double iMax = 0.0;
    for (int j=0;j<numStates;j++) iMax += kMatrix[i][j];
    kMatrix[i][i] = -iMax;
    // write out diagonal element of k matrix
    if (iMax > rowSumMax) {
      rowSumMax = iMax;
      maxRow = i;
    }
  }

  double h = 0.99999/(rowSumMax);
  std::cerr << "  rowSumMax = " << rowSumMax << " for row "<< maxRow << "   h=" << h << std::endl;

  std::ofstream pfile("tProb.mtx");
  std::ofstream pfile1("tProb.txt");
  pfile << numStates << " " << numStates << " " << elements + numStates << std::endl;
  kfile << numStates << " " << numStates << " " << elements + numStates << std::endl;
  for (int i=0;i<numStates;i++){
    double rowSum = 0.0;
    double kSum = 0.0;
    for (int j=0;j<numStates;j++){
      pMatrix[i][j] = h*kMatrix[i][j];
      if (i != j) rowSum += pMatrix[i][j];
      kSum += kMatrix[i][j];
      
    }
    //std::cerr << " sum of row " << i << " of kMatrix is " << kSum << std::endl;
    if (rowSum > 1.0){
      std::cerr << "Error:  h value is too large since rowSum = " << rowSum << ". Try reducing h." << std::endl;
      exit(EXIT_FAILURE);
    }
    pMatrix[i][i] = 1.0 - rowSum;
  }
  //  look at commitors to see if all states connected so have unique limit distribution
  MatrixXd T(numStates, numStates);
  for (int i=0;i<numStates;i++)
      for (int j=0;j<numStates;j++) {
         if (i==j) {
             T(i,j) = 1.0 - pMatrix[i][i];
         } else {
             T(i,j) = -pMatrix[i][j]; // row sums zero for conserved prob.
         }
      }

  removeRow(T, sink);
  removeColumn(T, sink);
  removeRow(T, 0); // source is first state
  removeColumn(T,0);

  for (int i=0;i<numStates-2;i++){
      double rowSum = 0.0;
      for (int j=0;j<numStates-2;j++) rowSum += T(i,j);
      if (rowSum != 1.0) std::cerr << " T matrix row sum for row " << i << " is " << rowSum << std::endl;
  }

  MatrixXd Tinv = T.inverse();
  for (int i=0;i<numStates-2;i++){
     double committor = 0.0; 
     for (int j=0;j<numStates-2;j++) committor += Tinv(i,j)*pMatrix[j+1][sink];
      if (committor == 0.0) std::cout << "  State " << i << " is disconnected from final state. " << std::endl;

  }
  



  // now output matrix
  for (int i=0;i<numStates;i++){
    double testSum = 0.0;
    double testSum2 = 0.0;
    for (int j=0; j < numStates;j++){
      if (pMatrix[i][j] != 0.0) {
	pfile << i+1 << " " << j+1 << " " << std::setprecision(16) << pMatrix[i][j] << std::endl;
	pfile1 << i << " " << j << " " << std::setprecision(16) << pMatrix[i][j] << std::endl;
      }

      if (kMatrix[i][j] != 0.0) {
		  kfile << i << " " << j << " " << std::setprecision(16) << kMatrix[i][j] << std::endl;
		  testSum2 += kMatrix[i][j];
      }
      testSum += pMatrix[i][j];
  
    }
    //std::cerr << "   Row " << i << " sums to " << testSum << std::endl;
     //std::cerr << "  Now row " << i << " sums to " << testSum2 << std::endl;


  }

  //
  // for (int i=0;i<numStates;i++)
  //   for (int j=0;j<numStates;j++){
  //     if (i==j or pMatrix[j][i] == 0.0) continue;
  //     double testBalance = pMatrix[i][j]/pMatrix[j][i];
  //     double ratio = state[j].prob/state[i].prob;
  //     std::cerr << "  Detailed balance test for T[" << i << "][" << j <<  "] " << state[i].prob * pMatrix[i][j] << " = " << state[j].prob * pMatrix[j][i] << std::endl;
  //   }


  //  now output info for disconnectivity
  int levels = 100;
  //double delta = 10.0/double(levels);
  double delta = 50.0/double(levels);
  std::ofstream infoFile("dinfo");
  double first = Fmin + levels*delta;
  if (maxTS < first) {
      first = maxTS;
      delta = (maxTS - Fmin)/double(levels);
  }
  //infoFile << "first " << std::fixed << std::setw(10) << std::setprecision(3) << maxTS << std::endl;
  infoFile << "first " << std::fixed << std::setw(10) << std::setprecision(3) 
          << first << std::endl;
  infoFile << "levels " << levels << std::endl;
  infoFile << "delta " << std::fixed << std::setw(10) << std::setprecision(3) << delta << std::endl;
  infoFile << "minima min.data" << std::endl;
  infoFile << "ts ts.data" << std::endl;
  //infoFile << "connectmin 1" << std::endl;
  //infoFile << "idmin " << 1024 << std::endl;
  infoFile << "idmin " << 1 << std::endl;
  infoFile << "idmin " << 512 << std::endl;
  //infoFile << "idmin " << 224 << std::endl;
  //infoFile << "idmin " << 892 << std::endl;
  //infoFile << "idmin " << 60 << std::endl;
  //infoFile << "identify" << std::endl;
  infoFile << "trval " << numStates << " bond.dat 1.0 false 0.0"  << std::endl;
  infoFile << "trvalrange 2 " << -int(state[numStates-1].energy) << std::endl;
  infoFile << "trvalscale"  << std::endl;
infoFile << "tsthresh "  << Fmin + levels*delta << std::endl;
//  infoFile << "tsthresh "  << 0 << std::endl;

  std::ofstream minFile("min.data");
  for (int i=0;i<numStates;i++){
    float zero = 0.0;
    minFile << std::fixed << std::setw(10) << std::setprecision(3) << state[i].freeEnergy << " ";
    minFile << std::fixed << std::setw(10) << std::setprecision(3) << zero << " 1 ";
    minFile << std::fixed << std::setw(10) << std::setprecision(3) << zero << " ";
    minFile << std::fixed << std::setw(10) << std::setprecision(3) << zero << " ";
    minFile << std::fixed << std::setw(10) << std::setprecision(3) << zero << std::endl;
  }

  std::ofstream bondFile("bond.dat");
  for (int i=0;i<numStates;i++) bondFile << -int(state[i].energy) << std::endl;

 

  

  delete[] freeEnergy;
  delete[] energy;
  delete[] entropy;

  return 0;
}

