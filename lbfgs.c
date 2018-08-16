#include "declarations.h"
#include "constants.h"

double* lbfgs(double *initial){


  // Copy the initial input, x_0, if needed by user later, to current x_i
  
  double* current = (double*)malloc(size * sizeof(double));
  int i,j;
  
  for(i=0; i<size; i++){
    current[i] = initial[i];
	  //printf("%12.8f\n",current[i]);
  }
  // next solution, x_i+1
  double* next = (double*)malloc(size * sizeof(double));
  
  // allocate the memory for the history, required by lbfgs
  // array of the differences of the past gradients, yks[i] = \del(f)(x_{i+1}) - \del(f)(x_i)
  double** ykS = (double**)malloc(historySize * sizeof(double*));
  for(i=0; i<historySize; i++){
    ykS[i] = (double*)malloc(size * sizeof(double));
  }
  
  //array of the differences of the past solution, sigmaS[i] = x_{i+1} - x_{i}
  double** sigmaS = (double**)malloc(historySize * sizeof(double*));
  for(i=0; i<historySize; i++){
    sigmaS[i] = (double*)malloc(size * sizeof(double));
  }
  
  // array of inverse of product of differences of gradient and solutions
  // rowS[i] = \frac{1}{\sigmaS[i] \cdot ykS[i]}
  double* rowS = (double*)malloc(historySize * sizeof(double));
  
  // gradient at current soultion
  double* gradient = (double*)malloc(size * sizeof(double));
  double* nextGradient = (double*)malloc(size * sizeof(double));
  double* q = (double*)malloc(size * sizeof(double));
  
  // to store local variables aS, see lfbgs
  double* aS = (double*)malloc(historySize * sizeof(double));
  
    
  int iterationNo = 0;
  myGradient(current,gradient);
  
  double fk = myFunction(current);
  double fk1,dxNorm,gradNorm,dfun;

  do{
    double b;
    // q = gradient;
    for(i=0; i<size; i++){
      q[i] = gradient[i];
    }
    
    for(i=0; i<historySize; i++){
      int k = (-i+iterationNo+historySize) % historySize;
      
      // a_i = \row_i \cdot sigma_i^T q
      aS[k] = rowS[k] * dot(sigmaS[k],q,size);
      
      // q = q - a_i \cdot y_i
      for(j=0; j<size; j++){
        q[j] = q[j] - aS[k] * ykS[k][j]; 
      }
    }
    // z = q
    for(i=0; i<historySize; i++){
      int k = (iterationNo+i+1) % historySize;
      // b = \row_i \cdot yk_i^T z
      b = rowS[k] * dot(ykS[k],q,size);
      
      // z = z + (a_i-b) sigma_i 
      for(j=0; j<size; j++){
        q[j] = q[j] + (aS[k]-b) * sigmaS[k][j]; 
      }
    }
    
    wolfeLineSearch(current,q,gradient,nextGradient,next);
    double fk1 = myFunction(next);
    int l = (iterationNo+1) % historySize;
    for(j=0; j<size; j++){
      sigmaS[l][j] = next[j] - current[j];
      ykS[l][j] = nextGradient[j] - gradient[j];
    }
    rowS[l] = 1/dot(sigmaS[l],ykS[l],size);

    
    dxNorm = sqrt(dot(sigmaS[l],sigmaS[l],size));
    gradNorm = sqrt(dot(gradient,gradient,size));
    dfun = fabs(fk-fk1)/(1+fabs(fk1));

    switchPointer(next,current);
    switchPointer(nextGradient,gradient);
    iterationNo++;
    
    
  }while((dfun > epsfun) && (dxNorm > epsx) && (gradNorm > epsgrad) && (iterationNo < maxIter));

  printf("solution at iteration %d\n",iterationNo);
  for(i=0; i<size; i++){
	printf("%.30lf ",next[i]);
  }
  printf("\n");    
  return next;

}