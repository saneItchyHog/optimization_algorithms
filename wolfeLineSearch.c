#include "declarations.h"
#include "constants.h"

void wolfeLineSearch(double* current,double* direction,double* gradient,double* nextGradient, double* next){

  double sDotjacU0 = dot(direction,gradient,size);
  if(sDotjacU0 < 0){
    printf("wolfeLineSearch input direction is not descent\n");
    exit(0);
  }
  
  double fU0 = myFunction(current);
  double fullstep = 1;
  int i;
    
  double left = 0, mid = fullstep, right = fullstep;
  
  while(right - left > wolfeTol){
   
    for(i=0; i<size; i++){
      next[i] = current[i] - mid * direction[i];
    }
    
    myGradient(next,nextGradient);
    if(myFunction(next) <= fU0 - c1 * mid * sDotjacU0){
      if((fabs(dot(direction,nextGradient,size)) <= c2 * fabs(sDotjacU0))  || mid == fullstep){
        return;
      }
      else{
         left = mid;
         mid = (mid+right)/2;
      }
    }
    else{
      right = mid;
      mid = (mid+left)/2;
    }
  }
  
}
