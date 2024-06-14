#ifndef __MDETECTOR_H__
#define __MDETECTOR_H__
#include "MSystem.h"  

// return the ITS hit vector according to the hit bit mask

vector<int> GetITSHitVector(int hitmask){
  std::vector<int> hitvector;
  // calculate the max number of hits
  int maxhit=0;
  for(int i=0;i<7;i++){
    if(hitmask&(1<<i)) hitvector.push_back(i);
  }
  return hitvector;
}

#endif