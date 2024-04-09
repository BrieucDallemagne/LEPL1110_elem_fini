/*
 *  main.c
 *  Projet 2023-2024
 *  Elasticite lineaire plane
 *
 *  Code de calcul
 *
 *  Copyright (C) 2023 UCL-IMMC : Vincent Legat, Miguel De Le Court
 *  All rights reserved.
 *
 */

#include "fem.h"
#include <time.h>
int main(void) {

  
  femGeo *theGeometry = geoGetGeometry();
  geoMeshRead("../data/mesh.txt");
  femProblem *theProblem = femElasticityRead(theGeometry, "../data/problem.txt");
  femElasticityPrint(theProblem);
  clock_t begin = clock();

  double *theSoluce = femElasticitySolve(theProblem);

  int nNodes = theGeometry->theNodes->nNodes;
  femSolutionWrite(nNodes, 2, theSoluce, "../data/UV.txt");
  femSolutionWrite(nNodes, 2, theSoluce, "../../ProjectPostProcessor/data/UV.txt");
  
  clock_t end = clock();
  double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
  printf("Time spent: %f seconds\n", time_spent);

  femElasticityFree(theProblem);
  geoFree();
  return 0;
}
