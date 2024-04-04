#include "fem.h"

// Il faut un fifrelin generaliser ce code.....
//  (1) Ajouter l'axisymétrique !    (mandatory) // Ok (à vérifier)
//  (2) Ajouter les conditions de Neumann !   (mandatory) // A faire
//  (3) Ajouter les conditions en normal et tangentiel !   (strongly advised) // A faire (à vérifier)
//  (4) Et remplacer le solveur plein par un truc un fifrelin plus subtil  (mandatory) // OK (à vérifier)

void femElasticityAssembleElements(femProblem *theProblem) {
  femFullSystem *theSystem = theProblem->system;
  femIntegration *theRule = theProblem->rule;
  femDiscrete *theSpace = theProblem->space;
  femGeo *theGeometry = theProblem->geometry;
  femNodes *theNodes = theGeometry->theNodes;
  femMesh *theMesh = theGeometry->theElements;
  double x[4], y[4], phi[4], dphidxsi[4], dphideta[4], dphidx[4], dphidy[4];
  int iElem, iInteg, iEdge, i, j, d, map[4], mapX[4], mapY[4];
  int nLocal = theMesh->nLocalNode;
  double a = theProblem->A;
  double b = theProblem->B;
  double c = theProblem->C;
  double rho = theProblem->rho;
  double gx = theProblem->gx;
  double gy = theProblem->gy;
  double **A = theSystem->A;
  double *B = theSystem->B;

  for (iElem = 0; iElem < theMesh->nElem; iElem++) {
    for (j = 0; j < nLocal; j++) {
      map[j] = theMesh->elem[iElem * nLocal + j];
      mapX[j] = 2 * map[j];
      mapY[j] = 2 * map[j] + 1;
      x[j] = theNodes->X[map[j]];
      y[j] = theNodes->Y[map[j]];
    }

    for (iInteg = 0; iInteg < theRule->n; iInteg++) {
      double xsi = theRule->xsi[iInteg];
      double eta = theRule->eta[iInteg];
      double weight = theRule->weight[iInteg];
      femDiscretePhi2(theSpace, xsi, eta, phi);
      femDiscreteDphi2(theSpace, xsi, eta, dphidxsi, dphideta);

      double dxdxsi = 0.0;
      double dxdeta = 0.0;
      double dydxsi = 0.0;
      double dydeta = 0.0;
      for (i = 0; i < theSpace->n; i++) {
        dxdxsi += x[i] * dphidxsi[i];
        dxdeta += x[i] * dphideta[i];
        dydxsi += y[i] * dphidxsi[i];
        dydeta += y[i] * dphideta[i];
      }
      double jac = dxdxsi * dydeta - dxdeta * dydxsi;
      if (jac < 0.0) printf("Negative jacobian! Your mesh is oriented in reverse. The normals will be wrong\n");
      jac = fabs(jac);

      for (i = 0; i < theSpace->n; i++) {
        dphidx[i] = (dphidxsi[i] * dydeta - dphideta[i] * dydxsi) / jac;
        dphidy[i] = (dphideta[i] * dxdxsi - dphidxsi[i] * dxdeta) / jac;
      }
      for (i = 0; i < theSpace->n; i++) {
        for (j = 0; j < theSpace->n; j++) {
          if (theProblem->planarStrainStress = AXISYM || PLANAR_STRESS){
             A[mapX[i]][mapX[j]] += (dphidx[i] * a * dphidx[j] * x[i] + dphidy[i] * c * dphidy[j] * x[i] + dphidx[i] * b * phi[j] + phi[i] * (b * dphidx[j] + a * phi[j]/x[i])) * jac * weight;
             A[mapX[i]][mapY[j]] += (dphidx[i] * b * dphidy[j] * x[i] + dphidy[i] * c * dphidx[j] * x[i] + phi[i] * b * dphidy[j]) * jac * weight;
             A[mapY[i]][mapX[j]] += (dphidy[i] * b * dphidx[j] * x[i] + dphidx[i] * c * dphidy[j] * x[i] + dphidy[i] * b * phi[j]) * jac * weight;
             A[mapY[i]][mapY[j]] += (dphidy[i] * a * dphidy[j] * x[i] + dphidx[i] * c * dphidx[j] * x[i]) * jac * weight;
          }else{
            A[mapX[i]][mapX[j]] += (dphidx[i] * a * dphidx[j] + dphidy[i] * c * dphidy[j]) * jac * weight;
            A[mapX[i]][mapY[j]] += (dphidx[i] * b * dphidy[j] + dphidy[i] * c * dphidx[j]) * jac * weight;
            A[mapY[i]][mapX[j]] += (dphidy[i] * b * dphidx[j] + dphidx[i] * c * dphidy[j]) * jac * weight;
            A[mapY[i]][mapY[j]] += (dphidy[i] * a * dphidy[j] + dphidx[i] * c * dphidx[j]) * jac * weight;
          }
          
        }
      }
      for (i = 0; i < theSpace->n; i++) {
        if (theProblem->planarStrainStress = AXISYM || PLANAR_STRESS){
          B[mapX[i]] += phi[i] * gx * rho * x[i] * jac * weight;
          B[mapY[i]] += phi[i] * gy * rho * x[i] * jac * weight;
        }else{
          B[mapX[i]] += phi[i] * gx * rho * jac * weight;
          B[mapY[i]] += phi[i] * gy * rho * jac * weight;
        }
      }
    }
  }
}

void femElasticityAssembleNeumann(femProblem *theProblem) {
  femFullSystem *theSystem = theProblem->system;
  femIntegration *theRule = theProblem->ruleEdge;
  femDiscrete *theSpace = theProblem->spaceEdge;
  femGeo *theGeometry = theProblem->geometry;
  femNodes *theNodes = theGeometry->theNodes;
  femMesh *theEdges = theGeometry->theEdges;
  double x[2], y[2], phi[2];
  int iBnd, iElem, iInteg, iEdge, i, j, d, map[2], mapU[2];
  int nLocal = 2;
  double *B = theSystem->B;

  for (iBnd = 0; iBnd < theProblem->nBoundaryConditions; iBnd++) {    
    femBoundaryCondition *theCondition = theProblem->conditions[iBnd];
    femBoundaryType type = theCondition->type;
    double value = theCondition->value1;
    
    int *elem = theCondition->domain->elem;
        int nElem = theCondition->domain->nElem;
        if(type == NEUMANN_X || type == NEUMANN_Y){
            for (int e=0; e<nElem; e++) {
                for (int j=0; j<nLocal; j++) {
                    int node = theCondition->domain->mesh->elem[2*elem[e]+j];
                    if(type == NEUMANN_X){
                        map[j] = node;
                        mapU[j] = 2*map[j];
                        x[j] = theNodes->X[map[j]];
                        y[j] = theNodes->Y[map[j]]; 
                    }
                    else if(type == NEUMANN_Y){
                        map[j] = node ;
                        mapU[j] = 2 * map[j] + 1;
                        x[j] = theNodes->X[map[j]];
                        y[j] = theNodes->Y[map[j]]; 
                    }
                }
                double jac = 0.0;
                for (iInteg=0; iInteg < theRule->n; iInteg++) {
                    double xsi = theRule->xsi[iInteg];
                    double weight = theRule->weight[iInteg];
                    femDiscretePhi(theSpace,xsi,phi);
                    jac = sqrt(pow(x[1]-x[0], 2) + pow(y[1]-y[0], 2))/2.0;
                    for (i = 0; i < nLocal; i++) {
                        B[mapU[i]] += phi[i] * value * jac * weight; 
                    }
                }
            }     
        }     

    //
    // A completer :-)
    // Attention, pour le normal tangent on calcule la normale (sortante) au SEGMENT, surtout PAS celle de constrainedNodes
    //
    // Une petite aide pour le calcul de la normale :-)
    // double tx = theNodes->X[node1] - theNodes->X[node0];
    // double ty = theNodes->Y[node1] - theNodes->Y[node0];
    // double nx = ty;
    // double ny = -tx;
  }
}

void femElasticityApplyDirichlet(femProblem *theProblem) {
  femFullSystem *theSystem = theProblem->system;
  femGeo *theGeometry = theProblem->geometry;
  femNodes *theNodes = theGeometry->theNodes;

  for (int node = 0; node < theNodes->nNodes; node++) {
    femConstrainedNode *theConstrainedNode = &theProblem->constrainedNodes[node];
    if (theConstrainedNode->type == UNDEFINED) continue;
    femBoundaryType type = theConstrainedNode->type;

    if (type == DIRICHLET_X) {
      double value = theConstrainedNode->value1;
      femFullSystemConstrain(theSystem, 2 * node, value);
    }
    if (type == DIRICHLET_Y) {
      double value = theConstrainedNode->value1;
      femFullSystemConstrain(theSystem, 2 * node + 1, value);
    }
    if (type == DIRICHLET_XY) {
      double value_x = theConstrainedNode->value1;
      double value_y = theConstrainedNode->value2;
      femFullSystemConstrain(theSystem, 2 * node, value_x);
      femFullSystemConstrain(theSystem, 2 * node + 1, value_y);
    }

    if (type == DIRICHLET_N) {
      double value = theConstrainedNode->value1;
      double nx = theConstrainedNode->nx;
      double ny = theConstrainedNode->ny;
      femFullSystemConstrain(theSystem, 2 * node, value * nx);
      femFullSystemConstrain(theSystem, 2 * node + 1, value * ny);
    }
    if (type == DIRICHLET_T) {
      double value = theConstrainedNode->value1;
      double nx = theConstrainedNode->nx;
      double ny = theConstrainedNode->ny;
      double tx = -ny;
      double ty = nx;
      femFullSystemConstrain(theSystem, 2 * node, value * tx);
      femFullSystemConstrain(theSystem, 2 * node + 1, value * ty);
    }
    if (type == DIRICHLET_NT) {
      double value_n = theConstrainedNode->value1;
      double value_t = theConstrainedNode->value2;
      double nx = theConstrainedNode->nx;
      double ny = theConstrainedNode->ny;
      double tx = -ny;
      double ty = nx;
      femFullSystemConstrain(theSystem, 2 * node, value_n * nx + value_t * tx);
      femFullSystemConstrain(theSystem, 2 * node + 1, value_n * ny + value_t * ty);
    }
  }
}

double *femElasticitySolve(femProblem *theProblem) {
  femElasticityAssembleElements(theProblem);
  femElasticityAssembleNeumann(theProblem);
  femElasticityApplyDirichlet(theProblem);

  // Implémentation du solver CG (conjugate gradient) {fait} + CSR (compressed sparse row) {à implémenter}

  // double* soluce = femFullSystemEliminate(theProblem->system);
  double* soluce = cgSolver(theProblem->system);
  memcpy(theProblem->soluce, soluce, theProblem->system->size * sizeof(double));
  return theProblem->soluce;
}
