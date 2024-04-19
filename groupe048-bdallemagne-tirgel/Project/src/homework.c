
#include "fem.h"

// Il faut un fifrelin generaliser ce code.....
//  (1) Ajouter l'axisymétrique ! (mandatory)
//  (2) Ajouter les conditions de Neumann ! (mandatory)
//  (3) Ajouter les conditions en normal et tangentiel !   (strongly advised)
//  (4) Et remplacer le solveur plein par un truc un fifrelin plus subtil  (mandatory)

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
  double r = 0.0;

  for (iElem = 0; iElem < theMesh->nElem; iElem++) {
    for (j = 0; j < nLocal; j++) {
      map[j] = theMesh->elem[iElem * nLocal + j];
      mapX[j] = 2 * map[j];
      mapY[j] = 2 * map[j] + 1;
      x[j] = theNodes->X[map[j]];
      y[j] = theNodes->Y[map[j]];
    }

    for (iInteg = 0; iInteg < theRule->n; iInteg++) {
      r = 0.0;
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
        r += x[i] * phi[i];
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
          if (theProblem->planarStrainStress == AXISYM){
             A[mapX[i]][mapX[j]] += (dphidx[i] * a * dphidx[j] * r + dphidy[i] * c * dphidy[j] * r + dphidx[i] * b * phi[j] + phi[i] * (b * dphidx[j] + a * phi[j]/r)) * jac * weight;
             A[mapX[i]][mapY[j]] += (dphidx[i] * b * dphidy[j] * r + dphidy[i] * c * dphidx[j] * r + phi[i] * b * dphidy[j]) * jac * weight;
             A[mapY[i]][mapX[j]] += (dphidy[i] * b * dphidx[j] * r + dphidx[i] * c * dphidy[j] * r + dphidy[i] * b * phi[j]) * jac * weight;
             A[mapY[i]][mapY[j]] += (dphidy[i] * a * dphidy[j] * r + dphidx[i] * c * dphidx[j] * r) * jac * weight;
          }else{
            A[mapX[i]][mapX[j]] += (dphidx[i] * a * dphidx[j] + dphidy[i] * c * dphidy[j]) * jac * weight;
            A[mapX[i]][mapY[j]] += (dphidx[i] * b * dphidy[j] + dphidy[i] * c * dphidx[j]) * jac * weight;
            A[mapY[i]][mapX[j]] += (dphidy[i] * b * dphidx[j] + dphidx[i] * c * dphidy[j]) * jac * weight;
            A[mapY[i]][mapY[j]] += (dphidy[i] * a * dphidy[j] + dphidx[i] * c * dphidx[j]) * jac * weight;
          }
        }
      }
      for (i = 0; i < theSpace->n; i++) {
        if (theProblem->planarStrainStress == AXISYM){
          B[mapX[i]] += phi[i] * gx * rho * r * jac * weight;
          B[mapY[i]] += phi[i] * gy * rho * r * jac * weight;
        }else{
          B[mapX[i]] += phi[i] * gx * rho * jac * weight;
          B[mapY[i]] += phi[i] * gy * rho * jac * weight;
        }
      }
    }
  }
}

void femElasticityAssembleNeumann(femProblem *theProblem) {
    femFullSystem  *theSystem = theProblem->system;
    femIntegration *theRule = theProblem->ruleEdge;
    femDiscrete    *theSpace = theProblem->spaceEdge;
    femGeo         *theGeometry = theProblem->geometry;
    femNodes       *theNodes = theGeometry->theNodes;
    femMesh        *theEdges = theGeometry->theEdges;
    double x[2],y[2],phi[2];
    int iBnd,iElem,iInteg,iEdge,i,j,d,map[2],mapU[2], mapUx[2], mapUy[2];
    int nLocal = 2;
    double *B  = theSystem->B;
    double r;

    for (iBnd = 0; iBnd < theProblem->nBoundaryConditions; iBnd++) {    
      femBoundaryCondition *theCondition = theProblem->conditions[iBnd];
      femBoundaryType type = theCondition->type;
      double value = theCondition->value1;
      int *elem = theCondition->domain->elem;
      int nElem = theCondition->domain->nElem;

      if(type == NEUMANN_X || type == NEUMANN_Y){
        for (int e=0; e<nElem; e++) {
          for (int j=0; j<nLocal; j++) {
            int node = theCondition->domain->mesh->elem[nLocal*elem[e]+j];
            map[j] = node;
            if(type == NEUMANN_X){
                mapU[j] = 2*map[j];
            }
            else if(type == NEUMANN_Y){
                mapU[j] = 2 * map[j] + 1;
            }
            x[j] = theNodes->X[map[j]];
            y[j] = theNodes->Y[map[j]]; 
          }
          double jac = sqrt(pow(x[1]-x[0], 2) + pow(y[1]-y[0], 2))/2.0;
          for (iInteg=0; iInteg < theRule->n; iInteg++) {
            r = 0.0; 
            for(int i=0; i<theSpace->n; i++) r += x[i]*phi[i];
            double xsi = theRule->xsi[iInteg];
            double weight = theRule->weight[iInteg];
            femDiscretePhi(theSpace,xsi,phi);
            if (type == AXISYM){
              for (i = 0; i < theSpace->n; i++) {
                B[mapU[i]] += phi[i] * value * r * jac * weight;
              }
            }else{
              for (i = 0; i < theSpace->n; i++) {
                B[mapU[i]] += phi[i] * value * jac * weight;
              }
            }
          }
        }     
      } 
      if(type == NEUMANN_N || type == NEUMANN_T){
          for (int e=0; e<nElem; e++) {
              for (int j=0; j<nLocal; j++) {
                int node = theCondition->domain->mesh->elem[nLocal*elem[e]+j];
                map[j] = node;
                mapUx[j] = 2*map[j];           
                mapUy[j] = 2 * map[j] + 1;
                x[j] = theNodes->X[map[j]];
                y[j] = theNodes->Y[map[j]]; 
              }
              
              double h = sqrt(pow(x[1]-x[0], 2) + pow(y[1]-y[0], 2));
              double jac = h/2.0;

              double tx = (x[1] - x[0]); 
              double ty = (y[1] - y[0]);
              
              for (iInteg=0; iInteg < theRule->n; iInteg++) {
                  r = 0.0;
                  for(int i=0; i<theSpace->n; i++) r += x[i]*phi[i];
                  double xsi = theRule->xsi[iInteg];
                  double weight = theRule->weight[iInteg];
                  femDiscretePhi(theSpace,xsi,phi); 
 
                  if(type == NEUMANN_T){
                    if(type == AXISYM){
                      for (i = 0; i < theSpace->n; i++) {
                        B[mapUx[i]] += phi[i] * value * jac * weight * tx * r; 
                        B[mapUy[i]] += phi[i] * value * jac * weight * ty * r; 
                      }  
                    }else{
                      for (i = 0; i < theSpace->n; i++) {
                        B[mapUx[i]] += phi[i] * value * jac * weight * tx; 
                        B[mapUy[i]] += phi[i] * value * jac * weight * ty; 
                      }
                    }
                  }
                  if(type == NEUMANN_N){
                    double nx = ty/h;
                    double ny = -tx/h;
                    if(type == AXISYM){
                      for (i = 0; i < theSpace->n; i++) {
                        B[mapUx[i]] += phi[i] * value * weight * jac * nx * r; 
                        B[mapUy[i]] += phi[i] * value * weight * jac * ny * r; 
                      }
                    }else{
                      for (i = 0; i < theSpace->n; i++) {
                        B[mapUx[i]] += phi[i] * value * weight * jac * nx; 
                        B[mapUy[i]] += phi[i] * value * weight * jac * ny;
                      }
                    }
                  }
              }
          }     
      }         
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
      femFullSystemConstrain(theSystem, 2 * node + 0, value);
    }
    if (type == DIRICHLET_Y) {
      double value = theConstrainedNode->value1;
      femFullSystemConstrain(theSystem, 2 * node + 1, value);
    }
    if (type == DIRICHLET_XY) {
      double value_x = theConstrainedNode->value1;
      double value_y = theConstrainedNode->value2;
      femFullSystemConstrain(theSystem, 2 * node + 0, value_x);
      femFullSystemConstrain(theSystem, 2 * node + 1, value_y);
    }

    if (type == DIRICHLET_N) {
      double value = theConstrainedNode->value1;
      double nx = theConstrainedNode->nx;
      double ny = theConstrainedNode->ny;
      femFullSystemConstrainN(theSystem, node, value, nx, ny);
    }
    if (type == DIRICHLET_T) {
      double value = theConstrainedNode->value1;
      double nx = theConstrainedNode->nx;
      double ny = theConstrainedNode->ny;
      femFullSystemConstrainT(theSystem, node, value, nx, ny);
    }
    if (type == DIRICHLET_NT) {
      double value_n = theConstrainedNode->value1;
      double value_t = theConstrainedNode->value2;
      double nx = theConstrainedNode->nx;
      double ny = theConstrainedNode->ny;
      double tx = -ny;
      double ty = nx;
      // femFullSystemConstrainN(theSystem, node, value_n, nx, ny);
      // femFullSystemConstrainT(theSystem, node, value_t, nx, ny);
      femFullSystemConstrain(theSystem, 2 * node + 0, value_n * nx + value_t * tx);
      femFullSystemConstrain(theSystem, 2 * node + 1, value_n * ny + value_t * ty);
    }
  }
}

double *femElasticitySolve(femProblem *theProblem) {
  femElasticityAssembleElements(theProblem);
  femElasticityAssembleNeumann(theProblem);
  femElasticityApplyDirichlet(theProblem);

  // double* soluce = femFullSystemEliminate(theProblem->system);
  double *soluce = cgSolver(theProblem->system);
  memcpy(theProblem->soluce, soluce, theProblem->system->size * sizeof(double));
  return theProblem->soluce;
}
