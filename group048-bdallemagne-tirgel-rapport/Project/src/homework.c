
#include "fem.h"

// Il faut un fifrelin generaliser ce code.....
//  (1) Ajouter l'axisymétrique !    (mandatory)
//  (2) Ajouter les conditions de Neumann !   (mandatory)
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
          if (theProblem->planarStrainStress = AXISYM){
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
        if (theProblem->planarStrainStress = AXISYM){
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
    int iBnd,iElem,iInteg,iEdge,i,j,d,map[2],mapU[2], mapX[2], mapY[2];
    int nLocal = 2;
    double *B  = theSystem->B;

  for (iBnd = 0; iBnd < theProblem->nBoundaryConditions; iBnd++) {    
    femBoundaryCondition *theCondition = theProblem->conditions[iBnd];
    femBoundaryType type = theCondition->type;
    double value = theCondition->value1;
    int *elem = theCondition->domain->elem;
    int nElem = theCondition->domain->nElem;
        if(type == NEUMANN_X || type == NEUMANN_Y){
        printf("nElem = %d\n", nElem);
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
            double jac = sqrt(pow(x[1]-x[0], 2) + pow(y[1]-y[0], 2))/2.0;
            for (iInteg=0; iInteg < theRule->n; iInteg++) {
                double xsi = theRule->xsi[iInteg];
                double weight = theRule->weight[iInteg];
                femDiscretePhi(theSpace,xsi,phi);
                for (i = 0; i < nLocal; i++) {
                    B[mapU[i]] += phi[i] * value * jac * weight; 
                }
            }
        }     
    } 

    if(type == NEUMANN_N || type == NEUMANN_T){
        for (int e=0; e<nElem; e++) {
            for (int j=0; j<nLocal; j++) {
                int node = theCondition->domain->mesh->elem[2*elem[e]+j];
                map[j] = node;
                mapUx[j] = 2*map[j];           
                mapUy[j] = 2 * map[j] + 1;
                x[j] = theNodes->X[map[j]];
                y[j] = theNodes->Y[map[j]]; 
            }
            double jac = sqrt(pow(x[1]-x[0], 2) + pow(y[1]-y[0], 2))/2.0;
            for (iInteg=0; iInteg < theRule->n; iInteg++) {
                double xsi = theRule->xsi[iInteg];
                double weight = theRule->weight[iInteg];
                femDiscretePhi(theSpace,xsi,phi); 

                double tx = (theNodes->X[map[1]] - theNodes->X[map[0]]) / jac * 2; 
                double ty = (theNodes->Y[map[1]] - theNodes->Y[map[0]]) / jac * 2;

                for (i = 0; i < nLocal; i++) {
                    if(type == NEUMANN_T){
                        B[mapUx[i]] += phi[i] * value * jac * weight * tx; 
                        B[mapUy[i]] += phi[i] * value * jac * weight * ty; 
                    }
                    else if(type == NEUMANN_N){
                        double nx = ty;
                        double ny = -tx;
                        B[mapUx[i]] += phi[i] * value * weight * jac * nx; 
                        B[mapUy[i]] += phi[i] * value * weight * jac * ny; 
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

  double* soluce = femFullSystemEliminate(theProblem->system);
  memcpy(theProblem->soluce, soluce, theProblem->system->size * sizeof(double));
  return theProblem->soluce;
}
