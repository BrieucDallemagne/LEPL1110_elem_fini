#include "fem.h"

# ifndef NOPOISSONCREATE

femPoissonProblem *femPoissonCreate(const char *filename)
{
    femGeo* theGeometry = geoMeshCreate(filename);
    femPoissonProblem *theProblem = malloc(sizeof(femPoissonProblem));
    theProblem->geo  = theGeometry;
    femMesh *theMesh = theGeometry->theElements;
    if (theMesh->nLocalNode == 4) {
        theProblem->space = femDiscreteCreate(4,FEM_QUAD);
        theProblem->rule = femIntegrationCreate(4,FEM_QUAD); }
    else if (theMesh->nLocalNode == 3) {
        theProblem->space = femDiscreteCreate(3,FEM_TRIANGLE);
        theProblem->rule = femIntegrationCreate(3,FEM_TRIANGLE); }
    theProblem->system = femFullSystemCreate(theMesh->nodes->nNodes);
    return theProblem;
}

# endif
# ifndef NOPOISSONBOUNDARY

void femPoissonFindBoundaryNodes(femPoissonProblem *theProblem)
{
    femGeo* theGeometry = theProblem->geo;  
    femMesh* theEdges = theGeometry->theEdges; 
    int nBoundary = 0;

    for (int i = 0; i < theEdges->nElem; i++) {
        if (theEdges->elem[i] == 0) {
            nBoundary++;
        }
    }

    femDomain *theBoundary = malloc(sizeof(femDomain));
    theGeometry->nDomains++;
    theGeometry->theDomains = realloc(theGeometry->theDomains,theGeometry->nDomains*sizeof(femDomain*));
    theGeometry->theDomains[theGeometry->nDomains-1] = theBoundary;
    theBoundary->nElem = nBoundary;
    theBoundary->elem = malloc(nBoundary*sizeof(int));
    theBoundary->mesh = NULL;
    sprintf(theBoundary->name,"Boundary");

    for(int i = 0; i < nBoundary; i++) {
        theBoundary->elem[i] = 0;
    }
 


}
    
# endif
# ifndef NOPOISSONFREE

void femPoissonFree(femPoissonProblem *theProblem)
{

    geoMeshFree(theProblem->geo);
    femDiscreteFree(theProblem->space);
    femIntegrationFree(theProblem->rule);
    femFullSystemFree(theProblem->system);
    free(theProblem);
}
    
# endif
# ifndef NOPOISSONLOCAL

void femPoissonLocal(femPoissonProblem *theProblem, const int iElem, int *map, double *x, double *y)
{
    femMesh *theMesh = theProblem->geo->theElements;
    femDiscrete *theSpace = theProblem->space;
    for (int i = 0; i < theMesh->nLocalNode; i++) {
        map[i] = theMesh->elem[iElem*theMesh->nLocalNode+i];
        x[i] = theMesh->nodes->X[map[i]];
        y[i] = theMesh->nodes->Y[map[i]];
    }
}

# endif
# ifndef NOPOISSONSOLVE

void femPoissonSolve(femPoissonProblem *theProblem)
{

    femMesh *theMesh = theProblem->geo->theElements;
    femDomain *theBoundary = geoGetDomain(theProblem->geo,"Boundary");
    femFullSystem *theSystem = theProblem->system;
    femIntegration *theRule = theProblem->rule;
    femDiscrete *theSpace = theProblem->space;
 
    if (theSpace->n > 4) Error("Unexpected discrete space size !");  
    double x[4],y[4],phi[4],dphidxsi[4],dphideta[4],dphidx[4],dphidy[4];
    int iElem,iInteg,iEdge,i,j,map[4];
    int nLocal = theMesh->nLocalNode;

    double **A = theSystem->A;
    double *B = theSystem->B;
    int size = theSystem->size;

    for (int Elem; Elem < theMesh->nElem; Elem++) {
        femPoissonLocal(theProblem,Elem,map,x,y);
        for (int i = 0; i < theRule->n; i++) {
            iInteg = i;
            theSpace->phi2(theRule->xsi[iInteg],theRule->eta[iInteg],phi);
            theSpace->dphi2dx(theRule->xsi[iInteg],theRule->eta[iInteg],dphidxsi,dphideta);
            double dxdxsi = 0.0, dydxsi = 0.0, dxdeta = 0.0, dydeta = 0.0;
            for (int j = 0; j < nLocal; j++) {
                dxdxsi += dphidxsi[j]*x[j];
                dydxsi += dphidxsi[j]*y[j];
                dxdeta += dphideta[j]*x[j];
                dydeta += dphideta[j]*y[j];
            }
            double jac = dxdxsi*dydeta - dydxsi*dxdeta;
            double invJac = 1.0/jac;
            double dphidx[4],dphidy[4];
            for (int j = 0; j < nLocal; j++) {
                dphidx[j] = invJac*(dphidxsi[j]*dydeta - dphideta[j]*dydxsi);
                dphidy[j] = invJac*(dphideta[j]*dxdxsi - dphidxsi[j]*dxdeta);
            }
            double **a = malloc(nLocal*sizeof(double*));
            for (int j = 0; j < nLocal; j++) {
                a[j] = malloc(nLocal*sizeof(double));
            }
            double *b = malloc(nLocal*sizeof(double));
            for (int j = 0; j < nLocal; j++) {
                b[j] = 0.0;
                for (int k = 0; k < nLocal; k++) {
                    a[j][k] = 0.0;
                }
            }
            for (int j = 0; j < nLocal; j++) {
                for (int k = 0; k < nLocal; k++) {
                    a[j][k] = a[j][k] + (dphidx[j]*dphidx[k] + dphidy[j]*dphidy[k])*jac*theRule->weight[iInteg];
                }
                b[j] = b[j] + 1.0*phi[j]*jac*theRule->weight[iInteg];
            }
            for (int j = 0; j < nLocal; j++) {
                for (int k = 0; k < nLocal; k++) {
                    A[map[j]][map[k]] = A[map[j]][map[k]] + a[j][k];
                }
                B[map[j]] = B[map[j]] + b[j];
            }
            for (int j = 0; j < nLocal; j++) {
                free(a[j]);
            }
            free(a);
            free(b);
        }

    }

    for (int i = 0; i < theBoundary->nElem; i++) {
        iElem = theBoundary->elem[i];
        for (int j = 0; j < nLocal; j++) {
            B[theMesh->elem[iElem*nLocal+j]] = 0.0;
            for (int k = 0; k < size; k++) {
                A[theMesh->elem[iElem*nLocal+j]][k] = 0.0;
                A[k][theMesh->elem[iElem*nLocal+j]] = 0.0;
            }
            A[theMesh->elem[iElem*nLocal+j]][theMesh->elem[iElem*nLocal+j]] = 1.0;
        }
    }


    
}

# endif



