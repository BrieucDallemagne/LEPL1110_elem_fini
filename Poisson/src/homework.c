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
    femPoissonFindBoundaryNodes(theProblem);
 
    if (theSpace->n > 4) Error("Unexpected discrete space size !");  
    double x[4],y[4],phi[4],dphidxsi[4],dphideta[4],dphidx[4],dphidy[4];
    int iElem,iInteg,iEdge,i,j,map[4];
    int nLocal = theMesh->nLocalNode;

    double **A = theSystem->A;
    double *B = theSystem->B;
    int size = theSystem->size;


    for (iElem = 0; iElem < theMesh->nElem; iElem++) {
        femPoissonLocal(theProblem,iElem,map,x,y);
        for (iInteg =0; iInteg < theRule->n; iInteg++){
            double xsi = theRule->xsi[iInteg];
            double eta = theRule->eta[iInteg];
            double weight = theRule->weight[iInteg];
            femDiscretePhi2(theSpace,xsi,eta, phi);
            femDiscreteDphi2 ( theSpace , xsi , eta , dphidxsi , dphideta ) ;
            double dxdxsi = 0.0, dxdeta = 0.0, dydxsi = 0.0, dydeta = 0.0;
            double xlocal = 0.0, ylocal = 0.0;
            for (i=0; i < theSpace->n; i++) {
                xlocal += phi[i]*x[i];
                ylocal += phi[i]*y[i];
                dxdxsi += dphidxsi[i]*x[i];
                dxdeta += dphideta[i]*x[i];
                dydxsi += dphidxsi[i]*y[i];
                dydeta += dphideta[i]*y[i];

            }
            double jacobian = fabs(dxdxsi*dydeta - dxdeta*dydxsi);
            for (i=0; i < theSpace->n; i++) {
                dphidx[i] = (dphidxsi[i]*dydeta - dphideta[i]*dydxsi);
                dphidy[i] = (-dphidxsi[i]*dxdeta + dphideta[i]*dxdxsi);
            }
            double wOverJ = weight/jacobian;
            double wJack = weight*jacobian;
            for ( i = 0; i < theSpace -> n ; i ++) {
                for ( j = 0; j < theSpace -> n ; j ++) {
                    A[map[i]][map[j]] += wOverJ*(dphidx[i]*dphidx[j] + dphidy[i]*dphidy[j]);
                }
                B[map[i]] += wJack*phi[i];
            }


        }
    }
    for (i = 0; i < theBoundary->nElem; i++) {
        iElem = theBoundary->elem[i];
        femPoissonLocal(theProblem,iElem,map,x,y);
        for (iEdge = 0; iEdge < theMesh->nLocalNode; iEdge++) {
            int iNode = map[iEdge];
            B[iNode] = 0.0;
            for (j = 0; j < size; j++) {
                A[iNode][j] = 0.0;
                A[j][iNode] = 0.0;
            }
            A[iNode][iNode] = 1.0;
        }
    }
    
    femFullSystemEliminate(theSystem);

}

# endif



