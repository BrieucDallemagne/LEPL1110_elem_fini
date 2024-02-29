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




    for (int nbElem = 0; nbElem < (theMesh->nElem); nbElem++){

        femPoissonLocal(theProblem,nbElem,map,x,y);

        for (int iInt = 0; iInt < (theRule -> n); iInt++){

            //on récup les fonctions de formes 
            double xsi = theRule->xsi[iInt];
            double eta = theRule->eta[iInt];
            double weight = theRule->weight[iInt];

            femDiscretePhi2(theSpace,xsi,eta,phi);
            femDiscreteDphi2(theSpace,xsi,eta,dphidxsi,dphideta);

            double dxdxsi = 0;
            double dxdeta = 0;
            double dydxsi = 0;
            double dydeta = 0;
            
            double xloc = 0;
            double yloc = 0;

            for (int i = 0; i < (theSpace->n); i++){

                xloc += x[i] * phi[i];
                yloc += y[i] * phi[i];
                dxdxsi += x[i] * dphidxsi[i];
                dxdeta += x[i] * dphideta[i];
                dydxsi += y[i] * dphidxsi[i];
                dydeta += y[i] * dphideta[i];
            }

           
            double jacobian = dxdxsi * dydeta - dxdeta * dydxsi;
            if(jacobian < 0){
                if(nLocal == 3){
                    int temp = theMesh->elem[nLocal*nbElem];
                    theMesh->elem[nLocal*nbElem] = theMesh->elem[nLocal*nbElem + 1];
                    theMesh->elem[nLocal*nbElem + 1] = temp;
                }else{
                    int temp = theMesh->elem[nLocal*nbElem];
                    theMesh->elem[nLocal*nbElem] = theMesh->elem[nLocal*nbElem + 1];
                    theMesh->elem[nLocal*nbElem + 1] = temp;
                    int temp2 = theMesh->elem[nLocal*nbElem+2];
                    theMesh->elem[nLocal*nbElem+2] = theMesh->elem[nLocal*nbElem + 3];
                    theMesh->elem[nLocal*nbElem + 3] = temp2;
                }
            }
            jacobian = fabs(jacobian);

            
            for (int i = 0; i < (theSpace->n); i++){

                dphidx[i] = (dphidxsi[i] * dydeta - dphideta[i] * dydxsi)/jacobian;
                dphidy[i] = (dphideta[i] * dxdxsi - dphidxsi[i] * dxdeta)/jacobian;
            }
 
            double Jw = jacobian * weight;    

            for (int i = 0; i < (theSpace->n); i++){
                for (int j = 0; j < (theSpace->n); j++){
                    theSystem -> A[map[i]][map[j]] += (dphidx[i] * dphidx[j] + dphidy[i] * dphidy[j]) * Jw ;
                }
            }
            for (int i = 0; i < (theSpace->n); i++){    
                theSystem -> B[map[i]] += phi[i] * Jw;
            }
        }
    }
    
    for (int nbEdge = 0; nbEdge < (theBoundary->nElem ); nbEdge++){
        if (theBoundary->elem[nbEdge] == -1){
            for(int i = 0; i < 2; i++){
                int node = theBoundary->mesh->elem[nbEdge*2+i];
                femFullSystemConstrain(theSystem, node, 0.0);
            }
        }
    }
    
    femFullSystemEliminate(theSystem);












    
}

# endif



