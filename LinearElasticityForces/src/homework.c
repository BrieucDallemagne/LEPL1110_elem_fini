#include "fem.h"


void femElasticityAssembleElements(femProblem *theProblem){
    femFullSystem  *theSystem = theProblem->system;
    femIntegration *theRule = theProblem->rule;
    femDiscrete    *theSpace = theProblem->space;
    femGeo         *theGeometry = theProblem->geometry;
    femNodes       *theNodes = theGeometry->theNodes;
    femMesh        *theMesh = theGeometry->theElements;
    femMesh        *theEdges = theGeometry->theEdges;
    double x[4],y[4],phi[4],dphidxsi[4],dphideta[4],dphidx[4],dphidy[4];
    int iElem,iInteg,iEdge,i,j,d,map[4],mapX[4],mapY[4];
    int nLocal = theMesh->nLocalNode;
    double a   = theProblem->A;
    double b   = theProblem->B;
    double c   = theProblem->C;      
    double rho = theProblem->rho;
    double g   = theProblem->g;
    double **A = theSystem->A;
    double *B  = theSystem->B;
    
    
    for (iElem = 0; iElem < theMesh->nElem; iElem++) {
        for (j=0; j < nLocal; j++) {
            map[j]  = theMesh->elem[iElem*nLocal+j];
            mapX[j] = 2*map[j];
            mapY[j] = 2*map[j] + 1;
            x[j]    = theNodes->X[map[j]];
            y[j]    = theNodes->Y[map[j]];
        } 
        
        for (iInteg=0; iInteg < theRule->n; iInteg++) {    
            double xsi    = theRule->xsi[iInteg];
            double eta    = theRule->eta[iInteg];
            double weight = theRule->weight[iInteg];  
            femDiscretePhi2(theSpace,xsi,eta,phi);
            femDiscreteDphi2(theSpace,xsi,eta,dphidxsi,dphideta);
            
            double dxdxsi = 0.0;
            double dxdeta = 0.0;
            double dydxsi = 0.0; 
            double dydeta = 0.0;
            for (i = 0; i < theSpace->n; i++) {  
                dxdxsi += x[i]*dphidxsi[i];       
                dxdeta += x[i]*dphideta[i];   
                dydxsi += y[i]*dphidxsi[i];   
                dydeta += y[i]*dphideta[i]; 
            }
            double jac = fabs(dxdxsi * dydeta - dxdeta * dydxsi);
            
            for (i = 0; i < theSpace->n; i++) {    
                dphidx[i] = (dphidxsi[i] * dydeta - dphideta[i] * dydxsi) / jac;       
                dphidy[i] = (dphideta[i] * dxdxsi - dphidxsi[i] * dxdeta) / jac; }            
            for (i = 0; i < theSpace->n; i++) { 
                for(j = 0; j < theSpace->n; j++) {
                    A[mapX[i]][mapX[j]] += (dphidx[i] * a * dphidx[j] + 
                                            dphidy[i] * c * dphidy[j]) * jac * weight;                                                                                            
                    A[mapX[i]][mapY[j]] += (dphidx[i] * b * dphidy[j] + 
                                            dphidy[i] * c * dphidx[j]) * jac * weight;                                                                                           
                    A[mapY[i]][mapX[j]] += (dphidy[i] * b * dphidx[j] + 
                                            dphidx[i] * c * dphidy[j]) * jac * weight;                                                                                            
                    A[mapY[i]][mapY[j]] += (dphidy[i] * a * dphidy[j] + 
                                            dphidx[i] * c * dphidx[j]) * jac * weight; }}
             for (i = 0; i < theSpace->n; i++) {
                B[mapY[i]] -= phi[i] * g * rho * jac * weight; }}} 
}


void femElasticityAssembleNeumann(femProblem *theProblem){
    femFullSystem  *theSystem = theProblem->system;
    femIntegration *theRule = theProblem->ruleEdge;
    femDiscrete    *theSpace = theProblem->spaceEdge;
    femGeo         *theGeometry = theProblem->geometry;
    femNodes       *theNodes = theGeometry->theNodes;
    femMesh        *theEdges = theGeometry->theEdges;
    double x[2],y[2],phi[2];
    int iBnd,iElem,iInteg,iEdge,i,j,d,map[2],mapU[2];
    int nLocal = 2;
    double *B  = theSystem->B;

    for(iBnd=0; iBnd < theProblem->nBoundaryConditions; iBnd++){
        femBoundaryCondition *theCondition = theProblem->conditions[iBnd];
        femBoundaryType type = theCondition->type;
        double value = theCondition->value;
        
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
                    jac = sqrt((x[1]-x[0])*(x[1]-x[0])+(y[1]-y[0])*(y[1]-y[0]))/2.0;
                    for (i = 0; i < nLocal; i++) {
                        B[mapU[i]] += phi[i] * value * jac * weight; 
                    }
                }
            }     
        }
        
    }
}

double *femElasticitySolve(femProblem *theProblem) {

    femFullSystem *theSystem = theProblem->system;

    femElasticityAssembleElements(theProblem);
    femElasticityAssembleNeumann(theProblem);
    int *theConstrainedNodes = theProblem->constrainedNodes;     
    for (int i=0; i < theSystem->size; i++) {
        if (theConstrainedNodes[i] != -1) {
            double value = theProblem->conditions[theConstrainedNodes[i]]->value;
            femFullSystemConstrain(theSystem,i,value); 
        }
    }

    memcpy(theProblem->soluce,femFullSystemEliminate(theSystem),theProblem->system->size*sizeof(double));
    return theProblem->soluce;

}

double * femElasticityForces(femProblem *theProblem){            
 
    double **A = theProblem->system->A; 
    double *B = theProblem->system->B;
    double *x = theProblem->soluce;
    double *residuals = theProblem->residuals;
    double size_prob = theProblem->system->size;

    femFullSystemInit(theProblem->system);
    femElasticityAssembleElements(theProblem); 
    femElasticityAssembleNeumann(theProblem);

    for(int i = 0; i < size_prob; i++){
        residuals[i]= 0.0;
        for(int j = 0; j < size_prob; j++){
            residuals[i] += A[i][j] * x[j];
        } 
        residuals[i] -= B[i];
    }
    return residuals;

}