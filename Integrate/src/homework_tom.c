#include <stdio.h>
#include <math.h>
#include "glfem.h"


double integrate(double x[3], double y[3], double (*f) (double, double))
{
    double I = 0;
    double xLoc[3];
    double yLoc[3];

    double xhi[3] = {1/6.0, 1/6.0, 2/3.0};
    double nu[3] = {1/6.0, 2/3.0, 1/6.0};
    double w[3] = {1/6.0, 1/6.0, 1/6.0};

    double phi[3];
    
    for(int j = 0; j<3; j++){

        for(int k = 0; k<3; k++){
            
            phi[0] = 1 - xhi[k] - nu[k];
            phi[1] = xhi[k];
            phi[2] = nu[k];

            for(int i = 0; i<3; i++){
                xLoc[k] += x[i] * phi[i];
                yLoc[k] += y[i] * phi[i];
            }
        }

        I += w[j] * f(xLoc[j], yLoc[j]);
    }

    glfemSetColor(GLFEM_BLACK); glfemDrawElement(x,y,3);
    glfemSetColor(GLFEM_BLUE);  glfemDrawNodes(x,y,3);
    glfemSetColor(GLFEM_RED);   glfemDrawNodes(xLoc,yLoc,3);

    return I;

//
// ... A modifier :-)
//
//
// Pour dessiner l'element, les sommets du triangle :-)
// Decommenter la ligne pour dessiner aussi les points d'integration
//


    
}

double integrateRecursive(double x[3], double y[3], double (*f)(double,double), int n)
{

//
// ... A modifier :-)
// y-compris la ligne juste en dessous :-)
//
    double I = integrate(x,y,f);
    
//
//
//    
     
    return I;
}
