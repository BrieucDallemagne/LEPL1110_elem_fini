#include <stdio.h>
#include <math.h>
#include "glfem.h"

double area(double x[3], double y[3]) {
    return 0.5 * fabs((x[0] - x[2]) * (y[1] - y[2]) - (x[1] - x[2]) * (y[0] - y[2]));
}
double integrate(double x[3], double y[3], double (*f) (double, double))
{
    double Jacobian = 2 * area(x, y);
    
    if ((x[0]-x[1])*(y[0]-y[2]) - (x[0]-x[2])*(y[0]-y[1]) > 0){
        Jacobian = (x[0]-x[1])*(y[0]-y[2]) - (x[0]-x[2])*(y[0]-y[1]);
    }
    else{
        Jacobian = (x[0]-x[2])*(y[0]-y[1]) - (x[0]-x[1])*(y[0]-y[2]);
    }

    double I = 0.0;
    double xLoc[3] = {0.0};
    double yLoc[3] = {0.0};

    double xhi[3] = {1/6.0, 1/6.0, 2/3.0};
    double nu[3] = {1/6.0, 2/3.0, 1/6.0};
    double w[3] = {1/6.0, 1/6.0, 1/6.0};

    double phi[3];
    
    for(int k = 0; k<3; k++){
        xLoc[k] = 0.0;
        yLoc[k] = 0.0;
        phi[0] = 1 - xhi[k] - nu[k];
        phi[1] = xhi[k];
        phi[2] = nu[k];

        for(int i = 0; i<3; i++){
            xLoc[k] += x[i] * phi[i];
            yLoc[k] += y[i] * phi[i];
        }

        I += w[k] * f(xLoc[k], yLoc[k]);
    }

    I *= Jacobian;
    
    glfemSetColor(GLFEM_BLACK); glfemDrawElement(x,y,3);
    glfemSetColor(GLFEM_BLUE);  glfemDrawNodes(x,y,3);
    glfemSetColor(GLFEM_RED);   glfemDrawNodes(xLoc,yLoc,3);
    
    return I;
}

double integrateRecursive(double x[3], double y[3], double (*f)(double , double),int n) {

    int points[4][3] = {{0 , 3 , 5} , {3 , 4 , 5}, {3 , 1 , 4} , {5 , 4 , 2}};
    double xi[6] = {0.0 ,1.0 ,0.0 ,0.5 ,0.5 ,0.0};
    double nu[6] = {0.0 ,0.0 ,1.0 ,0.0 ,0.5 ,0.5};
    double sum = 0.0;
    double xLoc[3];
    double yLoc[3];

    if (n <= 0) {
        return integrate(x, y, f);
    }
    else {
        for(int i = 0; i < 4; i++){

            for(int j = 0; j < 3; j++){
                xLoc[j] =  (x[0] *(1 - xi[points[i][j]] - nu[points[i][j]])) +(x[1]* nu[points[i][j]] )  +( x[2] *xi[points[i][j]] );
                yLoc[j] = (y[0]* (1 - xi[points[i][j]] - nu[points[i][j]]) )+( y[1]*  nu[points[i][j]]) + ( y[2]*xi[points[i][j]] );
            }

            sum += integrateRecursive(xLoc, yLoc, f, n-1);
        }
        return sum;
    }
}
