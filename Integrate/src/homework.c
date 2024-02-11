#include <stdio.h>
#include <math.h>
#include "glfem.h"

double integrhelp(double var[3], double xi, double nu){
    double I = 0;
    I = xi*var[0] + nu*var[1] + (1-xi-nu)*var[2];
    return I;
}

double integrate(double x[3], double y[3], double (*f) (double, double))
{
    double sum;

    double Jacobian;
    if ((x[0]-x[1])*(y[0]-y[2]) - (x[0]-x[2])*(y[0]-y[1]) > 0){
        Jacobian = (x[0]-x[1])*(y[0]-y[2]) - (x[0]-x[2])*(y[0]-y[1]);
    }
    else{
        Jacobian = (x[0]-x[2])*(y[0]-y[1]) - (x[0]-x[1])*(y[0]-y[2]);
    }

    const double xi[3] = {1.0/6.0, 2.0/3.0,1/6.0 };
    const double nu[3] = {1/6.0, 1/6.0,2/3.0};
    const double omega = 1/6.0;

    double xLoc [3]; 
    double yLoc [3];

    for (int i = 0; i < 3; ++i) {
        xLoc[i] =  x[0] * (1 - nu[i] -xi[i]) + nu[i] * x[1] + xi[i] *x[2];
        yLoc[i] =  y[0] * (1 -nu[i] -xi[i]) + nu[i] *y[1] +xi[i] * y[2];

    }

    for (int i = 0; i < 3; ++i) {
        sum += omega * f(xLoc[i], yLoc[i]);
    }

    sum *= Jacobian;

    glfemSetColor(GLFEM_BLACK); glfemDrawElement(x,y,3);
    glfemSetColor(GLFEM_BLUE);  glfemDrawNodes(x,y,3);
    glfemSetColor(GLFEM_RED);   glfemDrawNodes(xLoc,yLoc,3);
    
    return sum;
}



double integrateRecursive(double x[3], double y[3], double (*f)(double , double),int n) {

    int points[4][3] = {{0 , 3 , 5} , {3 , 4 , 5}, {3 , 1 , 4} , {5 , 4 , 2}};
    double xi[6] = {0.0 ,1.0 ,0.0 ,0.5 ,0.5 ,0.0};
    double nu[6] = {0.0 ,0.0 ,1.0 ,0.0 ,0.5 ,0.5};
    double sum;
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
