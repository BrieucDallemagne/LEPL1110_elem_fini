#include <stdio.h>
#include <math.h>
#include "glfem.h"


double integrate(double x[3], double y[3], double (*f) (double, double))
{
    double I = 0;

    double Jacobian = 0;
    if ((x[0]-x[1])*(y[0]-y[2]) - (x[0]-x[2])*(y[0]-y[1]) > 0){
        Jacobian = (x[0]-x[1])*(y[0]-y[2]) - (x[0]-x[2])*(y[0]-y[1]);
    }
    else{
        Jacobian = (x[0]-x[2])*(y[0]-y[1]) - (x[0]-x[1])*(y[0]-y[2]);
    }

    double xi[3] = {0.166666666666, 0.16666666666666, 0.666666666666};
    double nu[3] = {0.166666666666666, 0.66666666666666, 0.1666666666666666};
    double omega[3] = {0.16666666666, 0.166666666666666, 0.166666666666666};
    double XLoc[3];
    double YLoc[3];

    for (int i = 0; i < 3; i++){
        double xiLoc = xi[i];
        double nuLoc = nu[i];
        XLoc[i] = integrhelp(x, xiLoc, nuLoc);
        YLoc[i] = integrhelp(y, xiLoc, nuLoc);
        I += f(XLoc[i], YLoc[i]) * omega[i];
    }

    I = I * Jacobian;


  glfemSetColor(GLFEM_BLACK); glfemDrawElement(x,y,3);
  glfemSetColor(GLFEM_BLUE);  glfemDrawNodes(x,y,3);
//  glfemSetColor(GLFEM_RED);   glfemDrawNodes(xLoc,yLoc,3);
    


    return I;
}

double integrhelp(double var[3], double xi, double nu){
    double I = 0;
    I = xi*var[0] + nu*var[1] + (1-xi-nu)*var[2];
    return I;
}

double integrateRecursive(double x[3], double y[3], double (*f)(double,double), int n)
{

    if(n == 0){
        return integrate(x, y, f);
    }
    else{
        double I = 0;
        double x1[3] = {x[0],(x[0]+x[1])/2 ,(x[0]+x[2])/2};
        double x2[3] = {(x[0]+x[1])/2,x[1],(x[1]+x[2])/2};
        double x3[3] = {(x[0]+x[2])/2, (x[1]+x[2])/2, x[2]};
        double y1[3] = {y[0],(y[0]+y[1])/2, (y[0]+y[2])/2};
        double y2[3] = {(y[0]+y[1])/2, y[1], (y[1]+y[2])/2};
        double y3[3] = {(y[0]+y[2])/2,  (y[1]+y[2])/2, y[2]};
 
        I += integrateRecursive(x1, y1, f, n-1);
        I += integrateRecursive(x2, y2, f, n-1);
        I += integrateRecursive(x3, y3, f, n-1);
        return I;
    }

}
