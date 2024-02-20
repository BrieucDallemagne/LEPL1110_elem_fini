#include "fem.h"




double geoSize(double x, double y){

    femGeo* theGeometry = geoGetGeometry();
    
    double h = theGeometry->h;
    double x0 = theGeometry->xNotch;
    double y0 = theGeometry->yNotch;
    double r0 = theGeometry->rNotch;
    double h0 = theGeometry->hNotch;
    double d0 = theGeometry->dNotch;
  
    
    double x1 = theGeometry->xHole;
    double y1 = theGeometry->yHole;
    double r1 = theGeometry->rHole;
    double h1 = theGeometry->hHole;
    double d1 = theGeometry->dHole;

    double hfinal = h;
    double d = sqrt(pow(x-x0,2) + pow(y-y0,2)) - r0;
    if (d < d0) {
        double a = (-2*h + 2*h0)/pow(d0,3);
        double b = (3*h  - 3*h0)/pow(d0,2);
        double c = 0;
        hfinal = a*pow(d,3) + b*pow(d,2) + c*d + h0; }
        
    d = sqrt(pow(x-x1,2) + pow(y-y1,2)) - r1;
    if (d < d1) {
        double a = (-2*h + 2*h1)/(pow(d1,3));
        double b = (3*h  - 3*h1)/(pow(d1,2));
        double c = 0;
        hfinal = fmin(hfinal,a*pow(d,3) + b*pow(d,2) + c*d + h1); }
        
    return hfinal;
    
//   
// Your contribution ends here :-)
//

}


#define ___ 0

void geoMeshGenerate() {

    femGeo* theGeometry = geoGetGeometry();

    double w = theGeometry->LxPlate;
    double h = theGeometry->LyPlate;
     
    double x0 = theGeometry->xNotch;
    double y0 = theGeometry->yNotch;
    double r0 = theGeometry->rNotch;
    
    
    double x1 = theGeometry->xHole;
    double y1 = theGeometry->yHole;
    double r1 = theGeometry->rHole;
 
//
//  -1- Construction de la g�om�trie avec OpenCascade
//      On cr�e le rectangle
//      On cr�e les deux cercles
//      On soustrait les cercles du rectangle :-)
//
 
    int ierr;
    int idPlate = gmshModelOccAddRectangle(-w/2.0,-h/2.0,0.0,w,h,-1,0.0,&ierr);   
    ErrorGmsh(ierr);
    int idNotch = gmshModelOccAddDisk(x0,y0,0.0,r0,r0,-1,NULL,0,NULL,0,&ierr); 
    ErrorGmsh(ierr);
    int idHole  = gmshModelOccAddDisk(x1,y1,0.0,r1,r1,-1,NULL,0,NULL,0,&ierr);    
    ErrorGmsh(ierr);
    
    int *plate = (int*) malloc(2 * sizeof(int));
    int *notch = (int*) malloc(2 * sizeof(int));
    int *hole = (int*) malloc(2 * sizeof(int));
    plate[0] = 2;
    notch[0] = 2;
    hole[0]  = 2;
    plate[1] = idPlate;
    notch[1] = idNotch;
    hole[1]  = idHole;



    gmshModelOccCut(plate,2,notch,2,NULL,NULL,NULL,NULL,NULL,-1,1,1,&ierr); 
    ErrorGmsh(ierr);
    gmshModelOccCut(plate,2,hole,2,NULL,NULL,NULL,NULL,NULL,-1,1,1,&ierr); 
    ErrorGmsh(ierr);
 
//
//  -2- D�finition de la fonction callback pour la taille de r�f�rence
//      Synchronisation de OpenCascade avec gmsh
//      G�n�ration du maillage (avec l'option Mesh.SaveAll :-)
                  
   
    geoSetSizeCallback(geoSize);
                                  
    gmshModelOccSynchronize(&ierr);       
    gmshOptionSetNumber("Mesh.SaveAll", 1, &ierr);
    gmshModelMeshGenerate(2, &ierr);  

    free(plate);
    free(notch);
    free(hole);
       
//
//  Generation de quads :-)
//
//    gmshOptionSetNumber("Mesh.SaveAll", 1, &ierr);
//    gmshOptionSetNumber("Mesh.RecombineAll", 1, &ierr);
//    gmshOptionSetNumber("Mesh.Algorithm", 8, &ierr);  chk(ierr);
//    gmshOptionSetNumber("Mesh.RecombinationAlgorithm", 1.0, &ierr);  chk(ierr);
//    gmshModelGeoMeshSetRecombine(2,1,45,&ierr);  chk(ierr);
//    gmshModelMeshGenerate(2, &ierr);  
   
 
//
//  Plot of Fltk
//
//   gmshFltkInitialize(&ierr);
//   gmshFltkRun(&ierr);  chk(ierr);
//
    
}