#include "fem.h"

//
// Ici, vous pouvez définir votre géométrie :-)
//  (1) Raffiner intelligemment.... (yes )
//  (2) Construire la geometrie avec OpenCascade
//  (3) Construire la geometrie avec les outils de GMSH
//  (4) Obtenir la geometrie en lisant un fichier .geo de GMSH

double geoSize(double x, double y) {

  femGeo *theGeometry = geoGetGeometry();
  return theGeometry->h * (1.0 - 0.5 * x);

}


void geoMeshGenerate(void) {
  femGeo *theGeometry = geoGetGeometry();
  double Lx = 1.0;
  double Ly = 1.0;
  theGeometry->LxPlate = Lx;
  theGeometry->LyPlate = Ly;
  theGeometry->h = Lx * 0.05;
  theGeometry->elementType = FEM_QUAD;
  double lc = theGeometry->h;

  geoSetSizeCallback(geoSize);

  double w = theGeometry->LxPlate;
  double h = theGeometry->LyPlate;

  int ierr;
  double r = w / 4;
  int idRect = gmshModelOccAddRectangle(0.0, 0.0, 0.0, w, h, -1, 0.0, &ierr);

  int idp1 = gmshModelOccAddPoint(0.0, h, 0.0, lc, 1, &ierr);
  int idp2 = gmshModelOccAddPoint(w, 0.0, 0.0, lc, 2, &ierr);
  int idLine1 = gmshModelOccAddLine(idp1, idp2, -1, &ierr);
  //int idDisk = gmshModelOccAddDisk(w / 2.0, h / 2.0, 0.0, r, r, -1, NULL, 0, NULL, 0, &ierr);
  // int idSlit = gmshModelOccAddRectangle(w / 2.0, h / 2.0 - r, 0.0, w, 2.0 * r, -1, 0.0, &ierr);
  int rect[] = {2, idRect};
  // int disk[] = {2, idDisk};
  // int slit[] = {2, idSlit};

  // gmshModelOccCut(rect, 2, disk, 2, NULL, NULL, NULL, NULL, NULL, -1, 1, 1, &ierr);
  // gmshModelOccCut(rect, 2, slit, 2, NULL, NULL, NULL, NULL, NULL, -1, 1, 1, &ierr);
  gmshModelOccSynchronize(&ierr);

  if (theGeometry->elementType == FEM_QUAD) {
    gmshOptionSetNumber("Mesh.SaveAll", 1, &ierr);
    gmshOptionSetNumber("Mesh.RecombineAll", 1, &ierr);
    gmshOptionSetNumber("Mesh.Algorithm", 8, &ierr);
    gmshOptionSetNumber("Mesh.RecombinationAlgorithm", 1.0, &ierr);
    gmshModelGeoMeshSetRecombine(2, 1, 45, &ierr);
    gmshModelMeshGenerate(2, &ierr);
  }

  if (theGeometry->elementType == FEM_TRIANGLE) {
    gmshOptionSetNumber("Mesh.SaveAll", 1, &ierr);
    gmshModelMeshGenerate(2, &ierr);
  }

  return;
}

void geoMeshGenerateGeo(void) {
  femGeo *theGeometry = geoGetGeometry();
  double Lx = 0.25;
  double Ly = 1.0;
  theGeometry->LxPlate = Lx;
  theGeometry->LyPlate = Ly;
  theGeometry->h = Ly * 0.05;
  theGeometry->elementType = FEM_QUAD;

  // geoSetSizeCallback(geoSize);

  /*
  4 ------------------ 3
  |                    |
  |                    |
  5 ------- 6          |
             \         |
              )        |
             /         |
  8 ------- 7          |
  |                    |
  |                    |
  1 ------------------ 2
  */

  int ierr;
  double w = theGeometry->LxPlate;
  double h = theGeometry->LyPlate;
  double r = w / 4;
  double lc = theGeometry->h;
  
  // Geométrie de la plaque

  double big_inner_tri = 4*h/5;
  double thickness_rect = sqrt(pow(h/5, 2) + pow(h/5, 2));
  double hyp_lower_tri = sqrt(pow(thickness_rect, 2) + pow(thickness_rect, 2));

  int idp1 = gmshModelOccAddPoint(0.0, 4*h/5, 0.0, lc, 1, &ierr);
  int idp2 = gmshModelOccAddPoint(0.0, 0.0, 0.0, lc, 2, &ierr);
  int idp3 = gmshModelOccAddPoint(w, 0.0, 0.0, lc, 3, &ierr);
  int idp4 = gmshModelOccAddPoint(w, h, 0.0, lc, 4, &ierr);
  int idp5 = gmshModelOccAddPoint(0.0, h, 0.0, lc, 5, &ierr);

  int idp6 = gmshModelOccAddPoint(-h/5, h, 0.0, lc, 6, &ierr);
  int idp7 = gmshModelOccAddPoint(-(big_inner_tri + hyp_lower_tri), 0.0, 0.0, lc, 7, &ierr);
  int idp8 = gmshModelOccAddPoint(-big_inner_tri, 0.0, 0.0, lc, 8, &ierr);

  int l1 = gmshModelGeoAddLine(idp1, idp2, 1, &ierr);
  int l2 = gmshModelGeoAddLine(idp2, idp3, 2, &ierr);
  int l3 = gmshModelGeoAddLine(idp3, idp4, 3, &ierr);
  int l4 = gmshModelGeoAddLine(idp4, idp5, 4, &ierr);
  int l5 = gmshModelGeoAddLine(idp5, idp6, 5, &ierr);
  int l6 = gmshModelGeoAddLine(idp6, idp7, 6, &ierr);
  int l7 = gmshModelGeoAddLine(idp7, idp8, 7, &ierr);
  int l8 = gmshModelGeoAddLine(idp8, idp1, 8, &ierr);

  //int l4 = gmshModelGeoAddLine(idp4, idp1, 4, &ierr);
  //int idp8 = gmshModelOccAddPoint(0.0, big_inner_tri, 0.0, lc, 8, &ierr);

  // Tags
  //int c1[] = {1};
  //int lTags[] = {l1, l2, l3, l4, l6, l6, l7, l8}; // NB : "-l6" because the curve is reversed

  gmshModelGeoSynchronize(&ierr);

  // int p1 = gmshModelGeoAddPoint(-w / 2, -h / 2, 0., lc, 1, &ierr);
  // int p2 = gmshModelGeoAddPoint(w / 2, -h / 2, 0., lc, 2, &ierr);
  // int p3 = gmshModelGeoAddPoint(w / 2, h / 2, 0., lc, 3, &ierr);
  // int p4 = gmshModelGeoAddPoint(-w / 2, h / 2, 0., lc, 4, &ierr);
  // int p5 = gmshModelGeoAddPoint(-w / 2, r, 0., lc, 5, &ierr);
  // int p6 = gmshModelGeoAddPoint(0., r, 0., lc, 6, &ierr);
  // int p7 = gmshModelGeoAddPoint(0., -r, 0., lc, 7, &ierr);
  // int p8 = gmshModelGeoAddPoint(-w / 2, -r, 0., lc, 8, &ierr);
  // int p9 = gmshModelGeoAddPoint(0., 0., 0., lc, 9, &ierr); // center of circle

  // int l1 = gmshModelGeoAddLine(p1, p2, 1, &ierr);
  // int l2 = gmshModelGeoAddLine(p2, p3, 2, &ierr);
  // int l3 = gmshModelGeoAddLine(p3, p4, 3, &ierr);
  // int l4 = gmshModelGeoAddLine(p4, p5, 4, &ierr);
  // int l5 = gmshModelGeoAddLine(p5, p6, 5, &ierr);
  // int l6 = gmshModelGeoAddCircleArc(p7, p9, p6, 6, 0., 0., 0., &ierr); // NB : the direction of the curve is reversed
  // int l7 = gmshModelGeoAddLine(p7, p8, 7, &ierr);
  // int l8 = gmshModelGeoAddLine(p8, p1, 8, &ierr);

  // int lTags[] = {l1, l2, l3, l4, l5, -l6, l7, l8}; // NB : "-l6" because the curve is reversed
  // int c1[] = {1};
  // c1[0] = gmshModelGeoAddCurveLoop(lTags, 8, 1, 0, &ierr);
  // int s1 = gmshModelGeoAddPlaneSurface(c1, 1, 1, &ierr);

  geoSetSizeCallback(geoSize);
  //gmshModelGeoSynchronize(&ierr);

  if (theGeometry->elementType == FEM_QUAD) {
    gmshOptionSetNumber("Mesh.SaveAll", 1, &ierr);
    gmshOptionSetNumber("Mesh.RecombineAll", 1, &ierr);
    gmshOptionSetNumber("Mesh.Algorithm", 8, &ierr);
    gmshOptionSetNumber("Mesh.RecombinationAlgorithm", 1.0, &ierr);
    gmshModelGeoMeshSetRecombine(2, 1, 45, &ierr);
    gmshModelMeshGenerate(2, &ierr);
  }

  if (theGeometry->elementType == FEM_TRIANGLE) {
    gmshOptionSetNumber("Mesh.SaveAll", 1, &ierr);
    gmshModelMeshGenerate(2, &ierr);
  }
  
  //gmshFltkRun(&ierr);
  printf("ierr = %d\n", ierr);
}

void geoMeshGenerateGeoFile(const char *filename) {
  femGeo *theGeometry = geoGetGeometry();
  int ierr;
  gmshOpen(filename, &ierr);
  ErrorGmsh(ierr);
  if (theGeometry->elementType == FEM_QUAD) {
    gmshOptionSetNumber("Mesh.SaveAll", 1, &ierr);
    gmshOptionSetNumber("Mesh.RecombineAll", 1, &ierr);
    gmshOptionSetNumber("Mesh.Algorithm", 8, &ierr);
    gmshOptionSetNumber("Mesh.RecombinationAlgorithm", 1.0, &ierr);
    gmshModelGeoMeshSetRecombine(2, 1, 45, &ierr);
    gmshModelMeshGenerate(2, &ierr);
  }

  if (theGeometry->elementType == FEM_TRIANGLE) {
    gmshOptionSetNumber("Mesh.SaveAll", 1, &ierr);
    gmshModelMeshGenerate(2, &ierr);
  }
  return;
}

void geoMeshGenerateMshFile(const char *filename) {
  int ierr;
  gmshOpen(filename, &ierr);
  ErrorGmsh(ierr);
  return;
}
