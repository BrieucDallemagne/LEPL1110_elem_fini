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
  double Lx = 0.2;
  double Ly = 1.0;
  theGeometry->LxPlate = Lx;
  theGeometry->LyPlate = Ly;
  theGeometry->h = Ly * 0.05;
  theGeometry->elementType = FEM_QUAD;

  geoSetSizeCallback(geoSize);

  int ierr;
  double w = theGeometry->LxPlate;
  double h = theGeometry->LyPlate;
  double r = w / 4;
  double lc = theGeometry->h;
 
  // Points

  double big_inner_tri = 4*h/5;
  double thickness_rect = sqrt(pow(h/10, 2) + pow(h/5, 2));
  double hyp_lower_tri = sqrt(pow(thickness_rect, 2) + pow(thickness_rect, 2));

  int idp1 = gmshModelGeoAddPoint(0., 4*h/5, 0., lc, 1, &ierr);
  int idp2 = gmshModelGeoAddPoint(0., 0., 0., lc, 2, &ierr);
  int idp3 = gmshModelGeoAddPoint(w, 0., 0., lc, 3, &ierr);

  int idp4 = gmshModelGeoAddPoint(w, 2*h/5, 0., lc, 4, &ierr);
  int idp5 = gmshModelGeoAddPoint(w, 4*h/5, 0., lc, 5, &ierr);
  int idp6 = gmshModelGeoAddPoint(3*w, 2*h/5, 0., lc, 6, &ierr);

  int idp7 = gmshModelGeoAddPoint(5*w, 0., 0., lc, 7, &ierr);
  int idp8 = gmshModelGeoAddPoint(6*w, 0., 0., lc, 8, &ierr);
  int idp9 = gmshModelGeoAddPoint(6*w, h/5, 0., lc, 9, &ierr);

  int idp10 = gmshModelGeoAddPoint(5*w, h/5, 0., lc, 10, &ierr);
  int idp11 = gmshModelGeoAddPoint(5*w, 2*h/5, 0., lc, 11, &ierr);
  int idp12 = gmshModelGeoAddPoint(4*w, 2*h/5, 0., lc, 12, &ierr);

  int idp13 = gmshModelGeoAddPoint(w, h, 0., lc, 13, &ierr);
  int idp14 = gmshModelGeoAddPoint(-h/10, h, 0., lc, 14, &ierr);
  int idp15 = gmshModelGeoAddPoint(-(big_inner_tri + hyp_lower_tri), 0., 0., lc, 15, &ierr);
  int idp16 = gmshModelGeoAddPoint(-big_inner_tri, 0., 0., lc, 16, &ierr);

  gmshModelGeoSynchronize(&ierr);

  // Lines and circle arcs

  int l1 = gmshModelGeoAddLine(idp1, idp2, 1, &ierr);
  int l2 = gmshModelGeoAddLine(idp2, idp3, 2, &ierr);
  int l3 = gmshModelGeoAddLine(idp3, idp5, 3, &ierr);

  int arc1_bottom = gmshModelGeoAddCircleArc(idp5, idp4, idp6, 4, 0., 0., 0., &ierr);
  int arc2_bottom = gmshModelGeoAddCircleArc(idp6, idp11, idp7, 5, 0., 0., 0., &ierr);

  int l4 = gmshModelGeoAddLine(idp7, idp8, 6, &ierr);
  int l5 = gmshModelGeoAddLine(idp8, idp9, 7, &ierr);
  int l6 = gmshModelGeoAddLine(idp9, idp10, 8, &ierr);

  int arc3_up = gmshModelGeoAddCircleArc(idp10, idp11, idp12, 9, 0., 0., 0., &ierr);
  int arc4_up = gmshModelGeoAddCircleArc(idp12, idp4, idp13, 10, 0., 0., 0., &ierr);

  int l7 = gmshModelGeoAddLine(idp13, idp14, 11, &ierr);
  int l8 = gmshModelGeoAddLine(idp14, idp15, 12, &ierr);
  int l9 = gmshModelGeoAddLine(idp15, idp16, 13, &ierr);
  int l10 = gmshModelGeoAddLine(idp16, idp1, 14, &ierr);

  gmshModelGeoSynchronize(&ierr);

  // Tags

  int lTags[] = {l1, l2, l3, arc1_bottom, arc2_bottom, l4, l5, l6, arc3_up, arc4_up, l7, l8, l9, l10}; // NB : "-l6" because the curve is reversed
  int c1[] = {1};
  c1[0] = gmshModelGeoAddCurveLoop(lTags, 14, 1, 1, &ierr);
  int s1 = gmshModelGeoAddPlaneSurface(c1, 1, 1, &ierr);
  gmshModelGeoSynchronize(&ierr);

  if (theGeometry->elementType == FEM_QUAD) {
    gmshOptionSetNumber("Mesh.SaveAll", 1, &ierr);
    gmshOptionSetNumber("Mesh.RecombineAll", 1, &ierr);
    gmshOptionSetNumber("Mesh.Algorithm", 5, &ierr);
    gmshOptionSetNumber("Mesh.RecombinationAlgorithm", 1.0, &ierr);
    gmshModelGeoMeshSetRecombine(2, 1, 45, &ierr);
    gmshModelMeshGenerate(2, &ierr);
  }

  if (theGeometry->elementType == FEM_TRIANGLE) {
    gmshOptionSetNumber("Mesh.SaveAll", 1, &ierr);
    gmshModelMeshGenerate(2, &ierr);
  }

  // gmshFltkRun(&ierr);

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
