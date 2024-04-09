/*
 *  main.c
 *  Projet 2023-2024
 *  Elasticite lineaire plane
 *
 *  Preprocesseur
 *
 *  Copyright (C) 2023 UCL-IMMC : Vincent Legat, Miguel De Le Court
 *  All rights reserved.
 *
 */

#include "glfem.h"

int main(void) {

  //
  //  -1- Construction de la geometrie
  //
  geoInitialize();
  femGeo *theGeometry = geoGetGeometry();

  // OPTION 1 : Utilisation de GMSH avec OpenCascade
  // theGeometry->h = 0.05;
  // geoMeshGenerate();

  // OPTION 2 : Utilisation de GMSH directement
  theGeometry->h = 0.05;
  geoMeshGenerateGeo();

  // OPTION 3 : Lecture d'un fichier .geo
  // theGeometry->h = 0.05;
  // geoMeshGenerateGeoFile("../data/mesh.geo");

  // OPTION 4 : Lecture d'un fichier .msh
  // geoMeshGenerateMshFile("../data/mesh.msh");


  geoMeshImport();
  geoSetDomainName(0, "Rectangle left");
  geoSetDomainName(1, "Rectangle bottom");
  geoSetDomainName(2, "Rectangle right");
  geoSetDomainName(3, "Arc 1");
  geoSetDomainName(4, "Arc 2");
  geoSetDomainName(5, "Small rectangle bottom"); 
  geoSetDomainName(6, "Small rectangle right");
  geoSetDomainName(7, "Small rectangle top");
  geoSetDomainName(8, "Arc3");
  geoSetDomainName(9, "Arc4");
  geoSetDomainName(10, "Rectangle top");
  geoSetDomainName(11, "Leaning rectangle left");
  geoSetDomainName(12, "Leaning rectangle bottom");
  geoSetDomainName(13, "Leaning rectangle right");


  geoMeshWrite("../data/mesh.txt");
  geoMeshWrite("/element_fini/LEPL1110_elem_fini/groupe048-bdallemagne-tirgel/Project/data/mesh.txt");
  geoMeshWrite("/element_fini/LEPL1110_elem_fini/groupe048-bdallemagne-tirgel/ProjectPostProcessor/data/mesh.txt");
  

  //
  //  -2- Definition du probleme
  //

  double E = 211.e9;
  double nu = 0.3;
  double rho = 7.85e3;
  double gx = 0;
  double gy = -9.81;

  femProblem *theProblem = femElasticityCreate(theGeometry, E, nu, rho, gx, gy, PLANAR_STRAIN);

  femElasticityAddBoundaryCondition(theProblem, "Rectangle bottom", DIRICHLET_XY, 0.0, 0.0);
  femElasticityAddBoundaryCondition(theProblem, "Rectangle left", DIRICHLET_X, 0.0, NAN);
  // femElasticityAddBoundaryCondition(theProblem, "Rectangle right", DIRICHLET_X, 0.0, NAN);
  femElasticityAddBoundaryCondition(theProblem, "Leaning rectangle bottom", DIRICHLET_XY, 0.0, 0.0);
  // femElasticityAddBoundaryCondition(theProblem, "Leaning rectangle left", DIRICHLET_XY, 0.0, 0.0);
  // femElasticityAddBoundaryCondition(theProblem, "Leaning rectangle right", DIRICHLET_XY, 0.0, 0.0);
  femElasticityAddBoundaryCondition(theProblem, "Small rectangle bottom", DIRICHLET_XY, 0.0, 0.0);
  // femElasticityAddBoundaryCondition(theProblem, "Arc 1", DIRICHLET_XY, 0.0, 0.0);
  // femElasticityAddBoundaryCondition(theProblem, "Arc 2", DIRICHLET_XY, 0.0, 0.0);
  // femElasticityAddBoundaryCondition(theProblem, "Small rectangle top", DIRICHLET_XY, 0.0, 0.0);
  // femElasticityAddBoundaryCondition(theProblem, "Small rectangle right", DIRICHLET_XY, 0.0, 0.0);
  // femElasticityAddBoundaryCondition(theProblem, "Rectangle top", DIRICHLET_XY, 0.0, 0.0);
  femElasticityAddBoundaryCondition(theProblem, "Arc3", NEUMANN_Y, -1000*gy, NAN);
  
  femElasticityPrint(theProblem);
  femElasticityWrite(theProblem, "../data/problem.txt");
  femElasticityWrite(theProblem, "../../Project/data/problem.txt");
  femElasticityWrite(theProblem, "../../ProjectPostProcessor/data/problem.txt");

  //
  //  -3- Champ de la taille de référence du maillage (uniquement pour la visualisation)
  //

  double *meshSizeField = malloc(theGeometry->theNodes->nNodes * sizeof(double));
  femNodes *theNodes = theGeometry->theNodes;
  for (int i = 0; i < theNodes->nNodes; ++i)
    meshSizeField[i] = theGeometry->geoSize(theNodes->X[i], theNodes->Y[i]);
  double hMin = femMin(meshSizeField, theNodes->nNodes);
  double hMax = femMax(meshSizeField, theNodes->nNodes);
  printf(" ==== Global requested h : %14.7e \n", theGeometry->h);
  printf(" ==== Minimum h          : %14.7e \n", hMin);
  printf(" ==== Maximum h          : %14.7e \n", hMax);

  //
  //  -4- Visualisation
  //

  int mode = 1;
  int domain = 0;
  int freezingButton = FALSE;
  double t, told = 0;
  char theMessage[MAXNAME];

  GLFWwindow *window = glfemInit("EPL1110 : Project 2023-24 ");
  glfwMakeContextCurrent(window);
  glfwSetScrollCallback(window, scroll_callback);
  glfwSetMouseButtonCallback(window, mouse_button_callback);

  do {
    int w, h;
    glfwGetFramebufferSize(window, &w, &h);
    glfemReshapeWindows(window, theGeometry->theNodes, w, h);

    t = glfwGetTime();
    if (glfwGetKey(window, 'D') == GLFW_PRESS) {
      mode = 0;
    }
    if (glfwGetKey(window, 'V') == GLFW_PRESS) {
      mode = 1;
    }
    if (glfwGetKey(window, 'N') == GLFW_PRESS && freezingButton == FALSE) {
      domain++;
      freezingButton = TRUE;
      told = t;
    }

    if (t - told > 0.5) {
      freezingButton = FALSE;
    }
    if (mode == 1) {
      glfemPlotField(theGeometry->theElements, meshSizeField);
      glfemPlotMesh(theGeometry->theElements);
      sprintf(theMessage, "Number of elements : %d ", theGeometry->theElements->nElem);
      glColor3f(1.0, 0.0, 0.0);
      glfemMessage(theMessage);
    }
    if (mode == 0) {
      domain = domain % theGeometry->nDomains;
      glfemPlotDomain(theGeometry->theDomains[domain]);
      sprintf(theMessage, "%s : %d ", theGeometry->theDomains[domain]->name, domain);
      glColor3f(1.0, 0.0, 0.0);
      glfemMessage(theMessage);
    }

    glfwSwapBuffers(window);
    glfwPollEvents();
  } while (glfwGetKey(window, GLFW_KEY_ESCAPE) != GLFW_PRESS && glfwWindowShouldClose(window) != 1);

  // Check if the ESC key was pressed or the window was closed

  free(meshSizeField);
  femElasticityFree(theProblem);
  geoFree();
  glfwTerminate();

  exit(EXIT_SUCCESS);
  return 0;
}
