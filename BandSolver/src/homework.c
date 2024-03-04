

#include"fem.h"


#ifndef NORENUMBER 
double* theGlobalArrayPos;

int diffnode( const void * first , const void * second ) {
    int *One = (int*)first;
    int *Two = (int*)second;
    double diff = theGlobalArrayPos[*Two] - theGlobalArrayPos[*One];
    return diff;
}

void femMeshRenumber(femMesh *theMesh, femRenumType renumType)
{

    int i;
    int *inverse = (int*)malloc(sizeof(int)*theMesh->nodes->nNodes );

    for (i = 0; i < theMesh->nodes->nNodes; i++)
        inverse[i] = i;

    switch(renumType){
        case FEM_NO:
            break;
        case FEM_XNUM:
            theGlobalArrayPos = theMesh->nodes->X;
            qsort(inverse, theMesh->nodes->nNodes, sizeof(int), diffnode);
            break;
        case FEM_YNUM:
            theGlobalArrayPos = theMesh->nodes->Y;
            qsort(inverse, theMesh->nodes->nNodes, sizeof(int), diffnode);
            break;
    default:
        Error(" Unexpected renumbering option ");
    }

    for (i = 0; i < theMesh->nodes->nNodes; i++)
        theMesh->nodes->number[inverse[i]] = i;
    free(inverse);
        

}


#endif
#ifndef NOBAND 

int MAX(int a,int b){return ((a) > (b) ? (a) : (b));}
int MIN(int a,int b){return ((a) < (b) ? (a) : (b));}

int femMeshComputeBand(femMesh *theMesh)
{
    int nLocal = theMesh->nLocalNode;
    int lst[nLocal];
    int Band = 0;    
    int Max;
    int Min;

    for (int i = 0; i < theMesh -> nElem; i++)
    {
        for (int j = 0; j < nLocal; ++j){
            lst[j] = theMesh ->nodes-> number[theMesh -> elem[i * nLocal + j]];
            if(j == 0){
                Min = lst[0];
                Max = lst[0];
            }
            else{
                Max = MAX(lst[j], Max);
                Min = MIN(lst[j], Min);
            }  
        }
        Band = MAX(Max - Min, Band);
    }

    return (Band+1);
}


#endif
#ifndef NOBANDASSEMBLE


void femBandSystemAssemble(femBandSystem* myBandSystem, double *Aloc, double *Bloc, int *lst, int nLoc)
{

    for (int i = 0; i < nLoc; i++){
        int Row = lst[i];
        myBandSystem->B[Row] += Bloc[i];
        for (int j = 0; j < nLoc; j++){
            int Col = lst[j];
            if (Col >= Row){
                myBandSystem->A[Row][Col] += Aloc[i*nLoc+j];
            }
        }
    }
}

#endif
#ifndef NOBANDELIMINATE


double  *femBandSystemEliminate(femBandSystem *Band)
{
    double **A, *B, factor;
    int i, j, k, jend;
    int size = Band -> size;
    int band = Band -> band;
    A = Band -> A;
    B = Band -> B;


    for (k = 0; k < size; k++){
        if(k+band >= size){
            jend = size;
        }else{
            jend = k+band;
        }
        for (i = k + 1; i < jend; i++){
            factor = A[k][i] / A[k][k];
            for (j = i; j < jend; j++){
                A[i][j] -= A[k][j] * factor;
            }
            B[i] -= B[k] * factor;
        }
    }

    for (i = (size - 1); i >= 0; i--){
        factor = 0;
        if(i+band < size){
            jend = i+band;
        }else{
            jend = size;
        }
        for (j = i + 1; j < jend; j++){
            factor += A[i][j] * B[j];
        }
        B[i] = (B[i] - factor) / A[i][i];
    }

    return (Band -> B);
}


#endif



