/* This is a header file to store all the declarations for the wrapper*/

#ifndef GCTexture_h
#define GCTexture_h

/* Define general functions */
double getPairwiseTerm(int i, int j, int i1, int j1, int labelIJ, int labelIJ1, int*** ImageR,  int*** ImageG,  int*** ImageB, int distVersion);             
int createNodeNumber(int cols, int rows, int **Mask, int **NodeNumber, int NeighSystem);

/* Define the memory allocation utilities */
void destroyMatrixDouble(double **imPr, int mRows, int nCols);
double** buildMatrixDouble (int mRows, int nCols);

double*** buildArrayThreeDouble (int mRows, int nCols, int numIm);
void destroyArrayThreeDouble (double ***imPr, int mRows, int nCols, int numIm);

void destroyMatrixInt(int **imPr, int mRows, int nCols);
int** buildMatrixInt (int mRows, int nCols);

void destroyVectorInt(int *imPr, int mRows);
int* buildVectorInt (int mRows);

int*** buildArrayThreeInt (int mRows, int nCols, int numIm);
void destroyArrayThreeInt (int ***imPr, int mRows, int nCols, int numIm);

double**** buildArrayFourDouble (int mRows, int nCols, int fea, int numIm);
void destroyArrayFourDouble (double ****imPr, int mRows, int nCols, int fea, int
 numIm);

#endif
