/* This is a header file to store all the declarations for the wrapper*/

#ifndef minCut_h
#define minCut_h

/* Define the memory allocation utilities */
void destroyMatrixDouble(double **imPr, int mRows, int nCols);
double** buildMatrixDouble (int mRows, int nCols);

double*** buildArrayThreeDouble (int mRows, int nCols, int numIm);
void destroyArrayThreeDouble (double ***imPr, int mRows, int nCols, int numIm);

void destroyMatrixInt(int **imPr, int mRows, int nCols);
int** buildMatrixInt (int mRows, int nCols);

int*** buildArrayThreeInt (int mRows, int nCols, int numIm);
void destroyArrayThreeInt (int ***imPr, int mRows, int nCols, int numIm);

double**** buildArrayFourDouble (int mRows, int nCols, int fea, int numIm);
void destroyArrayFourDouble (double ****imPr, int mRows, int nCols, int fea, int
 numIm);

#endif
