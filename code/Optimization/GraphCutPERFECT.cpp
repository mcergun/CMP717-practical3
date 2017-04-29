//  This function implements a wrapper to compute the max-flow (or
//    min-cut) for binary label Graph
//    Mexify it with:  mex GraphCut.cpp graph.cpp maxflow.cpp


#include "mex.h"
#include "GraphCut.h"
#include "graph.h"
#include "energy.h"

#include <math.h>
#include <stdlib.h>

#define INFINITY 100000

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
  
  /* check for the proper no. of input and outputs */
  if (nrhs != 7)
    mexErrMsgTxt("7 input arguments are required");
  if (nlhs>1)
    mexErrMsgTxt("Too many outputs");
  
  /* Get the no. of dimentsions and size of each dimension in the input 
     array */
  int numDimsUnary = mxGetNumberOfDimensions(prhs[0]);
  const int *sizeDimUnary = mxGetDimensions(prhs[0]);
  
  if (numDimsUnary!=3)
  {
      mexErrMsgTxt("Check the input arrays");
  }
  
  /* Get the dimensions of the input image */
  int rows = sizeDimUnary[0]; 
  int cols = sizeDimUnary[1]; 
  int labelNum = 2;
  int Output = 1;
  
  /* Get the pointers to the input Data */
  double *UnaryPtr = mxGetPr(prhs[0]);
  double *lambda1Ptr = mxGetPr(prhs[1]);
  double *lambda2Ptr = mxGetPr(prhs[2]);
  double *regularisationPtr = mxGetPr(prhs[3]);
  double *NeighSystemPtr = mxGetPr(prhs[4]);
  double *ImagePtr = mxGetPr(prhs[5]);
  double *EdgeBetaPtr = mxGetPr(prhs[6]);
  
  double lambda1 = (double) *lambda1Ptr;
  double lambda2 = (double) *lambda2Ptr;
  double regular = (double) *regularisationPtr;
  int NeighSystem = (int) *NeighSystemPtr;
  double EdgeBeta = (double) *EdgeBetaPtr;
  
  if (Output)
    printf("Rows:%d, Cols:%d, lambda1: %e, lambda2: %e, regular: %e, NeighSystem 4(0),8(1): %d\n", rows,cols,lambda1,lambda2,regular,NeighSystem); 
  
  /* create matrices */
  int ***Image, **NodeNumber, **Label;
  double ***Unary;
  
  Image = buildArrayThreeInt(rows, cols, 3);
  Unary = buildArrayThreeDouble(rows, cols, labelNum);
  NodeNumber = buildMatrixInt(rows, cols);
  Label = buildMatrixInt(rows, cols);
  
  /* Assign the data from the */
  for (int k=0; k<3; k++){
       for (int j=0; j<cols; j++){
            for (int i=0; i<rows; i++){
                Image[i][j][k] = (int) (*ImagePtr * 255.0); // convert into 0-255
                // printf("Rows:%d, Cols:%d, Image %d, ImageR %d, ImageG %d, ImageB %d, Unary %d \n", i,j,k, ImageR[i][j][k], ImageG[i][j][k], ImageB[i][j][k], Unary[i][j][k]);
                ImagePtr++;                
            }
       }
  }
  
  /* Assign the data from the */
   for (int k=0; k<labelNum; k++){
       for (int j=0; j<cols; j++){
            for (int i=0; i<rows; i++){

            Unary[i][j][k] = (double) (*UnaryPtr); 
            
            // if (i==10 && j==10) 
            //     printf("Unary %e",Unary[i][j][k]);
            
            UnaryPtr++;
        }
    }
  }
    
  for (int j=0; j<cols; j++){
    for (int i=0; i<rows; i++){
       Label[i][j] = 0; 
       NodeNumber[i][j] = -1; // no number assigned

    }
  }
  
  /* Define the size of the output array */
  int *sizeLabelOut = new int[2];
  sizeLabelOut[0] = rows;
  sizeLabelOut[1] = cols;

  /* Create the outGoing Array  */
  plhs[0] = mxCreateNumericArray(2, sizeLabelOut, mxDOUBLE_CLASS, mxREAL);   
  double *labelOutPtr = mxGetPr(plhs[0]); 

  /*******************************************/
  /******** Define the Graph ******/
  /*******************************************/
    
  // extra parameters
  int nodeCount, maxNodes, neighR, neighC, neighCountNode;
  
  /********* Create Node Number from Mask **************/
  maxNodes = createNodeNumber(cols,rows,NodeNumber);
  
   
  /****** Create the Energy *************************/
  Energy::Var *varx = new Energy::Var [maxNodes];			
  Energy *e = new Energy();

  /****** Build Unary Term *************************/
  for (int j=0; j<cols; j++)
  {
      for (int i=0; i<rows; i++)
      {	
          nodeCount = NodeNumber[i][j];
          
          // add node
          varx[nodeCount] =  e -> add_variable();
          
          // add likelihood
          e -> add_term1(varx[nodeCount], (double)Unary[i][j][0], (double)Unary[i][j][1]);
          
          // add regularisation
          // e -> add_term1(varx[nodeCount], regular, 0.0);
      }
  }
  
  /****** Pairwise Terms of the Energy *************************/
  double term01;
  for (int j=0; j<cols; j++) // -1 not to step over the image
  {
       for (int i=0; i<rows; i++) // -1 not to step over the image
       {
            nodeCount = NodeNumber[i][j];
               
            // horizontal and vertical
            for(int whichNeigh=2; whichNeigh<=3; whichNeigh++) //  loop over the two neighbors          
            {
				if (whichNeigh==2) // 2 - bottom neighbor
                {
					neighR = i+1; 
					neighC = j; 
				}
				else // 3 - right neighbor
                {  
					neighR = i; 
					neighC = j+1;
				}
                
                if ( (neighR<rows) & (neighC<cols) )
                {
                    neighCountNode = NodeNumber[neighR][neighC];

                    // initialise
                    term01 = 0.0;

                    // get the Constant term
                    if (lambda1>0.0)
                    {
                        term01 = term01+lambda1;
                    }       

                    // get the Edge term
                    if (lambda2>0.0)
                    {
                        term01 = term01+lambda2*getEdgeTerm(i, j, neighR, neighC, Image, EdgeBeta);
                    }

                    e -> add_term2(varx[nodeCount],varx[neighCountNode],0.0,term01,term01,0.0);
                }   
            }
            
            // diagonal part
                if (NeighSystem)
                {
                    for(int whichNeigh=2; whichNeigh<=3; whichNeigh++) //  loop over the two neighbors          
                    { 
                        if (whichNeigh==2) // 2 - bottom, left neighbor
                        {
                            neighR = i+1; 
                            neighC = j-1; 
                        }
                        else // 3 - bottom, right neighbor
                        {  
                            neighR = i+1; 
                            neighC = j+1;
                        }
					
                        if ( (neighR<rows) & (neighC<cols) &  (neighC>-1) )
                        {
                            
                            neighCountNode = NodeNumber[neighR][neighC];
                        
                            // initialise
                            term01 = 0.0;

                            // get the Constant term
                            if (lambda1>0.0)
                            {
                                term01 = term01+lambda1;
                            }       

                            // get the Edge term
                            if (lambda2>0.0)
                            {
                                term01 = term01+lambda2*getEdgeTerm(i, j, neighR, neighC, Image, EdgeBeta);
                            }
                                
                            // add the pairwise term
                            e -> add_term2(varx[nodeCount],varx[neighCountNode],0.0,term01/sqrt(2.0),term01/sqrt(2.0),0.0);
                            
                        }
                    } 
                }
       }
  } 
                    
  /****** Perform graph cuts *************************/
  Energy::TotalValue Emin = e -> minimize();
                   
   // get the values of the nodes 
   for (int j=0; j<cols; j++)
   {
		for (int i=0; i<rows; i++)
        {			
             nodeCount = NodeNumber[i][j];
             if(e->get_var(varx[nodeCount]))
                Label[i][j] = 1;
             else
                Label[i][j] = 0;
        }
   }

   if (Output)
      printf("Energy: %f \n ", Emin);  
        
   // Delete the energy
   delete e;
   delete [] varx;	
  

    /*************** Assign the output array *************/
	
	for (int j=0; j<cols; j++){
        for (int i=0; i<rows; i++){
            *labelOutPtr =  (double) Label[i][j]; //convert into matlab indices
            //      mexPrintf("assoc1[%d][%d]: %f\n",i, j, assoc1[i][j]);
            labelOutPtr++;
        }
    }
    
    /************ GRAPH CUT CODE ENDS *****************/
   
    destroyArrayThreeInt(Image, rows, cols, 3);
    destroyArrayThreeDouble(Unary, rows, cols, labelNum);
    destroyMatrixInt(NodeNumber, rows, cols);
    destroyMatrixInt(Label, rows, cols);

} // mex function 
  
// General Functions 
double getEdgeTerm(int i, int j, int i1, int j1, int*** Image, double EdgeBeta)
{
        // potts model 
        double dist1 = pow( (((double) (Image[i][j][0]-Image[i1][j1][0]))/255.0), 2);
        dist1 += pow( (((double) (Image[i][j][1]-Image[i1][j1][1]))/255.0), 2);
        dist1 += pow( (((double) (Image[i][j][2]-Image[i1][j1][2]))/255.0), 2);
         
        // printf("Number of al % e \n", (1.0/ (1.0+sqrt(dist1))));
        return exp(-EdgeBeta*dist1);// (1.0/ (1.0+sqrt(dist1)));                      
}

int createNodeNumber(int cols, int rows, int **NodeNumber)
{
    int nodeCount=0;
    for (int j=0; j<cols; j++)
    {
		for (int i=0; i<rows; i++)
        {	
            NodeNumber[i][j] = nodeCount;
            nodeCount++;
        }
    }
    
    return  nodeCount;
}

  /****** UTILITY FUNCTIONS *****************/

int* buildVectorInt (int mRows){
  
  int *imPr = new int [mRows];
 
  return imPr;
}

void destroyVectorInt(int *imPr, int mRows){
  
    delete [] imPr ;
}

double** buildMatrixDouble (int mRows, int nCols){
  
  double **imPr = new double* [mRows];
  
  for (int i=0; i<mRows; i++){
    imPr[i] = new double[nCols] ;
  }
  return imPr;
}

void destroyMatrixDouble(double **imPr, int mRows, int nCols){
  
    for (int i=0; i<mRows; i++){
      delete [] imPr[i];
    }
    
    delete [] imPr ;
}
int** buildMatrixInt (int mRows, int nCols){
  
  int **imPr = new int* [mRows];
  
  for (int i=0; i<mRows; i++){
    imPr[i] = new int[nCols] ;
  }
  return imPr;
}

void destroyMatrixInt(int **imPr, int mRows, int nCols){
  
    for (int i=0; i<mRows; i++){
      delete [] imPr[i];
    }
    
    delete [] imPr ;
}

int*** buildArrayThreeInt (int mRows, int nCols, int numIm){

  int ***imPr = new int** [mRows];

  for (int i=0; i<mRows; i++){
    imPr[i] = new int* [nCols] ;
    for (int j=0; j<nCols; j++)
      imPr[i][j] = new int[numIm];
  }
  return imPr;
}

void destroyArrayThreeInt (int ***imPr, int mRows, int nCols, int numIm){


  for (int i=0; i<mRows; i++){
    for (int j=0; j<nCols; j++){
      delete [] imPr[i][j] ;
    }
    delete [] imPr[i];
  }
  delete [] imPr;

}

double*** buildArrayThreeDouble (int mRows, int nCols, int numIm){

  double ***imPr = new double** [mRows];

  for (int i=0; i<mRows; i++){
    imPr[i] = new double* [nCols] ;
    for (int j=0; j<nCols; j++)
      imPr[i][j] = new double[numIm];
  }
  return imPr;
}

void destroyArrayThreeDouble (double ***imPr, int mRows, int nCols, int numIm){

  for (int i=0; i<mRows; i++){
    for (int j=0; j<nCols; j++){
      delete [] imPr[i][j] ;
    }
    delete [] imPr[i];
  }
  delete [] imPr;

}