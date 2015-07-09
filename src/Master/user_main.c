
/*===========================================================================*/
/*                                                                           */
/* This file is part of the SYMPHONY Branch, Cut, and Price Library.         */
/*                                                                           */
/* SYMPHONY was jointly developed by Ted Ralphs (ted@lehigh.edu) and         */
/* Laci Ladanyi (ladanyi@us.ibm.com).                                        */
/*                                                                           */
/* (c) Copyright 2000-2007 Ted Ralphs. All Rights Reserved.                  */
/*                                                                           */
/* This software is licensed under the Eclipse Public License. Please see    */
/* accompanying file for terms.                                              */
/*                                                                           */
/*===========================================================================*/
/*===========================================================================*/

#define CALL_FUNCTION(f)					\
    if ((termcode = f) < 0){					\
	printf("Error detected: termcode = %i\n", termcode);	\
	printf("Exiting...\n\n");				\
	exit(termcode);						\
    }

/*===========================================================================*\
  This file contains the main() for the master process.

  Note that, if you want to use the OSI SYMPHONY interface, you should set the
  USE_OSI_INTERFACE flag and define the COINROOT path in the SYMPHONY 
  Makefile. Otherwise, the C callable library functions will be used by 
  default. See below for the usage.
  \*===========================================================================*/

#if defined(USE_OSI_INTERFACE)

#include "OsiSymSolverInterface.hpp"

int main(int argc, char **argv)
{
    OsiSymSolverInterface si;

    /* Parse the command line */
    si.parseCommandLine(argc, argv);
   
    /* Read in the problem */
    si.loadProblem();

    /* Find a priori problem bounds */
    si.findInitialBounds();

    /* Solve the problem */
    si.branchAndBound();
   
    return(0);
}

#else

#include "symphony.h"
#include "sym_master.h"
#include <stdlib.h>
//sa:start
#include "user.h"
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
using namespace std;
#include "sym_constants.h"
#include "sym_macros.h"
#include "sym_master_u.h"
//sa:end
int main(int argc, char **argv)
{
    /* Create the data structure for storing the problem instance.*/
    user_problem *prob = (user_problem *)calloc(1, sizeof(user_problem));
    cout<<"Please, enter the type of problem: 'G': General,'I': Interdiction"<<endl;
    cin>>prob->type;
    //prob->type='I';
    int termcode;
    sym_environment *env = sym_open_environment();

    sym_version();
    if (!env){
	printf("Error initializing environement\n");
	exit(0);
    }
    CALL_FUNCTION( sym_set_user_data(env, (void *)prob) );
    CALL_FUNCTION( sym_parse_command_line(env, argc, argv) );
    if(prob->type=='G'){
	CALL_FUNCTION( sym_load_problem(env) );
     
	CALL_FUNCTION( sym_find_initial_bounds(env) );
    }
    else{
	//char * infile;
	//char * infile2;
	CALL_FUNCTION( user_read_data(env, prob, prob->par.infile, prob->par.infile2) );
	//CALL_FUNCTION( user_read_data(prob, prob->par.infile) );
    }
    sym_set_int_param(env,"prep_level",-2);
    CALL_FUNCTION( sym_solve(env) );

    CALL_FUNCTION( sym_close_environment(env) );
    return(0);
}

#endif
/*===========================================================================*\
\*===========================================================================*/
//sa:start
int user_read_data(sym_environment *env, user_problem *prob, char *infile, char *infile2)
//int user_read_data(user_problem *prob, char *infile)
{
    FILE *f = NULL;
    CoinMpsIO *mps = new CoinMpsIO;
    if (mps->readMps(infile,"")){
	printf("\nMpp I/O: No problem data file specified\n\n");
	exit(1);
    }

    int i(0),j(0), k(0), y(0);
    double etol = 9.9999999999999995e-08;
    prob-> numvars = mps->getNumCols();
    prob-> numcons = mps->getNumRows();
    prob->IsInt = (int *) malloc(( prob->numvars)* ISIZE);
    prob->IsInt  = (int *) malloc(( prob->numvars)* ISIZE);
    for(i=0; i<prob->numvars; ++i){
	if(mps->isInteger(i)){
	    prob->IsInt[i]=1;
	}
	else{
	    prob->IsInt[i]=0;
	}
    }

    prob->TotalRhs = (double *) malloc(( prob->numcons) * DSIZE);
    memcpy(prob->TotalRhs, mps->getRightHandSide(), sizeof(double) * prob->numcons);
    ifstream data_stream(infile2);
    f = fopen(infile2, "r");
    if ((f = fopen(infile2, "r")) == NULL){
	printf("main(): user file %s can't be opened\n", infile);
	exit(1);
    }
    string key;
    int iValue(0);
    double dValue(0.0);
    i=0;
    prob->lowerColInd=NULL;
    prob->lowerRowInd=NULL;
    prob->obj=NULL;
    while (data_stream >> key){
	if(key == "N"){
	    data_stream >> iValue;
//n2: Number of lower variables
	    prob->n2=iValue;
	}
	else if(key == "M"){
	    data_stream >> iValue;
//m2: Number of lower cons
	    prob->m2=iValue;
	}
	else if(key == "LC"){
	    if(!prob->lowerColInd){
		//prob->lowerColInd=new int[prob->n2];
		prob->lowerColInd  = (int *) malloc(( prob->n2)*ISIZE);
	    }
	    data_stream >> iValue;
	    prob->lowerColInd[i] = iValue;
	    i++;
	}
	else if(key == "LR"){
	    if(!prob->lowerRowInd){
		prob->lowerRowInd  = (int *) malloc(( prob->m2)*ISIZE);
	    }
	    data_stream >> iValue;
	    prob->lowerRowInd[j] = iValue;
	    j++;
	}
	else if(key == "LO"){
	    if(!prob->obj)
		prob->obj  = (double *) malloc( prob->n2 * DSIZE);
	
	    data_stream >> dValue;
	    prob->obj[k] = dValue;
	    k++;
	}
//Note: I considered only Min
	else if(key == "OS"){
	    data_stream >> dValue;
	    prob->lowerObjSense= dValue; //1 min; -1 max
	}
	else if(key == "IC"){
	    if(!prob->intCosts){
		prob->intCosts  = (double *) malloc( prob->n2 * DSIZE);
	    }
       
	    data_stream >> dValue;
	    prob->intCosts[y] = dValue;
	    y++;
	}
	else if(key == "IB"){
	    data_stream >> dValue;
	    prob->InterBudget = dValue;
	}
    }
    fclose(f);
    int numInterdictNZ(0);
    int auxULRows(1);
    int interdictRows(prob->numvars);
    int auxRows(auxULRows + interdictRows);
    int numTotalCols(0), numTotalRows(0);
    for(i = 0; i < prob->numvars; i++){
	if((prob->intCosts[i] > etol) || (prob->intCosts[i] < - etol)){
	    numInterdictNZ++;
	}
    }
    numTotalCols = 2 * prob->numvars;
    numTotalRows = prob->numcons + auxRows;
    double *RelaxvarLB = (double *) malloc( numTotalCols * DSIZE);
    double *RelaxvarUB = (double *) malloc( numTotalCols * DSIZE);
    CoinDisjointCopyN(mps->getColLower(), prob->numvars, RelaxvarLB + prob->numvars);
    CoinDisjointCopyN(mps->getColUpper(), prob->numvars, RelaxvarUB + prob->numvars);       
    CoinFillN(RelaxvarLB, prob->numvars, 0.0); 
    CoinFillN(RelaxvarUB, prob->numvars, 1.0); 
    //CoinFillN(RelaxvarLB + 2 * prob->numvars, 0, 0.0); 
    //CoinFillN(RelaxvarUB + 2 * prob->numvars, 0, 1.0);
    double UpperobjSense = - 1.0;
    double *RelaxObj = (double *) calloc( numTotalCols, DSIZE);
    //CoinZeroN(objCoef, numTotalCols);
    const double *mpsObj =  mps->getObjCoefficients();
    if (UpperobjSense > 0.0) {
	for (j=0; j<prob->numvars; ++j){
	    RelaxObj[j+prob->numvars]=mpsObj[j];
	    prob->obj[j]=-mpsObj[j];
	}
    }
    else {
	for (j=0; j < prob->numvars; ++j) {
	    RelaxObj[j + prob->numvars]=- mpsObj[j];
	    prob->obj[j]=mpsObj[j];
	}
    }
    prob->nzt = mps->getNumElements();
    int Relaxnz (numInterdictNZ+prob->nzt+numTotalCols);
    CoinPackedMatrix * colMatrix = new CoinPackedMatrix();
    *colMatrix = *(mps->getMatrixByCol());
    const double *ColElements = NULL;
    const int *ColIndices = NULL;
    const int *ColLengths = NULL;
    const CoinBigIndex *ColStarts;
    // const ColStarts = colMatrix->getVectorStarts();
    ColElements = colMatrix->getElements();
    ColIndices = colMatrix->getIndices();
    ColLengths = colMatrix->getVectorLengths();
    ColStarts =  colMatrix->getVectorStarts();
    int *Relaxcolumn_starts=(int *) malloc((numTotalCols + 1) * ISIZE);
    int *Relaxmatrix_indices=(int *) malloc(Relaxnz * ISIZE);
    double *Relaxmatrix_values=(double *) malloc(Relaxnz * DSIZE);
    char *RelaxSense=(char *) malloc(numTotalRows * CSIZE);
    double *Relaxrngval=(double *) calloc(numTotalRows, DSIZE);
    char *Relaxis_int=(char *) malloc(numTotalCols * CSIZE);
    double *Relaxrhs=(double *) malloc(numTotalRows * DSIZE);
    Relaxrhs[0]=prob->InterBudget;
    CoinDisjointCopyN(prob->TotalRhs,prob->numcons,Relaxrhs+1);
    k=0;
    for(i=prob->numcons+1; i<numTotalRows; ++i){
	Relaxrhs[i]=RelaxvarUB[k];
	k++;
    }
    CoinFillN(RelaxSense,numTotalRows,'L');
    Relaxcolumn_starts[0]=0;
    Relaxis_int[0]=TRUE;
    for(i=1; i<prob->numvars; ++i){
	if((prob->intCosts[i] > etol) || (prob->intCosts[i] < - etol)){
	    Relaxcolumn_starts[i]=Relaxcolumn_starts[i-1]+2;
	}
	else{
	    Relaxcolumn_starts[i]=Relaxcolumn_starts[i-1]+1;
	}
	Relaxis_int[i]=TRUE;
    }
/*Note: I supposed that all lower variables have
  at least one nonzero coefficient*/
    Relaxcolumn_starts[prob->numvars]=numInterdictNZ+prob->numvars;
    Relaxis_int[prob->numvars]=TRUE;
    for(i=prob->numvars+1; i<numTotalCols; ++i){
	Relaxcolumn_starts[i]=Relaxcolumn_starts[i-1]+
	    ColLengths[i-prob->numvars-1]+1;
	Relaxis_int[i]=TRUE;
    }
    Relaxcolumn_starts[numTotalCols]=Relaxnz;
    j=0;
    for(i=0; i<prob->numvars; ++i){
	if((prob->intCosts[i] > etol) || (prob->intCosts[i] < - etol)){
	    Relaxmatrix_values[j]=prob->intCosts[i];
	    Relaxmatrix_indices[j]=0;
	    j++;
	}
	Relaxmatrix_values[j]=RelaxvarUB[i+prob->numvars];
	Relaxmatrix_indices[j]=i+prob->numcons+1;
	j++;
    }
    for(i=0; i<prob->numvars; ++i){
	for(k=ColStarts[i]; k<ColStarts[i]+ColLengths[i]; ++k){
	    Relaxmatrix_values[j]=ColElements[k];
	    Relaxmatrix_indices[j]=ColIndices[k]+1;
	    j++;
	}
        Relaxmatrix_values[j]=1;
	Relaxmatrix_indices[j]=i+prob->numcons+1;
	j++;
    }
    prob->n=prob->n2;
    int Lownz(prob->nzt+prob->numvars);
    prob->column_starts  = (int *) malloc(( prob->n+1) * ISIZE);
    prob->matrix_indices = (int *) malloc(Lownz * ISIZE);
    prob->matrix_values  = (double *) malloc(Lownz * DSIZE);
    prob->lb             = (double *) calloc(prob->n, DSIZE);
    prob->ub             = (double *) malloc( prob->n * DSIZE);
    prob->sense          = (char *) malloc( prob->m2 * CSIZE);
    prob->rngval         = (double *) calloc( prob->m2, DSIZE);
    prob->is_int         = (char *) malloc( prob->n * CSIZE);
    for(i=0; i<prob->numvars; ++i){
	prob->column_starts[i]=Relaxcolumn_starts[i+prob->numvars]-
	    numInterdictNZ-prob->numvars;
    }
    prob->column_starts[prob->n]=Lownz;
    CoinDisjointCopyN(Relaxmatrix_indices+numInterdictNZ+prob->numvars,Lownz,
		      prob->matrix_indices);
    for(i=0; i<Lownz; ++i){
	prob->matrix_indices[i] -=1;
    }
    CoinDisjointCopyN(Relaxmatrix_values+numInterdictNZ+prob->numvars,Lownz,
		      prob->matrix_values);
    CoinDisjointCopyN(RelaxvarLB+prob->numvars,prob->n,prob->lb);
    CoinDisjointCopyN(RelaxvarUB+prob->numvars,prob->n,prob->ub);
    CoinDisjointCopyN(RelaxSense+1,prob->m2,prob->sense);
    CoinDisjointCopyN(Relaxis_int+prob->numvars,prob->n,prob->is_int);
    sym_explicit_load_problem(env, numTotalCols, numTotalRows, Relaxcolumn_starts,
			      Relaxmatrix_indices, Relaxmatrix_values, RelaxvarLB,
			      RelaxvarUB, Relaxis_int, RelaxObj, 0, RelaxSense,
			      Relaxrhs, Relaxrngval, true);
    FREE(RelaxvarLB);
    FREE(RelaxvarUB);
    FREE(RelaxObj);
    FREE(Relaxcolumn_starts);
    FREE(Relaxmatrix_indices);
    FREE(Relaxmatrix_values);
    FREE(RelaxSense);
    FREE(Relaxrngval);
    FREE(Relaxis_int);
    FREE(Relaxrhs);
    return (FUNCTION_TERMINATED_NORMALLY);
}
