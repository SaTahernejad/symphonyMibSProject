//*===========================================================================*/
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

/* system include files */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

//sa:start
#include "symphony.h"
#include "sym_master.h"
#include <stdlib.h>
#include <fstream>
using namespace std;
//sa:end

/* SYMPHONY include files */
#include "sym_constants.h"
#include "sym_macros.h"
#include "sym_master_u.h"

/* User include files */
#include "user.h"
#ifdef COMPILE_IN_TM
#ifdef COMPILE_IN_LP
/* fill these in for sequential compilation if needed. */
#ifdef COMPILE_IN_CG
/* fill these in for sequential compilation if needed. */
#endif
#ifdef COMPILE_IN_CP
/* fill these in for sequential compilation if needed. */
#endif
#endif
#endif

/*===========================================================================*\
 * This file contains stubs for the user-written functions for the master 
 * process. The primary function that has to be filled in here is user_io(),
 * where the data for the instance is read in and the user data structure
 * that stores the instance data filled out (this data structure is defined 
 * in user.h). Other than that, the default routines should work fine.
\*===========================================================================*/

/*===========================================================================*/

/*===========================================================================*\
 * This function gives help on command-line switches defined by the user.
 * All user switches have capital letters by convention.
\*===========================================================================*/

void user_usage(void){
  printf("master [ -H ] [ -F file ] \n\n\t%s\n\t%s\n\t%s\n\t%s\n\n",
	 "-H: help (solver-specific switches)",
	 "-F model: model should be read in from file 'model'",
	 "          (MPS format is assumed unless -D is also present)",
	 "-D data: model is in AMPL format and data is in file 'data'");
}

/*===========================================================================*/

/*===========================================================================*\
 * Initialize user-defined data structures. This basically consists of 
 * allocating the memory. If you are using the default data structure,
 * nothing needs to be changed here.
\*===========================================================================*/

int user_initialize(void **user)
{
   /* Create the user's data structure and pass a pointer back to SYMPHONY. */
   user_problem *prob = (user_problem *) calloc(1, sizeof(user_problem));

   *user = prob;

   return(USER_SUCCESS);
}

/*===========================================================================*/

/*===========================================================================*\
 * Parse the user options and read in parameters from the parameter file 
 * given on the command line
\*===========================================================================*/

int user_readparams(void *user, char *filename, int argc, char **argv)
{
   FILE *f;
   char line[50], key[50], value[50], c, tmp;
   int i;
   /* This gives you access to the user data structure*/
   user_problem *prob = (user_problem *) user;
   user_parameters *par = &(prob->par);

   /* Here you can parse the command line for options. By convention, the
      users options should be capital letters */

   for (i = 1; i < argc; i++){
      sscanf(argv[i], "%c %c", &tmp, &c);
      if (tmp != '-')
	 continue;
      switch (c) {
       case 'H':
	 user_usage();
	 exit(0);
	 break;
       case 'F':
	 strncpy(par->infile, argv[++i], MAX_FILE_NAME_LENGTH);
	 break;
	 //sa:start	 
       case 'Z':
	 strncpy(par->infile2, argv[++i], MAX_FILE_NAME_LENGTH);
	 break;
	 //sa:end 
      };
   }
   
   // return(USER_SUCCESS);
    return(USER_DEFAULT);
}

/*===========================================================================*/

/*===========================================================================*\
 * Read in the data file, whose name was given in the parameter file.
 * This file contains instance data. Right now, this function is set up to 
 * read in just the number of columns and number of rows from the file.
 * Add more data as needed to describe the instance and set up the LP
 * relaxation.
\*===========================================================================*/

int user_io(void *user)
{
//sa:start
    user_problem *prob = (user_problem *) user;
    user_parameters *par = &(prob->par);
    char *infile = par->infile;
    char *infile2 = par->infile2;
    FILE *f = NULL;
    CoinMpsIO *mps = new CoinMpsIO;
    if (mps->readMps(infile,"")){
	printf("\nMpp I/O: No problem data file specified\n\n");
	exit(1);
    }
    int i(0),j(0), k(0), m(0);
    prob-> numvars = mps->getNumCols();
    prob-> numcons = mps->getNumRows();
    prob->IsInt  = (int *) malloc(( prob->numvars)* ISIZE);
    for(i=0; i<prob->numvars; ++i){
	if(mps->isInteger(i)){
	    prob->IsInt[i]=1;
	}
	else{
	    prob->IsInt[i]=0;
	}
    }	    
//TotalRhs: Determines all rhs (upper and lower)
    //prob->TotalRhs = new double [prob->numcons];
    prob->TotalRhs = (double *) malloc(( prob->numcons) * DSIZE);
    memcpy(prob->TotalRhs, mps->getRightHandSide(), sizeof(double) * prob->numcons);
    //prob->ColInd = new int[prob->numvars]();
    prob->ColInd  = (int *) calloc(( prob->numvars),ISIZE);
    //prob->RowInd = new int[prob->numcons]();
    prob->RowInd  = (int *) calloc(( prob->numcons),ISIZE);
//newIndex: The array of the new inices of lower cons
//  (without upper cons)
    //int *newIndex = new int[prob->numcons]();
    int *newIndex  = (int *) calloc(( prob->numcons),ISIZE);
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
	    prob->ColInd[iValue]=1;
	    i++;
	}
     else if(key == "LR"){
	 if(!prob->lowerRowInd){
	     //prob->lowerRowInd = new int[prob->m2];
	     prob->lowerRowInd  = (int *) malloc(( prob->m2)*ISIZE);
	 }
	data_stream >> iValue;
	prob->lowerRowInd[j] = iValue;
	prob->RowInd[iValue]=1;
	j++;
     }
     else if(key == "LO"){
       if(!prob->obj)
	   //prob->obj= new double[prob->n2];
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
	/*else if(key == "NZ"){
	data_stream >> dValue;
	prob->nzt= dValue;
	}*/
	/*else if(key == "NZL"){
	data_stream >> dValue;
	prob->nzlow= dValue;
	}*/
    }
    fclose(f);
    prob->nzt = mps->getNumElements() ;
    CoinPackedMatrix * colMatrix = new CoinPackedMatrix();
    *colMatrix = *(mps->getMatrixByCol());
    const double *ColElements = NULL;
    const int *ColIndices = NULL;
    const int *ColLengths = NULL;
    ColElements = colMatrix->getElements();
    ColIndices = colMatrix->getIndices();
    ColLengths = colMatrix->getVectorLengths();
    // ColStarts = colMatrix->getVectorStarts();
    CoinPackedMatrix * rowMatrix = new CoinPackedMatrix();
    *rowMatrix = *(mps->getMatrixByRow());
    const double *RowElements = NULL;
    const int *RowIndices = NULL;
    const int *RowLengths = NULL;
    RowElements = rowMatrix->getElements();
    RowIndices = rowMatrix->getIndices();
    RowLengths = rowMatrix->getVectorLengths();
    // RowStarts = rowMatrix->getVectorStarts();
    //double  *varLB = new double [prob->numvars];
    double *varLB = (double *) malloc( prob->numvars * DSIZE);
    //double  *varUB = new double [prob->numvars];
    double *varUB = (double *) malloc(prob->numvars * DSIZE);
    memcpy(varLB, mps->getColLower(), sizeof(double) * prob->numvars);
    memcpy(varUB, mps->getColUpper(), sizeof(double) * prob->numvars);
//Determine the indexes of lower level rows
    int sum(0);
    for(i = 0; i < prob->numcons; ++i){
	sum = 0;
	if(prob->RowInd[i]==1){
	    for(j = 0; j < i; ++j){
		if(prob->RowInd[j]==0){
		    sum++;
		}
		newIndex[i]=i-sum;
	    }
	}
    }
/*Determine the number of nonzero elements of lower problem
    corresponding to lower variables*/
    k=0;
    prob->nzlow=0;
    for(i = 0; i < prob->numcons; ++i){
	if(prob->RowInd[i]==1){
	    for(j = k; j < k+RowLengths[i]; ++j){
		if(prob->ColInd[RowIndices[j]]==1){
			prob->nzlow++;
		    }
	    }
	}
	    k = k+RowLengths[i];
	}
// set up the inital LP data
    prob-> n = prob->n2;
    prob-> m = prob->m2;
    prob-> nz = prob->nzlow;
// Allocate the arrays
    prob->column_starts  = (int *) malloc(( prob->n + 1) * ISIZE);
    prob->matrix_indices = (int *) malloc(( prob->nz) * ISIZE);
    prob->matrix_values  = (double *) malloc(( prob->nz) * DSIZE);
    prob->lb             = (double *) calloc( prob->n, DSIZE);
    prob->ub             = (double *) malloc( prob->n * DSIZE);
    prob->sense          = (char *) malloc( prob->m * CSIZE);
    prob->rngval         = (double *) calloc( prob->m, DSIZE);
    prob->is_int         = (char *) malloc( prob->n * CSIZE);
    k = 0;
    int i1(0),i2(0),i3(0),i4(0);
    for(i = 0; i < prob->numvars; ++i){
	if(prob->ColInd[i]==1){
	    for(j = k; j < k+ColLengths[i]; ++j){
		if(prob->RowInd[ColIndices[j]]==1){
		    prob->matrix_values[i1] = ColElements[j];
		    prob->matrix_indices[i1] = newIndex[ColIndices[j]];
		    i1++;
		    if(j==k){
			prob->column_starts[i2]=i1-1;
			i2++;
		    }
		}
	    }
	    prob->lb[i3] = varLB[i];
	    prob->ub[i3] = varUB[i];
	    //Note: it may be MIP
	    prob->is_int[i3] = TRUE;
	    i3++;
	}
	k = k+ColLengths[i];
    }
    prob->upcolumn_starts  = (int *) malloc((prob->m+1) * ISIZE);
//Note: Efficient memory allocation
    prob->matrix_upindices = (int *) malloc((prob->nzt-prob->nz+prob->m) * ISIZE);
    prob->matrix_upvalues  = (double *) malloc((prob->nzt-prob->nz+prob->m) * DSIZE);    
    k = 0;
    i1 = 0;
    i2 = 0;
    i3 = 0;
    prob->len_ul=0;
    for(i = 0; i < prob->numcons; ++i){
	if(prob->RowInd[i]==1){
	    i4 = 0;
	    for(j = k; j < k+RowLengths[i]; ++j){
		if(prob->ColInd[RowIndices[j]]==0){
		    i4 = 1;
		    prob->matrix_upvalues[i1] = RowElements[j];
		    prob->matrix_upindices[i1] = RowIndices[j];
		    i1++;
		    prob->len_ul++;
		    if(j==k){
			prob->upcolumn_starts[i2]=i1-1;
	                i2++;
		    }
		}
	    }
	    prob->sense[i3] = 'L';
	    i3++;
	    if(i4==0){
		prob->matrix_upvalues[i1] = 0;
		prob->matrix_upindices[i1] = 0;
		prob->upcolumn_starts[i2]=i1;
		i1++;
		i2++;
	    }
	}
	k = k+RowLengths[i];
    }
/*len_ul: Length of matrix_upvalues*/
    // prob->len_ul = i1;
    prob->column_starts[prob->n] = prob->nz;
    prob->upcolumn_starts[prob->m]=prob->len_ul;
    //sa:end//

    
    //delete [] newIndex;
    FREE(newIndex);
    delete [] ColElements;
    delete [] ColIndices;
    delete [] ColLengths;
    delete [] RowElements;
    delete [] RowIndices;
    delete [] RowLengths;
    //delete [] varLB;
    FREE(varLB);
    //delete [] varUB;
    FREE(varUB);
    return(USER_DEFAULT);
}
   
/*===========================================================================*/

/*===========================================================================*\
 * Here is where the heuristics are performed and an upper bound is calculated.
 * An upper bound can also be specified in the parameter file. This function
 * need not be filled in if no upper bounding is done.
\*===========================================================================*/

int user_start_heurs(void *user, double *ub, double *ub_estimate)
{
   /* This gives you access to the user data structure. */
   user_problem *prob = (user_problem *) user;

   return(USER_DEFAULT);
}

/*===========================================================================*/

/*===========================================================================*\
 * If graph drawing will be used, the user must initialize the drawing
 * window here. This function need not be filled in.
\*===========================================================================*/

int user_init_draw_graph(void *user, int dg_id)
{
   /* This gives you access to the user data structure. */
   user_problem *prob = (user_problem *) user;

   return(USER_DEFAULT);
}

/*===========================================================================*/

/*===========================================================================*\
 * This is the subroutine where the user specifies which variables are to be
 * in the base set and which variables are to be active in the root node but
 * not in the base set (these are called the "extra variables"). This is done
 * by listing the indices of the corresponding variables in arrays named
 * "basevars" and extravars below.
 *
 * The base set of variables form the core that is never removed from the LP
 * relaxation. Extra variables, on the other hand, can be removed if they are
 * fixed by reduced cost or by logic-based rules. Allowing the removal of
 * variables from the relaxation can lead to efficiencies, but there is a
 * price to pay in terms of extra bookkeeping. If possible, it is a good idea
 * to form a base set of variables consisting of those that are "likely" to be
 * present in some optimal solution. If this is not possible, the simplest
 * approach is just to put all the variables in the extra set, which allows
 * them all to be fixed by reduced cost if possible. This is implemented below
 * as an example.
 *
 * Note that each variable must have a unique user index by which the variable
 * can be identified later. Note also that it is possible to have variables
 * that are neither in the base set or active in the root node by using column
 * generation and filling out the function user_generate_column().
\*===========================================================================*/

int user_initialize_root_node(void *user, int *basevarnum, int **basevars,
			      int *basecutnum, int *extravarnum,
			      int **extravars, char *obj_sense,
			      double *obj_offset, char ***colnames,
			      int *colgen_strat)
{
   /* This gives you access to the user data structure. */
   user_problem *prob = (user_problem *) user;

   /* Since we don't know how to form a good set of base variables, we'll put all
      the variables in the extra set */

   /* Set the number of extra variables*/
   *extravarnum = prob->colnum;

#if 0
   /* This code is not really needed because this is the default, so it is 
      commented out and left for illustration. */

   /* Put all the variables in the extra set */
   vars = *extravars = (int *) malloc(varnum * ISIZE);
   for (i = 0; i < varnum; i++){
     vars[i] = i;
   }
#endif
   
   /* Set the number of rows in the initial formulation */
   *basecutnum = prob->rownum;

   /* The set of base variables will be empty */
   *basevarnum = 0;
   *basevars  = NULL;

   return(USER_DEFAULT);
}

/*===========================================================================*/

/*===========================================================================*\
 * Receive the feasible solution. Doesn't need to be filled in.
\*===========================================================================*/

int user_receive_feasible_solution(void *user, int msgtag, double cost,
				   int numvars, int *indices, double *values)
{
   /* This gives you access to the user data structure. */
   user_problem *prob = (user_problem *) user;

   return(USER_DEFAULT);
}

/*===========================================================================*/

/*===========================================================================*\
 * Here, we send the necessary data to the LP process. Notice that
 * there are two cases to deal with. If the LP or the TM are running
 * as separate processes, then we have to send the data by
 * message-passing. Otherwise, we can allocate the user-defined LP data
 * structure here and simply copy the necessary information. This is the
 * only place the user has to sorry about this distinction between
 * configurations. If running sequentially and using the default data
 * structure, nothing needs to be modified in here.
\*===========================================================================*/

int user_send_lp_data(void *user, void **user_lp)
{
   /* This gives you access to the user data structure. */
   user_problem *prob = (user_problem *) user;

#if defined(COMPILE_IN_TM) && defined(COMPILE_IN_LP)
   /* This is is the case when we are copying data directly because the LP is
      not running separately. The easiest thing to do here is just to use the
      same user data structure in both the master and the LP. Then this
      subroutine would simply consist of the line
      
      *user_lp = user;

      Otherwise, this code should be virtually
      identical to that of user_receive_lp_data() in the LP process.*/

   *user_lp = user;
#else
   /* Here, we send that data using message passing and the rest is
      done in user_receive_lp_data() in the LP process */
#endif
   return(USER_SUCCESS);
}

/*===========================================================================*/

/*===========================================================================*\
 * Here, we send the necessary data to the CG process. Notice that
 * there are two cases to deal with. If the CG, LP, or the TM are running
 * as separate processes, then we have to send the data by
 * message-passing. Otherwise, we can allocate the user-defined LP data
 * structure here and simply copy the necessary information. This is the
 * only place the user has to sorry about this distinction between
 * configurations. If running sequentially and using the default data
 * structure, nothing needs to be modified in here.
\*===========================================================================*/

int user_send_cg_data(void *user, void **user_cg)
{
   /* This gives you access to the user data structure. */
   user_problem *prob = (user_problem *) user;

#if defined(COMPILE_IN_TM) && defined(COMPILE_IN_LP) && defined (COMPILE_IN_CG)
   /* This is is the case when we are copying data directly because
      the CG is not running separately. The easiest thing to do here is just
      to use the same user data structure in both the master and the cut
      generator. Then this subroutine would simply consist of 
      
      *user_cg = user;

      Otherwise, this code should be virtually
      identical to that of user_receive_cg_data() in the CG process.*/

   *user_cg = user;
#ifdef CHECK_CUT_VALIDITY
   /* Send the feasible solution here */
#endif
#else
   /* Here, we send that data using message passing and the rest is
      done in user_receive_cg_data() in the CG process */
#ifdef CHECK_CUT_VALIDITY
   /* Send the feasible solution here */
#endif
#endif
   return(USER_SUCCESS);
}

/*===========================================================================*/

/*===========================================================================*\
 * Here, we send the necessary data to the CP process. Notice that
 * there are two cases to deal with. If the CP, LP, or the TM are running
 * as separate processes, then we have to send the data by
 * message-passing. Otherwise, we can allocate the user-defined LP data
 * structure here and simply copy the necessary information. This is the
 * only place the user has to sorry about this distinction between
 * configurations. If running sequentially and using the default data
 * structure, nothing needs to be modified in here.
\*===========================================================================*/

int user_send_cp_data(void *user, void **user_cp)
{
   /* This gives you access to the user data structure. */
   user_problem *prob = (user_problem *) user;

#if defined(COMPILE_IN_TM) && defined(COMPILE_IN_LP) && defined (COMPILE_IN_CP)
   /* This is is the case when we are copying data directly because
      the CP is not running separately. The easiest thing to do here is just
      to use the same user data structure in both the master and the cut
      pool. Then this subroutine would simply consist of 
      
      *user_cp = user;

      Otherwise, this code should be virtually
      identical to that of user_receive_cp_data() in the CP process.*/

   *user_cp = user;
#else
   /* Here, we send that data using message passing and the rest is
      done in user_receive_cp_data() in the CP process */
#endif
   return(USER_SUCCESS);
}

/*===========================================================================*/

/*===========================================================================*\
 * Generally, this function is not needed but you might find some use
 * for it. Someone did :).
\*===========================================================================*/

int user_process_own_messages(void *user, int msgtag)
{
   /* This gives you access to the user data structure. */
   user_problem *prob = (user_problem *) user;

   switch (msgtag){
    case 0:
    default:
      fprintf(stderr, "\nMaster: unknown message type %i!!!\n\n", msgtag);
      exit(1);
   }

   return(USER_DEFAULT);
}

/*===========================================================================*/

/*===========================================================================*\
 * This is the user's chance to display the solution in whatever
 * manner desired. A return value of USER_DEFAULT will cause the
 * default solution display routine to be executed, even if the user displays
 * the solution as well.
\*===========================================================================*/

int user_display_solution(void *user, double lpetol, int varnum, int *indices,
			  double *values, double objval)
{
   return(USER_DEFAULT);
}
   
/*===========================================================================*/

/*===========================================================================*\
 * This is a debugging feature which might
 * allow you to find out why a known feasible solution is being cut off.
\*===========================================================================*/

int user_send_feas_sol(void *user, int *feas_sol_size, int **feas_sol)
{
#ifdef TRACE_PATH

#endif
   return(USER_DEFAULT);
}   

/*===========================================================================*/

/*===========================================================================*\
 * This function frees everything.
\*===========================================================================*/

int user_free_master(void **user)
{
   user_problem *prob = (user_problem *) (*user);

   FREE(prob);

   return(USER_SUCCESS);
}
/*===========================================================================*/

/*===========================================================================*\
 * This function is used to lift the user created cuts during warm starting *
/*===========================================================================*/

int user_ws_update_cuts (void *user, int *size, char **coef, double * rhs, 
			 char *sense, char type, int new_col_num, 
			 int change_type)
{
   return(USER_DEFAULT);
}
/*===========================================================================*/

