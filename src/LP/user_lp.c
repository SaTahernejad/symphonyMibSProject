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
//Delete
/* system include files */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
using namespace std;
//Delete:end


/* system include files */
//#include <stdio.h>

/* SYMPHONY include files */
#include "sym_constants.h"
#include "sym_macros.h"
#include "sym_lp_u.h"

/* User include files */
#include "user.h"

/*===========================================================================*/

/*===========================================================================*\
 * This file contains the user-written functions for the LP process.
 \*===========================================================================*/

/*===========================================================================*\
 * Here is where the user must receive all of the data sent from
 * user_send_lp_data() and set up data structures. Note that this function is
 * only called if one of COMPILE_IN_LP or COMPILE_IN_TM is FALSE. For 
 * sequential computation, nothing is needed here.
 \*===========================================================================*/

int user_receive_lp_data(void **user)
{
    return(USER_DEFAULT);
}

/*===========================================================================*/

/*===========================================================================*\
 * Here is where the user must create the initial LP relaxation for
 * each search node. Basically, this involves constructing the base matrix in 
 * column ordered format. See the documentation for an explanation of how to 
 * fill out this function.
 \*===========================================================================*/

int user_create_subproblem(void *user, int *indices, MIPdesc *mip, 
			   int *maxn, int *maxm, int *maxnz)
{
    user_problem *prob = (user_problem *) user;

    return(USER_DEFAULT);
}      


/*===========================================================================*/

/*===========================================================================*\
 * This function takes an LP solution and checks it for feasibility. By 
 * default, SYMPHONY checks for integrality. If any integral solution for your 
 * problem is feasible, then nothing needs to be done here.
 \*===========================================================================*/

int user_is_feasible(void *user, double lpetol, int varnum, int *indices,
		     double *values, int *feasible, double *objval,
		     char branching, double *heur_solution)
{
    user_problem *prob = (user_problem *) user;
    if (prob->n==0){
	return(USER_DEFAULT);
    }
    if (prob->n>0){
	*feasible = IP_FEASIBLE;
	int i,j;
	for (i=0; i<varnum; i++){
	    //if(prob->IsInt[i]==1){
	    if ((fabs(values[i]-floor(values[i]))>lpetol) && (fabs(values[i]-ceil(values[i])))>lpetol){
		*feasible = IP_INFEASIBLE;
		return(USER_SUCCESS);
	    }
	    //}
	}
	sym_environment *env2 = sym_open_environment();
	if (!env2){
	    printf("Error initializing environement\n");
	    exit(0);
	}
	double *rhs1 = (double *) malloc(prob->m2 * DSIZE);
	double low_objval2 = 0 ;
	if(prob->type=='G'){
	    double *VarVal = (double *) calloc(prob->numvars, DSIZE);
	    for(i = 0; i < prob->numvars; ++i){
		for(j = 0; j < varnum; ++j){
		    if(indices[j]==i){
			VarVal[i]=values[j];
			break;
		    }
		}
	    }
	    int i1(0);
	    for(i = 0; i < prob->numcons; ++i){
		if(prob->RowInd[i]==1){
		    rhs1[i1] = prob->TotalRhs[i];
		    for(j = prob->upcolumn_starts[i1]; j < prob->upcolumn_starts[i1+1]; ++j){
			rhs1[i1] -= prob->matrix_upvalues[j]*VarVal[prob->matrix_upindices[j]];
		    }
		    i1++;
		}
	    }
	    j=0;
	    for( i= 0; i < prob->numvars; i++){
		if(prob->ColInd[i]==1){
		    low_objval2 = low_objval2+prob->obj[j]*VarVal[i];
		    j++;
		}
	    }
	    FREE(VarVal);
	}
	else{	    
	    double *VarVal = (double *) calloc(2*prob->numvars, DSIZE);
	    for(i = 0; i < 2*prob->numvars; ++i){
		for(j = 0; j < varnum; ++j){
		    if(indices[j]==i){
			VarVal[i]=values[j];
			break;
		    }
		}
	    }

	    for(i= prob->numvars; i < 2*prob->numvars; i++){
		low_objval2 = low_objval2+prob->obj[i-prob->numvars]*VarVal[i];
	    }
	    CoinDisjointCopyN(prob->TotalRhs, prob->numcons,rhs1); 
	    for(i=0; i<prob->numvars; ++i){
		if(VarVal[i]>lpetol){
		    rhs1[i+prob->numcons]=0;
		}
		else{
		    rhs1[i+prob->numcons]=prob->ub[i];
		}
	    }
	    FREE(VarVal);
	}
	/* Load the problem to SYMPHONY */
	sym_explicit_load_problem(env2, prob->n, prob->m2, prob->column_starts, prob->matrix_indices,
				  prob->matrix_values, prob->lb, prob->ub, prob->is_int, prob->obj, 0, prob->sense,
				  rhs1,prob->rngval, true);
	FREE(rhs1);
	sym_set_int_param(env2,"prep_level",-2);
	sym_set_int_param(env2,"verbosity",-2);
	if(varnum==13){
	for(i=0; i<varnum; i++){
	    cout<<indices[i]<<",";
	}
	cout<<endl;
	}
	sym_warm_solve(env2);
        double  objval2=0;
	sym_get_obj_val(env2,&objval2);
	if(fabs(objval2-low_objval2) > lpetol){
	    *feasible = IP_INFEASIBLE;
	}
	//delete [] VarVal;
	return(USER_SUCCESS);
    }//end/
}

/*===========================================================================*/

/*===========================================================================*\
 * Here, the user can specify a special routine for sending back the feasible
 * solution. This need not be used unless there is a special format the user
 * wants the solution in. This function is only called in sequential mode.
 \*===========================================================================*/

int user_send_feasible_solution(void *user, double lpetol, int varnum,
				int *indices, double *values)
{
    return(USER_DEFAULT);
}


/*===========================================================================*/

/*===========================================================================*\
 * This function graphically displays the current fractional solution
 * This is done using the Interactive Graph Drawing program, if it is used.
 \*===========================================================================*/

int user_display_lp_solution(void *user, int which_sol, int varnum,
			     int *indices, double *values)
{
    return(USER_DEFAULT);
}

/*===========================================================================*/

/*===========================================================================*\
 * You can add whatever information you want about a node to help you
 * recreate it. I don't have a use for it, but maybe you will.
 \*===========================================================================*/

int user_add_to_desc(void *user, int *desc_size, char **desc)
{
    return(USER_DEFAULT);
}

/*===========================================================================*/

/*===========================================================================*\
 * Compare cuts to see if they are the same. We use the default, which
 * is just comparing byte by byte.
 \*===========================================================================*/

int user_same_cuts(void *user, cut_data *cut1, cut_data *cut2, int *same_cuts)
{
    return(USER_DEFAULT);
}

/*===========================================================================*/

/*===========================================================================*\
 * This function receives a cut, unpacks it, and adds it to the set of
 * rows to be added to the LP. Only used if cutting planes are generated.
 \*===========================================================================*/

int user_unpack_cuts(void *user, int from, int type, int varnum,
		     var_desc **vars, int cutnum, cut_data **cuts,
		     int *new_row_num, waiting_row ***new_rows)
{
    /* This code is just here as a template for customization. Uncomment to use.*/
#if 0
    int j;
   
    user_problem *prob = (user_problem *) user;

    *new_row_num = 0;
    for (j = 0; j < cutnum; j++){
	switch (cuts[j]->type){
	 
	default:
	    printf("Unrecognized cut type!\n");
	}
    }
#endif
   
    return(USER_DEFAULT);
}

/*===========================================================================*/

/*===========================================================================*\
 * If the user wants to fill in a customized routine for sending and receiving
 * the LP solution, it can be done here. For most cases, the default routines
 * are fine.
 \*===========================================================================*/

int user_send_lp_solution(void *user, int varnum, var_desc **vars, double *x,
			  int where)
{
    return(USER_DEFAULT);
}

/*===========================================================================*/

/*===========================================================================*\
 * This routine does logical fixing of variables
 \*===========================================================================*/

int user_logical_fixing(void *user, int varnum, var_desc **vars, double *x,
			char *status, int *num_fixed)
{
    *num_fixed = 0;

    return(USER_DEFAULT);
}

/*===========================================================================*/

/*===========================================================================*\
 * This function generates the 'next' column. Only used for column generation.
 \*===========================================================================*/

int user_generate_column(void *user, int generate_what, int cutnum,
			 cut_data **cuts, int prevind, int nextind,
			 int *real_nextind, double *colval, int *colind,
			 int *collen, double *obj, double *lb, double *ub)
{
    switch (generate_what){
    case GENERATE_NEXTIND:
	/* Here we just have to generate the specified column. */
	break;
    case GENERATE_REAL_NEXTIND:
	/* In this case, we have to determine what the "real" next edge is*/
	break;
    }

    return(USER_DEFAULT);
}

/*===========================================================================*/

/*===========================================================================*\
 * You might want to print some statistics on the types and quantities
 * of cuts or something like that.
 \*===========================================================================*/

int user_print_stat_on_cuts_added(void *user, int rownum, waiting_row **rows)
{
    return(USER_DEFAULT);
}

/*===========================================================================*/

/*===========================================================================*\
 * You might want to eliminate rows from the local pool based on
 * knowledge of problem structure.
 \*===========================================================================*/

int user_purge_waiting_rows(void *user, int rownum, waiting_row **rows,
			    char *delete_rows)
{
    return(USER_DEFAULT);
}

/*===========================================================================*/

/*===========================================================================*\
 * The user might want to generate cuts in the LP using information
 * about the current tableau, etc. This is for advanced users only.
 \*===========================================================================*/

int user_generate_cuts_in_lp(void *user, LPdata *lp_data, int varnum,
			     var_desc **vars, double *x,
			     int *new_row_num, cut_data ***cuts)
{

    user_problem *prob = (user_problem *) user;
    *cuts = (cut_data **) malloc(sizeof(cut_data *));
    int VariNum=lp_data->n;
    if (prob->n==0){
	return(USER_DEFAULT);
	//return(GENERATE_CGL_CUTS);
    }
   
    if(prob->n>0){
	int i(0);
//Note:Does cut work for ISMP?
	for (i=0; i<varnum; i++){
	    // if(prob->IsInt[i]==1){
	    if ((fabs(x[i]-floor(x[i]))>lp_data->lpetol) && (fabs(x[i]-ceil(x[i])))>lp_data->lpetol){
		return(USER_DEFAULT);
		//return(GENERATE_CGL_CUTS);
	    }
	    // }
	}
	int j(0), index(0);
	int bindingCons(0);
	double slack(0.0);
	double upper(0.0);
	double cutrhs(0.0);
	double upSlack(0.0);
	double downSlack(0.0);
	double value(0.0);  
	int mult(0);
	int cutnum(0);
	const CoinPackedMatrix * matrix = lp_data->si->getMatrixByCol();
	const int numCols = lp_data->si->getNumCols();
	const int numRows = lp_data->si->getNumRows();
	const char * rowsense = lp_data->si->getRowSense();
	const double * rowrange = lp_data->si->getRowRange();
	const double * sol = lp_data->si->getColSolution();
	const double * rhs = lp_data->si->getRightHandSide();
	int *binding  = (int *) malloc((numRows + 2 * numCols) * ISIZE);
	double * slackVal = (double *) malloc( (numRows + 2 * numCols) * DSIZE);
	const double * rowActivity = lp_data->si->getRowActivity();
	const double * collb = lp_data->si->getColLower();
	const double * colub = lp_data->si->getColUpper();
	double * tempVals = (double *) malloc((numCols) * DSIZE);
	int  *indexList = (int *) malloc((numCols) * ISIZE); 
	double * valsList = (double *) malloc((numCols) * DSIZE);
//Binding
	for(i = 0; i < numRows; i++){
	    if((rowsense[i] =='R') && (rowrange[i] < 0))
		slackVal[i] = - rhs[i] + rowActivity[i];
	    else
		slackVal[i] = rhs[i] - rowActivity[i];
	}
	for(i = 0; i < numRows; i++){
	    slack = slackVal[i];
	    switch(rowsense[i]){
	    case 'L':
		if(slack > lp_data->lpetol){
		    binding[i] = 0;
		}
		else if(slack > - lp_data->lpetol){
		    binding[i] = 1;
		    upper += rhs[i];
		    bindingCons++;
		}
		break;
	    case 'G':
		if(slack < - lp_data->lpetol){
		    binding[i] = 0;
		}
		else if(slack < lp_data->lpetol){
		    binding[i] = 1;
		    upper -= rhs[i];
		    bindingCons++;
		}
		break;
	    case 'E':
		if((slack < - lp_data->lpetol) || (slack > lp_data->lpetol)){
		    binding[i] = 0;
		}
		else {
		    binding[i] = 1;
		    //upper += rhs[i];
		    bindingCons++;
		}
		break;
	    }
	}
	for(i = 0; i < numCols; i++){
	    index = numRows + i;
	    upSlack = colub[i]-sol[i];
	    downSlack = sol[i]-collb[i];
	    binding[index] = 0;
	    binding[index+numCols] = 0;
	    if((upSlack > - lp_data->lpetol) && (upSlack < lp_data->lpetol)){
		binding[index] = 1;
		upper += colub[i];
		bindingCons++;
	    }
	    else if((downSlack > - lp_data->lpetol) && (downSlack < lp_data->lpetol)){
		binding[index + numCols] = 1;
		upper -= collb[i];
		bindingCons++;
	    }
	}
	int nonzero(0);
	for(i = 0; i < numCols; i++){
	    tempVals[i]=0;
	    for(j = 0; j < numRows; j++){
		value = matrix->getCoefficient(j, i);
		switch(rowsense[j]){
		case 'L':
		    mult = 1;
		    break;
		case 'G':
		    mult = -1;
		    break;
		case 'E':
		    mult = 0;
		    break;
		}
		tempVals[i] += binding[j] * value * mult;
	    }
// variable upper bound
	    tempVals[i] += binding[numRows + i];
// variable lower bound
	    tempVals[i] += - binding[numRows + numCols + i];
	    if((tempVals[i] > lp_data->lpetol) || (tempVals[i] < - lp_data->lpetol)){
		indexList[nonzero]=i;
		valsList[nonzero]=tempVals[i];
		nonzero++;
	    }
	}
	cutrhs=upper-1;
	cuts[0][0]=create_explicit_cut(nonzero, indexList, valsList, cutrhs, 0, 'L',TRUE);
	cutnum++;
	*new_row_num =1;
	FREE(binding);
	FREE(tempVals);
	FREE(slackVal);
	FREE(indexList);
	FREE(valsList);
	return(USER_SUCCESS);
    }
}

/*===========================================================================*/

/*===========================================================================*\
 * Free all the user data structures
 \*===========================================================================*/

int user_free_lp(void **user)
{
    return(USER_DEFAULT);
}

/*===========================================================================*/

