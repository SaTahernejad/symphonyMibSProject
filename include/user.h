/*===========================================================================*/
/*                                                                           */
/* This file is part of the SYMPHONY Branch, Cut, and Price Library.         */
/*                                                                           */
/* SYMPHONY was jointly developed by Ted Ralphs (ted@lehigh.edu) and         */
/* Laci Ladanyi (ladanyi@us.ibm.com).                                        */
/*                                                                           */
/* (c) Copyright 2000-2013 Ted Ralphs. All Rights Reserved.                  */
/*                                                                           */
/* This software is licensed under the Eclipse Public License. Please see    */
/* accompanying file for terms.                                              */
/*                                                                           */
/*===========================================================================*/

#ifndef _USER_H
#define _USER_H
/*sa: start*/
#include "sym_master.h"
#include "sym_macros.h"
//sa:end


/*---------------------------------------------------------------------------*\
 * Use this data structure to store the value of any run-time parameters.
\*---------------------------------------------------------------------------*/

typedef struct USER_PARAMETERS{
   /* Name of file containingthe instance data */
   char             infile[MAX_FILE_NAME_LENGTH + 1];
  //sa:start
   char             infile2[MAX_FILE_NAME_LENGTH + 1];
  //sa:end
}user_parameters;

/*---------------------------------------------------------------------------*\
 * Use this data structure to store the instance data after it is read in.
\*---------------------------------------------------------------------------*/

typedef struct USER_PROBLEM{
   int              colnum; /* Number of rows in base matrix */
   int              rownum; /* Number of columns in base matrix */
   user_parameters  par;    /* Parameters */
//sa:start
    int              numvars;
    int              numcons;
    double           *TotalRhs;
    int              *ColInd;
    int              *RowInd;
    int              n2;
    int              m2;
    int              *lowerColInd;
    int              *lowerRowInd;
    double           *obj;
    int              lowerObjSense;
    int              nzt;
    int              nzlow;
    int              n,m,nz;
    int              *column_starts;
    int              *matrix_indices;
    double           *matrix_values;
    double           *lb;
    double           *ub;
    char             *sense;
    double           *rngval;
    char             *is_int;
    int              *upcolumn_starts;
    int              *matrix_upindices;
    double           *matrix_upvalues;
    int              len_ul;
    int              *IsInt;
    char             type;
    double           *intCosts;
    double           InterBudget;
   //sa:end

}user_problem;
//int user_read_data PROTO((user_problem *prob, char *infile));
int user_read_data PROTO((sym_environment *env,user_problem *prob, char *infile, char *infile2));

#endif
