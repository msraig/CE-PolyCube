#pragma warning( disable : 4477 4018 4267 4244 4838)
#include "ConvexQuadOptimization.h"
#include <stdio.h>
bool solveConvexQuadPorgramming_mosek(std::vector<MSKboundkeye>& bkc, std::vector<double>& blc, std::vector<double>& buc,/* Bounds on constraints. */
	std::vector<MSKboundkeye>& bkx, std::vector<double>& blx, std::vector<double>& bux,/* Bounds on variables. */
	std::vector<MSKlidxt>& aptrb, std::vector<MSKidxt>& asub, std::vector<double>& aval, 
	std::vector<MSKidxt>& qsubi, std::vector<MSKidxt>& qsubj, std::vector<double>& qval, std::vector<double>& c,
	std::vector<double>& XX)
{
	int NUMCON = bkc.size(); int NUMVAR, NUMANZ;
	if (NUMCON == 0)
	{
		NUMVAR = c.size();
		NUMANZ = 0;
	}
	else
	{
		NUMVAR = aptrb.size() - 1;
		NUMANZ = aptrb.back();
	}
	int NUMQNZ = qval.size();

	MSKidxt       i,j;
	XX.resize(NUMVAR);

	MSKenv_t      env;
	MSKtask_t     task;
	MSKrescodee   r;

	r = MSK_makeenv(&env,NULL);
	if ( r==MSK_RES_OK )
	{
		/* Directs the log stream to the 'printstr' function. */
		MSK_linkfunctoenvstream(env,
			MSK_STREAM_LOG,
			NULL,
			printstr);
	}

	bool s = false;
	/* Initialize the environment. */   
	//r = MSK_initenv(env); //diabled by yangliu, MOSKE 8 removes this function.
	if ( r == MSK_RES_OK )
	{ 
		/* Create the optimization task. */
		r = MSK_maketask(env,NUMCON,NUMVAR,&task);

		if ( r==MSK_RES_OK )
		{
			//r = MSK_linkfunctotaskstream(task,MSK_STREAM_LOG,NULL,printstr);
      
			/* Give MOSEK an estimate of the size of the input data. 
			 This is done to increase the speed of inputting data. 
			 However, it is optional. */
			if (r == MSK_RES_OK)
				r = MSK_putmaxnumvar(task,NUMVAR);
    
			if (r == MSK_RES_OK)
				r = MSK_putmaxnumcon(task,NUMCON);
      
			if (r == MSK_RES_OK)
				r = MSK_putmaxnumanz(task,NUMANZ);
  
			/* Append 'NUMCON' empty constraints.
			The constraints will initially have no bounds. */
			if ( r == MSK_RES_OK )
				r = MSK_appendcons(task, NUMCON);
  
			/* Append 'NUMVAR' variables.
			The variables will initially be fixed at zero (x=0). */
			if ( r == MSK_RES_OK )
				r = MSK_appendvars(task, NUMVAR);
  
			/* Optionally add a constant term to the objective. */
			if ( r ==MSK_RES_OK )
				r = MSK_putcfix(task,0.0);

			if (r == MSK_RES_OK)
				r = MSK_putintparam(task, MSK_IPAR_NUM_THREADS, 8);

			/*if (r == MSK_RES_OK)
			r = MSK_putdouparam(task,MSK_DPAR_INTPNT_CO_TOL_INFEAS ,1e-8);

			if (r == MSK_RES_OK)
			r = MSK_putdouparam(task,MSK_DPAR_INTPNT_CO_TOL_MU_RED ,1e-20);*/
			for(j=0; j<NUMVAR && r == MSK_RES_OK; ++j)
			{
				/* Set the linear term c_j in the objective.*/  
				if(r == MSK_RES_OK)
					r = MSK_putcj(task,j,c[j]);

				/* Set the bounds on variable j.
				blx[j] <= x_j <= bux[j] */
				if(r == MSK_RES_OK)
					r = MSK_putvarbound(task,
										j,           /* Index of variable.*/
										bkx[j],      /* Bound key.*/
										blx[j],      /* Numerical value of lower bound.*/
										bux[j]);     /* Numerical value of upper bound.*/
  
				if (NUMCON > 0)
				{
					/* Input column j of A */
					if (r == MSK_RES_OK)
						r = MSK_putacol(task,
						j,                 /* Variable (column) index.*/
						aptrb[j + 1] - aptrb[j], /* Number of non-zeros in column j.*/
						&asub[0] + aptrb[j],     /* Pointer to row indexes of column j.*/
						&aval[0] + aptrb[j]);    /* Pointer to Values of column j.*/
				}
			}
  
			/* Set the bounds on constraints.
			for i=1, ...,NUMCON : blc[i] <= constraint i <= buc[i] */
			for(i=0; i<NUMCON && r==MSK_RES_OK; ++i)
			{
				r = MSK_putconbound(task,
									i,           /* Index of constraint.*/
									bkc[i],      /* Bound key.*/
									blc[i],      /* Numerical value of lower bound.*/
									buc[i]);     /* Numerical value of upper bound.*/
			}

			if ( r==MSK_RES_OK )
			{
				/*
					* The lower triangular part of the Q
					* matrix in the objective is specified.
					*/

				/* Input the Q for the objective. */

				r = MSK_putqobj(task,NUMQNZ,&qsubi[0],&qsubj[0],&qval[0]);
			}

		
			//MSK_IPAR_INTPNT_NUM_THREADS

			if ( r==MSK_RES_OK )
			{
				MSKrescodee trmcode;

				/* Run optimizer */
				r = MSK_optimizetrm(task,&trmcode);

				/* Print a summary containing information
					about the solution for debugging purposes*/
				MSK_solutionsummary (task,MSK_STREAM_LOG);
        
				if ( r==MSK_RES_OK )
				{
					MSKsolstae solsta;
					int j;
          
					MSK_getsolsta(task,MSK_SOL_ITR,&solsta);

					switch (solsta)
					{
					case MSK_SOL_STA_OPTIMAL:
					case MSK_SOL_STA_NEAR_OPTIMAL:
						MSK_getxx(task,
							MSK_SOL_ITR,    /* Request the interior solution. */
							&XX[0]);

						printf("Optimal primal solution\n");
						s = true;
						break;
					case MSK_SOL_STA_DUAL_INFEAS_CER:
					case MSK_SOL_STA_PRIM_INFEAS_CER:
					case MSK_SOL_STA_NEAR_DUAL_INFEAS_CER:
					case MSK_SOL_STA_NEAR_PRIM_INFEAS_CER:
						printf("Primal or dual infeasibility certificate found.\n");
						break;

					case MSK_SOL_STA_UNKNOWN:
						printf("The status of the solution could not be determined.\n");
						break;
					default:
						printf("Other solution status.");
						break;
					}
				}
				else
				{
					printf("Error while optimizing.\n");
				}
			}

			if (r != MSK_RES_OK)
			{
				/* In case of an error print error code and description. */      
				char symname[MSK_MAX_STR_LEN];
				char desc[MSK_MAX_STR_LEN];

				printf("An error occurred while optimizing.\n");     
				MSK_getcodedesc (r,
									symname,
									desc);
				printf("Error %s - '%s'\n",symname,desc);
			}
		}
	}
	MSK_deletetask(&task);
	MSK_deleteenv(&env);
	return s;
}

bool solveConvexQuadPorgramming_mosek_integer(std::vector<MSKboundkeye>& bkc, std::vector<double>& blc, std::vector<double>& buc,/* Bounds on constraints. */
	std::vector<MSKboundkeye>& bkx, std::vector<double>& blx, std::vector<double>& bux,/* Bounds on variables. */
	std::vector<MSKlidxt>& aptrb, std::vector<MSKidxt>& asub, std::vector<double>& aval,
	std::vector<MSKidxt>& qsubi, std::vector<MSKidxt>& qsubj, std::vector<double>& qval, std::vector<double>& c,
	std::vector<double>& XX)
{
	int NUMCON = bkc.size(); int NUMVAR, NUMANZ;
	if (NUMCON == 0)
	{
		NUMVAR = c.size();
		NUMANZ = 0;
	}
	else
	{
		NUMVAR = aptrb.size() - 1;
		NUMANZ = aptrb.back();
	}
	int NUMQNZ = qval.size();

	MSKidxt       i, j;
	XX.resize(NUMVAR);

	MSKenv_t      env;
	MSKtask_t     task;
	MSKrescodee   r;

	r = MSK_makeenv(&env, NULL);
	if (r == MSK_RES_OK)
	{
		/* Directs the log stream to the 'printstr' function. */
		MSK_linkfunctoenvstream(env,
			MSK_STREAM_LOG,
			NULL,
			printstr);
	}

	bool s = false;
	/* Initialize the environment. */
	//r = MSK_initenv(env); //diabled by yangliu, MOSKE 8 removes this function.
	if (r == MSK_RES_OK)
	{
		/* Create the optimization task. */
		r = MSK_maketask(env, NUMCON, NUMVAR, &task);

		if (r == MSK_RES_OK)
		{
			r = MSK_linkfunctotaskstream(task,MSK_STREAM_LOG,NULL,printstr);

			/* Give MOSEK an estimate of the size of the input data.
			 This is done to increase the speed of inputting data.
			 However, it is optional. */
			if (r == MSK_RES_OK)
				r = MSK_putmaxnumvar(task, NUMVAR);

			if (r == MSK_RES_OK)
				r = MSK_putmaxnumcon(task, NUMCON);

			if (r == MSK_RES_OK)
				r = MSK_putmaxnumanz(task, NUMANZ);

			/* Append 'NUMCON' empty constraints.
			The constraints will initially have no bounds. */
			if (r == MSK_RES_OK)
				r = MSK_appendcons(task, NUMCON);

			/* Append 'NUMVAR' variables.
			The variables will initially be fixed at zero (x=0). */
			if (r == MSK_RES_OK)
				r = MSK_appendvars(task, NUMVAR);

			/* Optionally add a constant term to the objective. */
			if (r == MSK_RES_OK)
				r = MSK_putcfix(task, 0.0);

			if (r == MSK_RES_OK)
				r = MSK_putintparam(task, MSK_IPAR_NUM_THREADS, 8);

			/*if (r == MSK_RES_OK)
			r = MSK_putdouparam(task,MSK_DPAR_INTPNT_CO_TOL_INFEAS ,1e-8);

			if (r == MSK_RES_OK)
			r = MSK_putdouparam(task,MSK_DPAR_INTPNT_CO_TOL_MU_RED ,1e-20);*/
			for (j = 0; j < NUMVAR && r == MSK_RES_OK; ++j)
			{
				/* Set the linear term c_j in the objective.*/
				if (r == MSK_RES_OK)
					r = MSK_putcj(task, j, c[j]);

				/* Set the bounds on variable j.
				blx[j] <= x_j <= bux[j] */
				if (r == MSK_RES_OK)
					r = MSK_putvarbound(task,
						j,           /* Index of variable.*/
						bkx[j],      /* Bound key.*/
						blx[j],      /* Numerical value of lower bound.*/
						bux[j]);     /* Numerical value of upper bound.*/

				if (NUMCON > 0)
				{
					/* Input column j of A */
					if (r == MSK_RES_OK)
						r = MSK_putacol(task,
							j,                 /* Variable (column) index.*/
							aptrb[j + 1] - aptrb[j], /* Number of non-zeros in column j.*/
							&asub[0] + aptrb[j],     /* Pointer to row indexes of column j.*/
							&aval[0] + aptrb[j]);    /* Pointer to Values of column j.*/
				}
			}

			/* Set the bounds on constraints.
			for i=1, ...,NUMCON : blc[i] <= constraint i <= buc[i] */
			for (i = 0; i < NUMCON && r == MSK_RES_OK; ++i)
			{
				r = MSK_putconbound(task,
					i,           /* Index of constraint.*/
					bkc[i],      /* Bound key.*/
					blc[i],      /* Numerical value of lower bound.*/
					buc[i]);     /* Numerical value of upper bound.*/
			}


			

			

			


			if (r == MSK_RES_OK)
			{
				/*
					* The lower triangular part of the Q
					* matrix in the objective is specified.
					*/

					/* Input the Q for the objective. */

				r = MSK_putqobj(task, NUMQNZ, &qsubi[0], &qsubj[0], &qval[0]);
			}

			//specify integer variables
			for (i = 0; i < NUMVAR && r == MSK_RES_OK; i++)
			{
				r = MSK_putvartype(task, i, MSK_VAR_TYPE_INT);
			}

			/* Construct an initial feasible solution from the
			   values of the integer variables specified */
			if (r == MSK_RES_OK)
			r = MSK_putintparam(task, MSK_IPAR_MIO_CONSTRUCT_SOL, MSK_ON);

			if (r == MSK_RES_OK)
			{
			//double xx_initial[] = { 0.0, 2.0, 0.0 };
			double *xx_initial = c.data();

			/* Assign values 0,2,0 to integer variables */
			r = MSK_putxxslice(task, MSK_SOL_ITG, 0, c.size(), xx_initial);
			}

			if (r == MSK_RES_OK)
				r = MSK_putobjsense(task,
					MSK_OBJECTIVE_SENSE_MINIMIZE);

			if (r == MSK_RES_OK)
				/* Set max solution time */
				r = MSK_putdouparam(task,
					MSK_DPAR_MIO_MAX_TIME,
					60.0);


			//MSK_IPAR_INTPNT_NUM_THREADS


			if (r == MSK_RES_OK)
			{
				MSKrescodee trmcode;

				/* Run optimizer */
				r = MSK_optimizetrm(task, &trmcode);

				/* Print a summary containing information
				   about the solution for debugging purposes*/
				MSK_solutionsummary(task, MSK_STREAM_MSG);

				if (r == MSK_RES_OK)
				{
					MSKint32t  j;
					MSKsolstae solsta;
					//double     *xx = NULL;

					MSK_getsolsta(task, MSK_SOL_ITG, &solsta);
					
					//xx = (double *)calloc(NUMVAR, sizeof(double));
					if (r == MSK_RES_OK)
					{
						switch (solsta)
						{

						case MSK_SOL_STA_OPTIMAL:
						case MSK_SOL_STA_NEAR_OPTIMAL:
							MSK_getxx(task,
								MSK_SOL_ITR,    /* Request the interior solution. */
								&XX[0]);

							printf("Optimal primal solution\n");
							/*for (j = 0; j < NUMVAR; ++j)
								printf("x[%d]: %e\n", j, xx[j]);*/

							break;
						case MSK_SOL_STA_INTEGER_OPTIMAL:
						case MSK_SOL_STA_NEAR_INTEGER_OPTIMAL:
							MSK_getxx(task,
								MSK_SOL_ITG,    /* Request the integer solution. */
								&XX[0]);

							printf("Optimal integer solution.\n");
							s = true;
							/*for (j = 0; j < NUMVAR; ++j)
								printf("x[%d]: %e\n", j, xx[j]);*/
							break;
						case MSK_SOL_STA_PRIM_FEAS:
							/* A feasible but not necessarily optimal solution was located. */
							MSK_getxx(task, MSK_SOL_ITG, &XX[0]);

							printf("Feasible solution.\n");
							/*for (j = 0; j < NUMVAR; ++j)
								printf("x[%d]: %e\n", j, xx[j]);*/
							break;
						case MSK_SOL_STA_UNKNOWN:
						{
							MSKprostae prosta;
							MSK_getprosta(task, MSK_SOL_ITG, &prosta);
							switch (prosta)
							{
							case MSK_PRO_STA_PRIM_INFEAS_OR_UNBOUNDED:
								printf("Problem status Infeasible or unbounded\n");
								break;
							case MSK_PRO_STA_PRIM_INFEAS:
								printf("Problem status Infeasible.\n");
								break;
							case MSK_PRO_STA_UNKNOWN:
								printf("Problem status unknown.\n");
								break;
							default:
								printf("Other problem status.");
								break;
							}
						}
						break;
						default:
							printf("Other solution status.");
							break;
						}
					}
					else
					{
						r = MSK_RES_ERR_SPACE;
					}
					//free(xx);
				}
			}


			//if (r == MSK_RES_OK)
			//{
			//	MSKrescodee trmcode;

			//	/* Run optimizer */
			//	r = MSK_optimizetrm(task, &trmcode);

			//	/* Print a summary containing information
			//		about the solution for debugging purposes*/
			//	MSK_solutionsummary(task, MSK_STREAM_LOG);

			//	if (r == MSK_RES_OK)
			//	{
			//		MSKsolstae solsta;
			//		int j;

			//		MSK_getsolsta(task, MSK_SOL_ITR, &solsta);

			//		switch (solsta)
			//		{
			//		case MSK_SOL_STA_OPTIMAL:
			//		case MSK_SOL_STA_NEAR_OPTIMAL:
			//			MSK_getxx(task,
			//				MSK_SOL_ITR,    /* Request the interior solution. */
			//				&XX[0]);

			//			printf("Optimal primal solution\n");
			//			s = true;
			//			break;
			//		//case MSK_SOL_STA_DUAL_INFEAS_CER:
			//		//case MSK_SOL_STA_PRIM_INFEAS_CER:
			//		//case MSK_SOL_STA_NEAR_DUAL_INFEAS_CER:
			//		//case MSK_SOL_STA_NEAR_PRIM_INFEAS_CER:
			//		//	printf("Primal or dual infeasibility certificate found.\n");
			//		//	break;

			//		//case MSK_SOL_STA_UNKNOWN:
			//		//	printf("The status of the solution could not be determined.\n");
			//		//	break;
			//		//default:
			//		//	printf("Other solution status.\n");
			//		//	break;

			//		case MSK_SOL_STA_INTEGER_OPTIMAL:
			//		case MSK_SOL_STA_NEAR_INTEGER_OPTIMAL:
			//			MSK_getxx(task,
			//				MSK_SOL_ITR,    /* Request the interior solution. */
			//				&XX[0]);

			//			printf("Optimal primal integer solution\n");
			//			s = true;
			//			break;
			//		case MSK_SOL_STA_PRIM_FEAS:
			//			/* A feasible but not necessarily optimal solution was located. */
			//			MSK_getxx(task, MSK_SOL_ITG, &XX[0]);

			//			printf("Feasible solution.\n");
			//			/*for (j = 0; j < numvar; ++j)
			//				printf("x[%d]: %e\n", j, xx[j]);*/
			//			break;
			//		case MSK_SOL_STA_UNKNOWN:
			//		{
			//			printf("Status Unknown. \n");
			//			MSKprostae prosta;
			//			MSK_getprosta(task, MSK_SOL_ITG, &prosta);
			//			switch (prosta)
			//			{
			//			case MSK_PRO_STA_PRIM_INFEAS_OR_UNBOUNDED:
			//				printf("Problem status Infeasible or unbounded\n");
			//				break;
			//			case MSK_PRO_STA_PRIM_INFEAS:
			//				printf("Problem status Infeasible.\n");
			//				break;
			//			case MSK_PRO_STA_UNKNOWN:
			//				printf("Problem status unknown.\n");
			//				break;
			//			default:
			//				printf("Other problem status.");
			//				break;
			//			}
			//		}
			//		break;
			//		default:
			//			printf("Other solution status.");
			//			break;
			//		
			//		}
			//	}
			//	else
			//	{
			//		printf("Error while optimizing.\n");
			//	}
			//}

			if (r != MSK_RES_OK)
			{
				/* In case of an error print error code and description. */
				char symname[MSK_MAX_STR_LEN];
				char desc[MSK_MAX_STR_LEN];

				printf("An error occurred while optimizing.\n");
				MSK_getcodedesc(r,
					symname,
					desc);
				printf("Error %s - '%s'\n", symname, desc);
			}
		}
	}
	MSK_deletetask(&task);
	MSK_deleteenv(&env);
	return s;
}


//static void MSKAPI printstr(void *handle,
//	const char str[])
//{
//	printf("%s", str);
//} /* printstr */

bool mip_test()
{
	const MSKint32t numvar = 2,
		numcon = 2;

	double       c[] = { 1.0, 0.64 };
	MSKboundkeye bkc[] = { MSK_BK_UP,    MSK_BK_LO };
	double       blc[] = { -MSK_INFINITY, -4.0 };
	double       buc[] = { 250.0,        MSK_INFINITY };

	MSKboundkeye bkx[] = { MSK_BK_LO,    MSK_BK_LO };
	double       blx[] = { 0.0,          0.0 };
	double       bux[] = { MSK_INFINITY, MSK_INFINITY };


	MSKint32t    aptrb[] = { 0, 2 },
		aptre[] = { 2, 4 },
		asub[] = { 0,    1,   0,    1 };
	double       aval[] = { 50.0, 3.0, 31.0, -2.0 };
	MSKint32t    i, j;
	MSKint32t     qsubi[2];
	MSKint32t     qsubj[2];
	double        qval[2];


	MSKenv_t     env = NULL;
	MSKtask_t    task = NULL;
	MSKrescodee  r;

	/* Create the mosek environment. */
	r = MSK_makeenv(&env, NULL);

	/* Check if return code is ok. */
	if (r == MSK_RES_OK)
	{
		/* Create the optimization task. */
		r = MSK_maketask(env, 0, 0, &task);

		if (r == MSK_RES_OK)
			r = MSK_linkfunctotaskstream(task, MSK_STREAM_LOG, NULL, printstr);

		/* Append 'numcon' empty constraints.
		 The constraints will initially have no bounds. */
		if (r == MSK_RES_OK)
			r = MSK_appendcons(task, numcon);

		/* Append 'numvar' variables.
		 The variables will initially be fixed at zero (x=0). */
		if (r == MSK_RES_OK)
			r = MSK_appendvars(task, numvar);

		/* Optionally add a constant term to the objective. */
		if (r == MSK_RES_OK)
			r = MSK_putcfix(task, 0.0);
		for (j = 0; j < numvar && r == MSK_RES_OK; ++j)
		{
			/* Set the linear term c_j in the objective.*/
			/*if (r == MSK_RES_OK)
				r = MSK_putcj(task, j, c[j]);*/

				/* Set the bounds on variable j.
				 blx[j] <= x_j <= bux[j] */
			if (r == MSK_RES_OK)
				r = MSK_putvarbound(task,
					j,           /* Index of variable.*/
					bkx[j],      /* Bound key.*/
					blx[j],      /* Numerical value of lower bound.*/
					bux[j]);     /* Numerical value of upper bound.*/

/* Input column j of A */
			if (r == MSK_RES_OK)
				r = MSK_putacol(task,
					j,                 /* Variable (column) index.*/
					aptre[j] - aptrb[j], /* Number of non-zeros in column j.*/
					asub + aptrb[j],   /* Pointer to row indexes of column j.*/
					aval + aptrb[j]);  /* Pointer to Values of column j.*/

		}

		/* Set the bounds on constraints.
		   for i=1, ...,numcon : blc[i] <= constraint i <= buc[i] */
		for (i = 0; i < numcon && r == MSK_RES_OK; ++i)
			r = MSK_putconbound(task,
				i,           /* Index of constraint.*/
				bkc[i],      /* Bound key.*/
				blc[i],      /* Numerical value of lower bound.*/
				buc[i]);     /* Numerical value of upper bound.*/



		if (r == MSK_RES_OK)
		{
			/*
			 * The lower triangular part of the Q
			 * matrix in the objective is specified.
			 */

			qsubi[0] = 0;   qsubj[0] = 0;  qval[0] = 2.0;
			qsubi[1] = 1;   qsubj[1] = 1;  qval[1] = 2.0;
			/*qsubi[2] = 2;   qsubj[2] = 0;  qval[2] = -1.0;
			qsubi[3] = 2;   qsubj[3] = 2;  qval[3] = 2.0;*/

			/* Input the Q for the objective. */

			r = MSK_putqobj(task, 2, qsubi, qsubj, qval);
		}


		/* Specify integer variables. */
		for (j = 0; j < numvar && r == MSK_RES_OK; ++j)
			r = MSK_putvartype(task, j, MSK_VAR_TYPE_INT);

		if (r == MSK_RES_OK)
			r = MSK_putobjsense(task,
				MSK_OBJECTIVE_SENSE_MINIMIZE);

		if (r == MSK_RES_OK)
			/* Set max solution time */
			r = MSK_putdouparam(task,
				MSK_DPAR_MIO_MAX_TIME,
				60.0);

		if (r == MSK_RES_OK)
		{
			MSKrescodee trmcode;

			/* Run optimizer */
			r = MSK_optimizetrm(task, &trmcode);

			/* Print a summary containing information
			   about the solution for debugging purposes*/
			MSK_solutionsummary(task, MSK_STREAM_MSG);

			if (r == MSK_RES_OK)
			{
				MSKint32t  j;
				MSKsolstae solsta;
				double     *xx = NULL;

				MSK_getsolsta(task, MSK_SOL_ITG, &solsta);

				xx = (double *)calloc(numvar, sizeof(double));
				if (xx)
				{
					switch (solsta)
					{

					case MSK_SOL_STA_OPTIMAL:
					case MSK_SOL_STA_NEAR_OPTIMAL:
						MSK_getxx(task,
							MSK_SOL_ITR,    /* Request the interior solution. */
							xx);

						printf("Optimal primal solution\n");
						for (j = 0; j < numvar; ++j)
							printf("x[%d]: %e\n", j, xx[j]);

						break;
					case MSK_SOL_STA_INTEGER_OPTIMAL:
					case MSK_SOL_STA_NEAR_INTEGER_OPTIMAL:
						MSK_getxx(task,
							MSK_SOL_ITG,    /* Request the integer solution. */
							xx);

						printf("Optimal integer solution.\n");
						for (j = 0; j < numvar; ++j)
							printf("x[%d]: %e\n", j, xx[j]);
						break;
					case MSK_SOL_STA_PRIM_FEAS:
						/* A feasible but not necessarily optimal solution was located. */
						MSK_getxx(task, MSK_SOL_ITG, xx);

						printf("Feasible solution.\n");
						for (j = 0; j < numvar; ++j)
							printf("x[%d]: %e\n", j, xx[j]);
						break;
					case MSK_SOL_STA_UNKNOWN:
					{
						MSKprostae prosta;
						MSK_getprosta(task, MSK_SOL_ITG, &prosta);
						switch (prosta)
						{
						case MSK_PRO_STA_PRIM_INFEAS_OR_UNBOUNDED:
							printf("Problem status Infeasible or unbounded\n");
							break;
						case MSK_PRO_STA_PRIM_INFEAS:
							printf("Problem status Infeasible.\n");
							break;
						case MSK_PRO_STA_UNKNOWN:
							printf("Problem status unknown.\n");
							break;
						default:
							printf("Other problem status.");
							break;
						}
					}
					break;
					default:
						printf("Other solution status.");
						break;
					}
				}
				else
				{
					r = MSK_RES_ERR_SPACE;
				}
				free(xx);
			}
		}

		if (r != MSK_RES_OK)
		{
			/* In case of an error print error code and description. */
			char symname[MSK_MAX_STR_LEN];
			char desc[MSK_MAX_STR_LEN];

			printf("An error occurred while optimizing.\n");
			MSK_getcodedesc(r,
				symname,
				desc);
			printf("Error %s - '%s'\n", symname, desc);
		}

		MSK_deletetask(&task);
	}
	MSK_deleteenv(&env);

	printf("Return code: %d.\n", r);
	return (r);
} /* main */
