/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2014 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   probdata_mochila.c
 * @brief  Problem data for mochila problem
 * @author Edna Hoshino (based on template codified by Timo Berthold and Stefan Heinz)
 *
 * This file handles the main problem data used in that project. For more details see \ref PROBLEMDATA page.
 *
 * @page PROBLEMDATA Main problem data
 *
 * The problem data is accessible in all plugins. The function SCIPgetProbData() returns the pointer to that
 * structure. We use this data structure to store all the information of the mochila problem. Since this structure is
 * not visible in the other plugins, we implemented setter and getter functions to access this data. The problem data
 * structure SCIP_ProbData is shown below.
 *
 * \code
 *  ** @brief Problem data which is accessible in all places
 *  *
 *  *   This problem data is used to store the input of the mochila instance, all variables which are created, and all
 *  *   constraints.
 *  *
 * struct SCIP_ProbData
 * {
 *    const_char*           probname;     **< problem name *
 *    SCIP_VAR**            vars;         **< all variables of the problem *
 *    SCIP_CONS**           conss;        **< all constraints  *
 *    int                   nvars;        **< size of vars *
 *    int                   ncons;        **< number of constraints *
 *    instanceT*            I;            **< instance of knapsack *
 * };
 * \endcode
 *
 * The function SCIPprobdataCreate(), which is called in the \ref loadProblem.c after the input file was
 * parsed, initializes the problem data structure and creates the problem in the SCIP environment. 
 * This code can be easily modified to dealt with other problem. Just take a look in the code marked with "TODO".
 * See the body of the function SCIPprobdataCreate() for more details.
 *
 * The following problem refers to the Knapsack problem
 *
 *  \f[
 *  \begin{array}[t]{rll}
 *       \max & \displaystyle \sum_{i \in I} v_i x_i \\
 *        & \\
 *        subject \ to & \displaystyle \sum_{i \in I} w_i x_i <= C \\
 *        & \\
 *        & x_i \in \{0,1\} & \quad \forall i \in I \\
 *  \end{array}
 * \f]
 *
 * where \f$ x_i \f$ for \f$i\in I\f$ are binary variables and \f$w_i\f$ and \f$v_i\f$ are the weight and value of item \f$ i\f$, respectively.
 *
 * A list of all interface methods can be found in probdata_mochila.h.
 **/
/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
/*
#define SCIP_DEBUG 
*/

#include <assert.h>
#include <string.h>

#include "scip/scipdefplugins.h"

#include "probdata_mochila.h"

/**@name Local methods
 *
 * @{
 */

/** creates problem data */
static
SCIP_RETCODE probdataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA**       probdata,           /**< pointer to problem data */
   const char*           probname,           /**< problem name */
   SCIP_VAR**            vars,               /**< array of variables */
   SCIP_CONS**           conss,              /**< all constraints */
   int                   nvars,              /**< size of vars */
   int                   ncons,              /**< number of constraints */
   instanceT*            I                   /**< pointer to the instance data */
   )
{
   assert(scip != NULL);
   assert(probdata != NULL);

   /* allocate memory */
   SCIP_CALL( SCIPallocMemory(scip, probdata) );

   if( nvars > 0 )
   {
      /* copy variable array */
      (*probdata)->vars=vars;
      SCIP_CALL( SCIPduplicateMemoryArray(scip, &(*probdata)->vars, vars, nvars) ); /* NEEDED for transformed problem*/
   }
   else
      (*probdata)->vars = NULL;
   /* duplicate arrays */
   SCIP_CALL( SCIPduplicateMemoryArray(scip, &(*probdata)->conss, conss, ncons) ); /* NEEDED for transformed problem */

   (*probdata)->I=I;
   (*probdata)->nvars = nvars;
   (*probdata)->ncons = ncons;
   (*probdata)->probname = probname;

   return SCIP_OKAY;
}

/** frees the memory of the given problem data */
static
SCIP_RETCODE probdataFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA**       probdata,            /**< pointer to problem data */
   int transformed
   )
{
   int i;

   assert(scip != NULL);
   assert(probdata != NULL);

   /* release all variables */
   for( i = 0; i < (*probdata)->nvars; ++i )
   {
      SCIP_CALL( SCIPreleaseVar(scip, &(*probdata)->vars[i]) );
   }
  
   /* release all constraints */
   for( i = 0; i < (*probdata)->ncons; ++i )
   {
      SCIP_CALL( SCIPreleaseCons(scip, &(*probdata)->conss[i]) );
   }

   /* free memory of arrays */
   SCIPfreeMemoryArray(scip, &(*probdata)->conss);
   if(!transformed){
     //SCIPfreeMemoryArray(scip, &(*probdata)->vars);
     freeInstance((*probdata)->I);
   }
   /* free probdata */
   SCIPfreeMemory(scip, probdata);

   return SCIP_OKAY;
}


/**@} */

/**@name Callback methods of problem data
 *
 * @{
 */

/** frees user data of original problem (called when the original problem is freed) */
static
SCIP_DECL_PROBDELORIG(probdelorigMochila)
{
   SCIPdebugMessage("free original problem data\n");

   SCIP_CALL( probdataFree(scip, probdata, 0) );

   return SCIP_OKAY;
}

/** creates user data of transformed problem by transforming the original user problem data
 *  (called after problem was transformed) @@@ ?*/
static
SCIP_DECL_PROBTRANS(probtransMochila)
{
   /* create transform probdata */
   SCIP_CALL( probdataCreate(scip, targetdata, sourcedata->probname, sourcedata->vars, sourcedata->conss, sourcedata->nvars, sourcedata->ncons, sourcedata->I) );

   /* transform all constraints */
   SCIP_CALL( SCIPtransformConss(scip, (*targetdata)->ncons, (sourcedata)->conss, (*targetdata)->conss) );
   /* transform all variables */
   SCIP_CALL( SCIPtransformVars(scip, (*targetdata)->nvars, (sourcedata)->vars, (*targetdata)->vars) );

   return SCIP_OKAY;
}

/** frees user data of transformed problem (called when the transformed problem is freed) */
static
SCIP_DECL_PROBDELTRANS(probdeltransMochila)
{
   SCIPdebugMessage("free transformed problem data\n");

   SCIP_CALL( probdataFree(scip, probdata,1) );

   return SCIP_OKAY;
}

/** solving process initialization method of transformed data (called before the branch and bound process begins) */
static
SCIP_DECL_PROBINITSOL(probinitsolMochila)
{
   assert(probdata != NULL);

   return SCIP_OKAY;
}

/** solving process deinitialization method of transformed data (called before the branch and bound data is freed) */
static
SCIP_DECL_PROBEXITSOL(probexitsolMochila)
{
   assert(probdata != NULL);

   return SCIP_OKAY;
}

/**@} */


/**@name Interface methods
 *
 * @{
 */

/** sets up the problem data 
 * TODO: specific for the problem
*/
SCIP_RETCODE SCIPprobdataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           probname,           /**< problem name */
   instanceT*            I                   /**< instance of knapsack */
   )
{
   SCIP_PROBDATA* probdata;
   SCIP_CONS** conss;   // esse cara aqui é um vetor de restricoes
   SCIP_VAR** vars, * var;

   char name[SCIP_MAXSTRLEN];
   int i, j;   // i itens e j mochilas
   int ncons;  // n restricoes
   int nvars;  // n variaveis
   
   assert(scip != NULL);

   /* create problem in SCIP and add non-NULL callbacks via setter functions */
   SCIP_CALL( SCIPcreateProbBasic(scip, probname) );

   SCIP_CALL( SCIPsetProbDelorig(scip, probdelorigMochila) );
   SCIP_CALL( SCIPsetProbTrans(scip, probtransMochila) );
   SCIP_CALL( SCIPsetProbDeltrans(scip, probdeltransMochila) );
   SCIP_CALL( SCIPsetProbInitsol(scip, probinitsolMochila) );
   SCIP_CALL( SCIPsetProbExitsol(scip, probexitsolMochila) );

   /* set objective sense */
   SCIP_CALL( SCIPsetObjsense(scip, SCIP_OBJSENSE_MAXIMIZE) );

   /* TODO: tell SCIP that the objective will be always integral (it depends on the problem) */
   SCIP_CALL( SCIPsetObjIntegral(scip) );
  
   // alloc memory to create vars and cons - it is necessary to probdatacreate(), that will make a copy of them.
   SCIP_CALL( SCIPallocBufferArray(scip, &conss, I->m + I->n) );   // m restricoes. uma para cada mochila + n restricoes. uma para cada item (para garantir que cada item so seja colocado em uma mochila)
   SCIP_CALL( SCIPallocBufferArray(scip, &vars, (I->n) * (I->m)) );   // alocando agr n*m variaveis ao inves de n(pq agr a variavel eh x_i_j)
   
   ncons=0;  // incializando o numero de restricoes do modelo
   nvars=0;  // inicializando o numero de variaveis do modelo
   // TODO: configure vars and constraints ....   CONFIGURANDO AS VARIAVEIS E AS RESTRICOES

   // criando as restricoes de capacidade para cada uma das j mochilas (mochila multipla)
   for( j = 0; j < I->m; j++){
      // name = nome da restricao
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "capacity_%d", j);
      
      // essa funcao que cria a restricao                                             // esse trecho aqui eh equivalente a <=
      SCIP_CALL( SCIPcreateConsBasicLinear (scip, &conss[ncons], name, 0, NULL, NULL, -SCIPinfinity(scip), (double) I->C[j]) );
      SCIP_CALL( SCIPaddCons(scip, conss[ncons]) );
      ncons++; // eai o numero de restricoes aumenta para cada iteracao desse for
   }

     // criando a restrição de que cada item so pode ser colocado em uma unica mochila
   for( i = 0; i < I->n; i++){
      // name = nome da restricao
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "item_%d", i);
      //                                                                              // x_i_j <= 1
      SCIP_CALL( SCIPcreateConsBasicLinear (scip, &conss[ncons], name, 0, NULL, NULL, -SCIPinfinity(scip), 1.0) );
      SCIP_CALL( SCIPaddCons(scip, conss[ncons]) );
      ncons++;  
   }
   
   /* create one variable for each item i */
   // criando n*m variaveis, pq agr eu tenho j variaveis para cada item i
   for( i = 0; i < I->n; nvars++, ++i )
   {
      // com esse segundo for aqui eu consigo tanto fazer a variavel ter dois indices
      // quanto colocar cada variavel i na restricao de capacidade da mochila j
      for(j = 0; j < I->m; j++){

         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "x_%d_%d", i, j);

         /* create a basic variable object */
         // criando uma variavel do tipo basica
         SCIP_CALL( SCIPcreateVarBasic(scip, &var, name, 0.0, 1.0, (double) I->item[i].value, SCIP_VARTYPE_BINARY));

         assert(var != NULL);
         /* save the pointer to the created var */
         vars[nvars]=var;

         /* add variable to the problem */
         SCIP_CALL( SCIPaddVar(scip, var) );
      
         // colocando a variavel com aquele coeficiente (peso do item i) nas restriçoes (mochila multipla)
         // adicionando a variaval a restricao de capacidade j (mochila multipla)
         SCIP_CALL( SCIPaddCoefLinear(scip, conss[j], var, (double) I->item[i].weight) );

         // adicionando a variaval a restricao do item i (mochila multipla).
         SCIP_CALL( SCIPaddCoefLinear(scip, conss[I->m+i], var, 1.0 ));  // item i so pode ser colocado em uma mochila
         // pulando m restricoes para assim chegar na restricao correta
      }
      
   }   
      
   // TODO: ... after vars and constraints have been created, nothing more is necessary. Just do exactly as follows:
   /* create problem data */
   SCIP_CALL( probdataCreate(scip, &probdata, probname, vars, conss, nvars, ncons, I) );

#ifdef DEBUG_PROBDATA
   SCIP_CALL( SCIPwriteOrigProblem(scip, "mochila.lp", "lp", FALSE) ); /* save the problem in a file */
#endif
   /* set user problem data */
   SCIP_CALL( SCIPsetProbData(scip, probdata) );


   /* free local buffer arrays */
   for( i = 0; i < nvars; ++i )
     SCIP_CALL( SCIPreleaseVar(scip, &(vars[i])) );
   for( i = 0; i < ncons; ++i )
     SCIP_CALL( SCIPreleaseCons(scip, &(conss[i])) );
   SCIPfreeBufferArray(scip, &conss);
   SCIPfreeBufferArray(scip, &vars);

   return SCIP_OKAY;
}

instanceT* SCIPprobdataGetInstance(
   SCIP_PROBDATA*        probdata
   )
{
   return probdata->I;
}

/** returns array of all variables itemed in the way they got generated */
SCIP_VAR** SCIPprobdataGetVars(
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   return probdata->vars;
}

/** returns number of variables */
int SCIPprobdataGetNVars(
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   return probdata->nvars;
}
/** returns array of set partitioning constrains */
SCIP_CONS** SCIPprobdataGetConss(
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   return probdata->conss;
}

/** returns array of set partitioning constrains */
int SCIPprobdataGetNcons(
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   return probdata->ncons;
}
/** returns Probname of the instance */
const char* SCIPprobdataGetProbname(
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   return probdata->probname;
}

/**@} */
