/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2016 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   heur_aleatoria.c
 * @brief  aleatoria primal heuristic
 * @author Edna Hoshino (based on template provided by Tobias Achterberg)
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <time.h>

#include "probdata_mochila.h"
#include "parameters_mochila.h"
#include "heur_aleatoria.h"
#include "heur_problem.h"

//#define DEBUG_ALEATORIA 1
/* configuracao da heuristica */
#define HEUR_NAME             "aleatoria"
#define HEUR_DESC             "primal heuristic template"
#define HEUR_DISPCHAR         'a'
#define HEUR_PRIORITY         3 /* comeca pelas heuristicas de maior prioridade */
#define HEUR_FREQ             1 /* a cada 1 nivel da arvore de B&B */
#define HEUR_FREQOFS          0 /* comecando do nivel 0 */
#define HEUR_MAXDEPTH         -1 /* nivel max para chamar a heuristica. -1 = sem limites */
#define HEUR_TIMING           SCIP_HEURTIMING_AFTERNODE //SCIP_HEURTIMING_DURINGLPLOOP // SCIP_HEURTIMING_AFTERNODE /* chamado depois que o LP resolvido */
#define HEUR_USESSUBSCIP      FALSE  /**< does the heuristic use a secondary SCIP instance? */

#ifdef DEBUG
   #define PRINTF(...) printf(__VA_ARGS__)
#else
   #define PRINTF(...) 
#endif

/*
 * Data structures
 */

/* TODO: fill in the necessary primal heuristic data */

/** primal heuristic data */
/*struct SCIP_HeurData
{
};
*/

/*
 * Local methods
 */

/* put your local methods here, and declare them static */

/*
 * Callback methods of primal heuristic
 */

/* TODO: Implement all necessary primal heuristic methods. The methods with an #if 0 ... #else #define ... are optional */

/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
static
SCIP_DECL_HEURCOPY(heurCopyAleatoria)
{  /*lint --e{715}*/

   return SCIP_OKAY;
}

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeAleatoria)
{  /*lint --e{715}*/

   return SCIP_OKAY;
}


/** initialization method of primal heuristic (called after problem was transformed) */
static
SCIP_DECL_HEURINIT(heurInitAleatoria)
{  /*lint --e{715}*/


   return SCIP_OKAY;
}


/** deinitialization method of primal heuristic (called before transformed problem is freed) */
static
SCIP_DECL_HEUREXIT(heurExitAleatoria)
{  /*lint --e{715}*/

   return SCIP_OKAY;
}


/** solving process initialization method of primal heuristic (called when branch and bound process is about to begin) */
static
SCIP_DECL_HEURINITSOL(heurInitsolAleatoria)
{  /*lint --e{715}*/

   return SCIP_OKAY;
}


/** solving process deinitialization method of primal heuristic (called before branch and bound process data is freed) */
static
SCIP_DECL_HEUREXITSOL(heurExitsolAleatoria)
{  /*lint --e{715}*/

   return SCIP_OKAY;
}

static int numero_aleatorio(int n_cand){
   return rand() % n_cand;
}



/**
 * @brief Core of the aleatoria heuristic: it builds one solution for the problem by aleatoria procedure.
 *
 * @param scip problem
 * @param sol pointer to the solution structure where the solution wil be saved
 * @param heur pointer to the aleatoria heuristic handle (to contabilize statistics)
 * @return int 1 if solutions is found, 0 otherwise.
 */
int aleatoria(SCIP* scip, SCIP_SOL** sol, SCIP_HEUR* heur)
{
   int found, infeasible, nInSolution;
   unsigned int stored;
   int nvars;
   int *covered, n, m, nCovered = 0, custo, **cand, *nCands, s;
   SCIP_VAR *var, **solution, **varlist;
   //  SCIP* scip_cp;
   SCIP_Real valor, bestUb;
   SCIP_PROBDATA* probdata;
   int i, *residual;
   instanceT* I;
   
   found = 0;
   infeasible = 0;
   
#ifdef DEBUG_ALEATORIA
   printf("\n============== New aleatoria heur at node: %lld\n", SCIPnodeGetNumber(SCIPgetCurrentNode(scip)));
#endif

   /* recupera os dados do problema original*/
   probdata=SCIPgetProbData(scip);
   assert(probdata != NULL);

   nvars = SCIPprobdataGetNVars(probdata);
   varlist = SCIPprobdataGetVars(probdata);
   I = SCIPprobdataGetInstance(probdata);
   n = I->n;
   m = I->m;

   // // print ordered items
   // printf("\nGULOSA: Ordered items by value/weight\n");
   // for(i=0;i<n;i++){
   //   printf("\nGULOSA: Item %d: value=%d weight=%d value/weight=%lf", I->item[i].label, I->item[i].value, I->item[i].weight, (double) I->item[i].value/I->item[i].weight);
   //   SCIPprintVar(scip, varlist[i], NULL);
   //   printf("\n%f\n", SCIPvarGetLbLocal(varlist[i]));
   // }
   
    
   solution = (SCIP_VAR**) malloc(sizeof(SCIP_VAR*)*n);
   covered = (int*) calloc(n,sizeof(int));
   cand = (int**) malloc(m*sizeof(int*));

   for(int k = 0; k < m; k++){
      cand[k] = (int *) malloc(n*sizeof(int));
   }

   residual = (int *) malloc(m*sizeof(int));
   nCands = (int *) calloc(m,sizeof(int));  // inicializa todas as posicoes com 0


   for(int k = 0; k < m; k++){
      residual[k] = I->C[k];
   }

   nInSolution = 0;
   custo = 0;
   
   // first, select all variables already fixed in 1.0
  for(i=0;i<nvars;i++){
      var = varlist[i];
      
      if(SCIPvarGetLbLocal(var) > 1.0 - EPSILON){ // var >= 1.0
         printf("\n\nTESTEEEEEEEEEEEEEEEEEEEEEEEE\n");
        solution[nInSolution++]=var;
        //exit(0);
        residual[i%m] -= I->item[i].weight;
        covered[i%n]=1;
        nCovered++;
        custo += I->item[i].value;
        infeasible = residual < 0?1:0;
#ifdef DEBUG_ALEATORIA
        printf("\nSelected fixed var= %s. TotalItems=%d value=%d residual=%d infeasible=%d", SCIPvarGetName(var), nInSolution, custo, residual, infeasible);
#endif
      }
      else{ // discard items fixed in 0.0
        if(SCIPvarGetUbLocal(var) < EPSILON){ // var fixed in 0.0
        }

        else{
          if (i < n){ // include item i in cand list
            for(int k = 0; k < m; k++){

               if(I->item[i].weight <= residual[k] && !covered[i%n]){
                  cand[k][nCands[k]] = i;
                  nCands[k] += 1;
               }
            }
          }
        }
      }
   }


   // complete solution using items not fixed (not covered)
   for(i=0; i < m && nCovered < n; i++){
      // enquanto a capacidade atual da mochila i for > 0 E a mochila i ainda estiver itens candidatos
      while(residual[i] > 0 && nCands[i] > 0){
         printf("\nLISTA DOS ITENS COBERTOS:\n");
         for(int k = 0; k < n; k++){
            printf("%d ", covered[k]);
         }
         s = numero_aleatorio(nCands[i]);
         int aux = cand[i][s];
         printf("\n\nNUMERO ESCOLHIDO DA MOCHILA %d: %d\n", i, s);
         printf("AUXILIAR: %d\n", aux);
         printf("\nLISTA DE CANDIDATOS DA MOCHILA %d:\n", i);
         for(int j = 0; j < nCands[i]; j++){
            printf("%d ", cand[i][j]);
         }
        
         cand[i][s] = cand[i][nCands[i]-1];
         nCands[i] -= 1;

         printf("\nLISTA DE CANDIDATOS DA MOCHILA %d DEPOIS DO ITEM ESCOLHIDO:\n", i);
         for(int j = 0; j < nCands[i]; j++){
            printf("%d ", cand[i][j]);
         }

         
         if(!covered[aux] && I->item[aux].weight <= residual[i]){
            printf("\nITEM %d COLODADO NA MOCHILA: %d\n", aux, i);
            valor = I->item[aux].value;
            var = varlist[aux];

            solution[nInSolution++] = var;

            residual[i] -= I->item[aux].weight;
            covered[aux] = 1;
            nCovered++;
            custo += valor;
            infeasible = residual<0?1:0;
            #ifdef DEBUG_ALEATORIA
               printf("\n\nSelected var= %s. TotalItems=%d value item=%d value = %d residual=%d infeasible=%d\n", SCIPvarGetName(var), nInSolution, valor, custo, residual, infeasible);
            #endif
         }
         else{// desconsidere o item
            #ifdef DEBUG_ALEATORIA
               printf("\n\nNOT selected var= %s. TotalItems=%d value item=%d value=%d residual=%d infeasible=%d\n", SCIPvarGetName(varlist[selected]), nInSolution, valor, custo, residual, infeasible);
            #endif
         }
      }
   }
   if(!infeasible){
      /* create SCIP solution structure sol */
      SCIP_CALL( SCIPcreateSol(scip, sol, heur) );
      // save found solution in sol
      for(i=0;i<nInSolution;i++){
         var = solution[i];
         SCIP_CALL( SCIPsetSolVal(scip, *sol, var, 1.0) );
      }
      bestUb = SCIPgetPrimalbound(scip);
#ifdef DEBUG_ALEATORIA
      printf("\nFound solution...\n");
      //      SCIP_CALL( SCIPprintSol(scip, *sol, NULL, FALSE) );
      printf("\ninfeasible=%d value = %lf > bestUb = %lf? %d\n\n", infeasible, valor, bestUb, valor > bestUb + EPSILON);
#endif
      if(!infeasible && custo > bestUb + EPSILON){
#ifdef DEBUG_ALEATORIA
         printf("\nBest solution found...\n");
         SCIP_CALL( SCIPprintSol(scip, *sol, NULL, FALSE) );
#endif
         
         /* verificar se a solucao eh viavel e armazena */
         SCIP_CALL( SCIPtrySolMine(scip, *sol, TRUE, TRUE, FALSE, TRUE, &stored) );
         if( stored )
         {
#ifdef DEBUG_PRIMAL
            printf("\nSolution is feasible and was saved! Total of items = %d", nInSolution);
            SCIPdebugMessage("found feasible aleatoria solution:\n");
            SCIP_CALL( SCIPprintSol(scip, *sol, NULL, FALSE) );
#endif
            found = 1;
         }
         else{
           found = 0;
#ifdef DEBUG_ALEATORIA
           printf("\nCould not found\n. BestUb=%lf", bestUb);
#endif
         }
      }
   }
#ifdef DEBUG_ALEATORIA
   getchar();
#endif
   free(varlist);
   free(solution);
   free(covered);
   free(cand);
   return found;
}

/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecAleatoria)
{  /*lint --e{715}*/
   SCIP_SOL*             sol;                /**< solution to round */
   int nlpcands;

   assert(result != NULL);
   //   assert(SCIPhasCurrentNodeLP(scip));

   *result = SCIP_DIDNOTRUN;

   /* continue only if the LP is finished */
   if ( SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL )
      return SCIP_OKAY;

   /* continue only of the LP value is less than the cutoff bound */
   if( SCIPisGE(scip, SCIPgetLPObjval(scip), SCIPgetCutoffbound(scip)) )
      return SCIP_OKAY;


   /* check if there exists integer variables with fractionary values in the LP */
   SCIP_CALL( SCIPgetLPBranchCands(scip, NULL, NULL, NULL, &nlpcands, NULL, NULL) );
   //Fractional implicit integer variables are stored at the positions *nlpcands to *nlpcands + *nfrac - 1
  
   /* stop if the LP solution is already integer   */
   if ( nlpcands == 0 )
     return SCIP_OKAY;

   /* solve aleatoria */
   if(aleatoria(scip, &sol, heur)){
     *result = SCIP_FOUNDSOL;
   }
   else{
#ifdef DEBUG_PRIMAL
     printf("\nAleatoria could not find feasible solution!");      
#endif
   }
   return SCIP_OKAY;
}


/*
 * primal heuristic specific interface methods
 */

/** creates the aleatoria_crtp primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurAleatoria(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_HEUR* heur;

   /* create aleatoria primal heuristic data */
   heurdata = NULL;

   heur = NULL;

   /* include primal heuristic */
#if 0
   /* use SCIPincludeHeur() if you want to set all callbacks explicitly and realize (by getting compiler errors) when
    * new callbacks are added in future SCIP versions
    */
   SCIP_CALL( SCIPincludeHeur(scip, HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, param.heur_freq, param.heur_freqofs,
         param.heur_maxdepth, HEUR_TIMING, HEUR_USESSUBSCIP,
         heurCopyAleatoria, heurFreeAleatoria, heurInitAleatoria, heurExitAleatoria, heurInitsolAleatoria, heurExitsolAleatoria, heurExecAleatoria,
         heurdata) );
#else
   /* use SCIPincludeHeurBasic() plus setter functions if you want to set callbacks one-by-one and your code should
    * compile independent of new callbacks being added in future SCIP versions
    */
   SCIP_CALL( SCIPincludeHeurBasic(scip, &heur,
         HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, param.heur_round_freq, param.heur_round_freqofs,
         param.heur_round_maxdepth, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecAleatoria, heurdata) );

   assert(heur != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetHeurCopy(scip, heur, heurCopyAleatoria) );
   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreeAleatoria) );
   SCIP_CALL( SCIPsetHeurInit(scip, heur, heurInitAleatoria) );
   SCIP_CALL( SCIPsetHeurExit(scip, heur, heurExitAleatoria) );
   SCIP_CALL( SCIPsetHeurInitsol(scip, heur, heurInitsolAleatoria) );
   SCIP_CALL( SCIPsetHeurExitsol(scip, heur, heurExitsolAleatoria) );
#endif

   /* add aleatoria primal heuristic parameters */
   /* TODO: (optional) add primal heuristic specific parameters with SCIPaddTypeParam() here */

   return SCIP_OKAY;
}
