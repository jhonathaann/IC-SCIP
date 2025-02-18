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

/**@file   heur_grasp.c
 * @brief  grasp primal heuristic
 * @author Edna Hoshino (based on template provided by Tobias Achterberg)
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <time.h> 

#include "probdata_mochila.h"
#include "parameters_mochila.h"
#include "heur_grasp.h"
#include "heur_problem.h"
#include "problem.h"   // eu que inclui isso.

//#define DEBUG_GRASP 1
/* configuracao da heuristica */
#define HEUR_NAME             "grasp"
#define HEUR_DESC             "primal heuristic template"
#define HEUR_DISPCHAR         'g'
#define HEUR_PRIORITY         3 /* comeca pelas heuristicas de maior prioridade */
#define HEUR_FREQ             1 /* a cada 1 nivel da arvore de B&B */
#define HEUR_FREQOFS          0 /* comecando do nivel 0 */
#define HEUR_MAXDEPTH         10 /* nivel max para chamar a heuristica. -1 = sem limites */
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
SCIP_DECL_HEURCOPY(heurCopyGrasp)
{  /*lint --e{715}*/

   return SCIP_OKAY;
}

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeGrasp)
{  /*lint --e{715}*/

   return SCIP_OKAY;
}


/** initialization method of primal heuristic (called after problem was transformed) */
static
SCIP_DECL_HEURINIT(heurInitGrasp)
{  /*lint --e{715}*/


   return SCIP_OKAY;
}


/** deinitialization method of primal heuristic (called before transformed problem is freed) */
static
SCIP_DECL_HEUREXIT(heurExitGrasp)
{  /*lint --e{715}*/

   return SCIP_OKAY;
}


/** solving process initialization method of primal heuristic (called when branch and bound process is about to begin) */
static
SCIP_DECL_HEURINITSOL(heurInitsolGrasp)
{  /*lint --e{715}*/

   return SCIP_OKAY;
}


/** solving process deinitialization method of primal heuristic (called before branch and bound process data is freed) */
static
SCIP_DECL_HEUREXITSOL(heurExitsolGrasp)
{  /*lint --e{715}*/

   return SCIP_OKAY;
}

void min_max(itemType *candidatos, int n, int *max, int *min);

void cria_RCL(itemType *candidatos, itemType *RCL, int minimo, int maximo, int alpha, int n, int *n_RCL);

int numero_aleatorio(int n_cand);

void atualiza_candidatos(itemType *candidatos, int rotulo, int *n_cand, int capacidade_atual);

int verifica_solucao(solucao *itens_solucao, int label, int n);

/**
 * @brief Core of the grasp heuristic: it builds one solution for the problem by grasp procedure.
 *
 * @param scip problem
 * @param sol pointer to the solution structure where the solution wil be saved
 * @param heur pointer to the grasp heuristic handle (to contabilize statistics)
 * @return int 1 if solutions is found, 0 otherwise.
 */
int grasp(SCIP* scip, SCIP_SOL** sol, SCIP_HEUR* heur)
{
   int found, infeasible, nInSolution;
   unsigned int stored;
  int nvars;
   int *covered, n, m, custo, *cand, nCands, selected, s;
   SCIP_VAR *var, **solution, **varlist;
   //  SCIP* scip_cp;
   SCIP_Real valor, bestUb;
   SCIP_PROBDATA* probdata;
   int i, residual;
   instanceT* I;

   itemType *candidatos, *RCL;
   solucao *itens_solucao;

   int max_iteracoes, maximo, minimo, alpha = 0.7, n_cand = 0;
   int posicao_item_escolhido, n_RCL = 0, rotulo, n_solucao;
   
   found = 0;
   infeasible = 0;
   
#ifdef DEBUG_GRASP
   printf("\n============== New grasp heur at node: %lld\n", SCIPnodeGetNumber(SCIPgetCurrentNode(scip)));
#endif

   /* recupera os dados do problema original*/
   probdata=SCIPgetProbData(scip);
   assert(probdata != NULL);

   //nvars = SCIPprobdataGetNVars(probdata);
   varlist = SCIPprobdataGetVars(probdata);
   I = SCIPprobdataGetInstance(probdata);
   n = I->n;
   m = I->m;
    
   solution = (SCIP_VAR**) malloc(sizeof(SCIP_VAR*)*n);
   covered = (int*) calloc(n,sizeof(int));
   cand = (int*) malloc(n*sizeof(int));
   nCands = 0;
   nInSolution = 0;
   custo = 0;
   residual = I->C[0];

   // first, select all variables already fixed in 1.0
   for(i=0;i<nvars;i++){
      var = varlist[i];
      if(SCIPvarGetLbLocal(var) > 1.0 - EPSILON){ // var >= 1.0
        solution[nInSolution++]=var;
        // update residual capacity
        residual -= I->item[i].weight;
        covered[i]=1;
        // update solution value
        custo += I->item[i].value;
        infeasible = residual < 0?1:0;
#ifdef DEBUG_GRASP
        printf("\nSelected fixed var= %s. TotalItems=%d value=%d residual=%d infeasible=%d", SCIPvarGetName(var), nInSolution, custo, residual, infeasible);
#endif
      }
      else{ // discard items fixed in 0.0
        if(SCIPvarGetUbLocal(var) < EPSILON){ // var fixed in 0.0
        }
        else{
          if (i < n){ // include item i in cand list
           cand[nCands++]=i;
          }
        }
      }
   }

   // alocando o vetor de itens que foram colocados na solucao
   itens_solucao = (solucao *) malloc(sizeof(solucao)*n);
   max_iteracoes = m;
   // quantidade de iteracoes = quantidade de mochilas q eu tenho
    for( i = 0; i < max_iteracoes; i++){
      candidatos = (itemType*) malloc(sizeof(itemType)*I->n);  // alocando o vetor de candidatos

      for (int j = 0; j < n; j++){
         // colocando no vetor candidatos todos os itens que tem peso <= capacidade da mochila i E que nao foram colocados ainda na solucao
         if (I->item[j].weight <= I->C[i] && (verifica_solucao(itens_solucao, I->item[j].label, n_solucao) == 0))
         {
            candidatos[n_cand].label = I->item[j].label; // assim?
            candidatos[n_cand].weight = I->item[j].weight;
            candidatos[n_cand].value = I->item[j].value;
            n_cand++;
         }
      }

      // a partir daqui, candidatos possui todos os itens com peso <= capacidade da mochila i

      int  capacidade_atual = I->C[i];

      while(capacidade_atual > 0.0 && n_cand >= 1){
         min_max(candidatos, n_cand, &maximo, &minimo);
         RCL = (itemType*) malloc(sizeof(itemType)*I->n);
         cria_RCL(candidatos, RCL, minimo, maximo, alpha, n_cand, &n_RCL);

         posicao_item_escolhido = numero_aleatorio(n_RCL);
         rotulo = RCL[posicao_item_escolhido].label;

         // atualzia o custo e diminui a capacidade atual da mochila
         custo += RCL[posicao_item_escolhido].value;
         capacidade_atual -= RCL[posicao_item_escolhido].weight;

         itens_solucao[n_solucao].label = RCL[posicao_item_escolhido].label;
         itens_solucao[n_solucao].mochila = i+1;
         n_solucao++;

         atualiza_candidatos(candidatos, rotulo, &n_cand, capacidade_atual);
   
         free(RCL);
      }

      free(candidatos);

    }

   // complete solution using items not fixed (not covered)
   while(nCands > 0 && residual>0){
      s = randomIntegerB (0, nCands-1);
      selected = cand[s]; // selected candidate
      cand[s] = cand[--nCands]; // remove selected candidate
      // only accept the item if not covered yet and not exceed the capacity
      if(!covered[selected] && I->item[selected].weight <= residual){
         // compute the real value
         valor = I->item[selected].value;
         var = varlist[selected];
         // include selected var in the solution
         solution[nInSolution++]=var;
         // update residual capacity
         residual -= I->item[selected].weight;
         // update covered
         covered[selected] = 1;
         // update the solution value
         custo += valor;
         infeasible = residual<0?1:0;
#ifdef DEBUG_GRASP
         printf("\n\nSelected var= %s. TotalItems=%d value item=%d value = %d residual=%d infeasible=%d\n", SCIPvarGetName(var), nInSolution, valor, custo, residual, infeasible);
#endif
      }
      else{
        // desconsidere o item
#ifdef DEBUG_GRASP
        printf("\n\nNOT selected var= %s. TotalItems=%d value item=%d value=%d residual=%d infeasible=%d\n", SCIPvarGetName(varlist[selected]), nInSolution, valor, custo, residual, infeasible);
#endif
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
#ifdef DEBUG_GRASP
      printf("\nFound solution...\n");
      //      SCIP_CALL( SCIPprintSol(scip, *sol, NULL, FALSE) );
      printf("\ninfeasible=%d value = %lf > bestUb = %lf? %d\n\n", infeasible, valor, bestUb, valor > bestUb + EPSILON);
#endif
      if(!infeasible && custo > bestUb + EPSILON){
#ifdef DEBUG_GRASP
         printf("\nBest solution found...\n");
         SCIP_CALL( SCIPprintSol(scip, *sol, NULL, FALSE) );
#endif
         
         /* verificar se a solucao eh viavel e armazena */
         SCIP_CALL( SCIPtrySolMine(scip, *sol, TRUE, TRUE, FALSE, TRUE, &stored) );
         if( stored )
         {
#ifdef DEBUG_PRIMAL
            printf("\nSolution is feasible and was saved! Total of items = %d", nInSolution);
            SCIPdebugMessage("found feasible GRASP solution:\n");
            SCIP_CALL( SCIPprintSol(scip, *sol, NULL, FALSE) );
#endif
            found = 1;
         }
         else{
           found = 0;
#ifdef DEBUG_GRASP
           printf("\nCould not found\n. BestUb=%lf", bestUb);
#endif
         }
      }
   }
#ifdef DEBUG_GRASP
   getchar();
#endif
   free(varlist);
   free(solution);
   free(covered);
   free(cand);
   return found;


   // =============== FIM DA FUNCAO DA HEURISTICA
}

/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecGrasp)
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

   /* solve grasp */
   if(grasp(scip, &sol, heur)){
     *result = SCIP_FOUNDSOL;
   }
   else{
#ifdef DEBUG_PRIMAL
     printf("\nGrasp could not find feasible solution!");      
#endif
   }
   return SCIP_OKAY;
}

/*
 * primal heuristic specific interface methods
 */

/** creates the grasp_crtp primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurGrasp(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_HEUR* heur;

   /* create grasp primal heuristic data */
   heurdata = NULL;

   heur = NULL;

   /* include primal heuristic */
#if 0
   /* use SCIPincludeHeur() if you want to set all callbacks explicitly and realize (by getting compiler errors) when
    * new callbacks are added in future SCIP versions
    */
   SCIP_CALL( SCIPincludeHeur(scip, HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, param.heur_freq, param.heur_freqofs,
         param.heur_maxdepth, HEUR_TIMING, HEUR_USESSUBSCIP,
         heurCopyGrasp, heurFreeGrasp, heurInitGrasp, heurExitGrasp, heurInitsolGrasp, heurExitsolGrasp, heurExecGrasp,
         heurdata) );
#else
   /* use SCIPincludeHeurBasic() plus setter functions if you want to set callbacks one-by-one and your code should
    * compile independent of new callbacks being added in future SCIP versions
    */
   SCIP_CALL( SCIPincludeHeurBasic(scip, &heur,
         HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, param.heur_round_freq, param.heur_round_freqofs,
         param.heur_round_maxdepth, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecGrasp, heurdata) );

   assert(heur != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetHeurCopy(scip, heur, heurCopyGrasp) );
   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreeGrasp) );
   SCIP_CALL( SCIPsetHeurInit(scip, heur, heurInitGrasp) );
   SCIP_CALL( SCIPsetHeurExit(scip, heur, heurExitGrasp) );
   SCIP_CALL( SCIPsetHeurInitsol(scip, heur, heurInitsolGrasp) );
   SCIP_CALL( SCIPsetHeurExitsol(scip, heur, heurExitsolGrasp) );
#endif

   /* add grasp primal heuristic parameters */
   /* TODO: (optional) add primal heuristic specific parameters with SCIPaddTypeParam() here */

   return SCIP_OKAY;
}

void min_max(itemType *candidatos, int n, int *max, int *min)
{
   *max = candidatos[0].value;
   *min = candidatos[0].value;

   for (int i = 1; i < n; i++)
   {
      if (candidatos[i].value >= *max)
      {
         *max = candidatos[i].value;
      }

      if (candidatos[i].value <= *min)
      {
         
         *min = candidatos[i].value;
      }
   }

}

void cria_RCL(itemType *candidatos, itemType *RCL, int minimo, int maximo, int alpha, int n, int *n_RCL)
{

   for (int i = 0; i < n; i++)
   {
      if (candidatos[i].value >= minimo + alpha * (maximo - minimo))
      {
         RCL[*n_RCL].label = candidatos[i].label;
         RCL[*n_RCL].weight = candidatos[i].weight;
         RCL[*n_RCL].value = candidatos[i].value;
         (*n_RCL) += 1;
      }
   }
}

int numero_aleatorio(int n_cand){
   return rand() % n_cand;
}

void atualiza_candidatos(itemType *candidatos, int rotulo, int *n_cand, int capacidade_atual)
{

   for(int i = 0; i < *n_cand; i++){
      if(candidatos[i].label == rotulo){
         candidatos[i] = candidatos[--(*n_cand)];
      }
   }

   for (int i = 0; i < *n_cand; i++)
   {

      // se isso acontecer, eu preciso remover esse item
      if (candidatos[i].weight > capacidade_atual)
      {

         while (candidatos[(*n_cand) - 1].weight > capacidade_atual && *n_cand > i)
         {
            (*n_cand) -= 1;
         }

         if (*n_cand > i)
         {
            candidatos[i] = candidatos[--(*n_cand)];
         }
      }
   }

}

int verifica_solucao(solucao *itens_solucao, int label, int n){
   if(n == 0){
      return 0;
   }else{
      for(int i = 0; i < n; i++){
         if(itens_solucao[i].label == label){
            return 1;
         }
      }
   }
   return 0; 
}