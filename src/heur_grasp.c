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
//#include "problem.h"   // eu que inclui isso.

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

void min_max(int *candidatos, instanceT* I, int n, int *max, int *min);

void cria_RCL(instanceT* I,int *candidatos, int *RCL, int minimo, int maximo, int alpha, int n, int *n_RCL);

static int numero_aleatorio(int n_cand);

void atualiza_candidatos(instanceT* I, int *candidatos, int rotulo, int *n_cand, int capacidade_atual);

int verifica_solucao(int *covered, int label);

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
   int found = 0, infeasible = 0, nInSolution = 0;
   unsigned int stored;
   int nvars;
   int *covered, n, m, custo = 0, **cand, *nCands;
   SCIP_VAR *var, **solution, **varlist;
   //  SCIP* scip_cp;
   SCIP_Real valor, bestUb;
   SCIP_PROBDATA* probdata;
   int i, *residual;
   int *RCL;
   instanceT* I;

   int max_iteracoes, maximo, minimo, alpha = 0.7;
   int posicao_item_escolhido, n_RCL = 0;
   
   
#ifdef DEBUG_GRASP
   printf("\n============== New grasp heur at node: %lld\n", SCIPnodeGetNumber(SCIPgetCurrentNode(scip)));
#endif

   /* recupera os dados do problema original*/
   probdata=SCIPgetProbData(scip);
   assert(probdata != NULL);

   nvars = SCIPprobdataGetNVars(probdata);
   varlist = SCIPprobdataGetVars(probdata);
   I = SCIPprobdataGetInstance(probdata);
   n = I->n;
   m = I->m;
   int flag = 0;
    
   solution = (SCIP_VAR**) malloc(sizeof(SCIP_VAR*)*n);
   covered = (int*) calloc(n,sizeof(int));
   cand = (int**) malloc(m*sizeof(int*));  // sao m mochilas, ent eh uma lista de candidatos para cada uma

   for(i = 0; i < m; i++){
      cand[i] = (int *) malloc(n*sizeof(int));
   }

   residual = (int*) malloc(m*sizeof(int));
   nCands = (int *) calloc(m,sizeof(int));  // agr eu preciso diferenciar o numero de candidatos para cada mochila

   for(int k = 0; k < m; k++){
      residual[k] = I->C[k];
   }

   // first, select all variables already fixed in 1.0
   // esse for aqui percorre todas as n*m variaveis do problema, certo?
   for(i=0;i<nvars;i++){

      var = varlist[i];  // o que isso aqui faz exatamente??? pega as variaveis da linha i do sistema de equacoes?

      if(SCIPvarGetLbLocal(var) > 1.0 - EPSILON){
         solution[nInSolution++] = var;
         // atualiza a capacidade atual da mochila na qual o item foi colocado
         residual[i%m] -= I->item[i].weight;

         // marca aquele item como coberto
         covered[i%n] = 1;

         // atualiza o custo
         custo += I->item[i].value;
         flag++;
         infeasible = residual < 0?1:0;
#ifdef DEBUG_GRASP
         printf("\nSelected fixed var= %s. TotalItems=%d value=%d residual=%d infeasible=%d", SCIPvarGetName(var), nInSolution, custo, residual, infeasible);
 #endif
      }
      else{// discard items fixed in 0.0
         if(SCIPvarGetUbLocal(var) < EPSILON){ // var fixed in 0.0
         }
         else{
            // aqui ele (item ou variavel?) esta fracionario??
           if (i < n){ // include item i in cand list
            for(int k = 0; k < m; k++){
               // verificando se o item i cabe na mochila k
               if(I->item[i].weight <= residual[k]){
                  cand[k][nCands[k]] = i;
                  nCands[k] += 1;
               }
           }
         }
      }
   }
}

   max_iteracoes = m;
   // quantidade de iteracoes = quantidade de mochilas q eu tenho
   for( i = 0; i < max_iteracoes; i++){

      for (int j = 0; j < n; j++){
         // colocando no vetor candidatos todos os itens que tem peso <= capacidade da mochila i E que nao foram colocados ainda na solucao
         if (I->item[j].weight <= I->C[i] && (verifica_solucao(covered, I->item[j].label) == 0))
         {
            cand[i][nCands[i]] = I->item[j].label;
            /*candidatos[n_cand].label = I->item[j].label; // assim?
            candidatos[n_cand].weight = I->item[j].weight;
            candidatos[n_cand].value = I->item[j].value;
            n_cand++;*/
            nCands[i] += 1;
         }
      }

      // a partir daqui, candidatos possui todos os itens com peso <= capacidade da mochila i

      int capacidade_atual = I->C[i];
      int j = 0;
      while(capacidade_atual > 0.0 && nCands[i] >= 1){
         // passa a lista de candidatos da mochila i, o numero de candidatos que ela tem
         min_max(cand[i],I, nCands[i], &maximo, &minimo);
         RCL = (int*) malloc(sizeof(int)*nCands[i]);
         n_RCL = 0;
         cria_RCL(I, cand[i], RCL, minimo, maximo, alpha, nCands[i], &n_RCL);

         posicao_item_escolhido = numero_aleatorio(n_RCL); // 0 ate n-1

         // atualzia o custo e diminui a capacidade atual da mochila
         valor = I->item[RCL[posicao_item_escolhido]].value;
         custo += valor;
         capacidade_atual -= I->item[RCL[posicao_item_escolhido]].weight;

         var = varlist[posicao_item_escolhido]; 
         solution[nInSolution++] = var;
         covered[posicao_item_escolhido] = 1;

         infeasible = capacidade_atual <0?1:0;
         

         atualiza_candidatos(I, cand[i], posicao_item_escolhido, &nCands[i], capacidade_atual);
         
         // se ao final, a RCL nao estiver vazia, eu "reseto" ela novamente para a proxima iteracao
         if(n_RCL != 0){
            free(RCL);
            RCL = NULL;
         }
         
         
         
         j++;
      }

      free(cand[i]);

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

void min_max(int *candidatos, instanceT* I , int n, int *max, int *min)
{
   *max = I->item[candidatos[0]].value;
   *min = I->item[candidatos[0]].value;

   for (int i = 1; i < n; i++)
   {
      if (I->item[candidatos[i]].value >= *max)
      {
         *max = I->item[candidatos[i]].value;
      }

      if (I->item[candidatos[i]].value <= *min)
      {
         
         *min = I->item[candidatos[i]].value;
      }
   }

}

void cria_RCL(instanceT* I, int *candidatos, int *RCL, int minimo, int maximo, int alpha, int n, int *n_RCL)
{

   for (int i = 0; i < n; i++)
   {
      if (I->item[candidatos[i]].value >= minimo + alpha * (maximo - minimo))
      {
         RCL[*n_RCL] = I->item[candidatos[i]].label-1;
         (*n_RCL) += 1;
      }
   }
}

int numero_aleatorio(int n_cand){
   return rand() % n_cand;
}

void atualiza_candidatos(instanceT* I, int *candidatos, int rotulo, int *n_cand, int capacidade_atual)
{

   for(int i = 0; i < *n_cand; i++){
      if(candidatos[i] == rotulo-1){
         candidatos[i] = candidatos[--(*n_cand)];
      }
   }

   for (int i = 0; i < *n_cand; i++)
   {

      // se isso acontecer, eu preciso remover esse item
      if (I->item[candidatos[i]].weight > capacidade_atual)
      {

         while (I->item[candidatos[(*n_cand) - 1]].weight > capacidade_atual && *n_cand > i)
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

int verifica_solucao(int *covered, int label){

   if(covered[label] == 1){
      // item ja foi colocado em alguma mochila
      return 1;
   }else{
      // esse item n foi colocado em nenhuma mochila
      return 0;
   }

}