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

/**@file   heur_aleatoria.h
 * @ingroup PRIMALHEURISTICS
 * @brief  aleatoria primal heuristic
 * @author Edna Hoshino (based on template provided by Tobias Achterberg)
 *
 * template file for primal heuristic plugins
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_HEUR_ALEATORIA_H__
#define __SCIP_HEUR_ALEATORIA_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif
  
//int printVarLpSol(SCIP* scip);
//int incluiCandidatoB(candidatoT* candidatos, int total, int *nCands, int c, SCIP_Real solval);
//int randomIntegerB (int low, int high);
//int updateSolution(SCIP_VAR* var, int n, int* covered, int *nCovered, int *custo);
//SCIP_Real createSolution(SCIP* scip, SCIP_SOL* sol, SCIP_VAR** solution, int nSolution, int *infeasible);
//SCIP_RETCODE selectCand(SCIP* scip, SCIP_VAR** solution, int nInSolution, int custo, SCIP_VAR** pvar, SCIP_VAR** varlist, int n1, int nfracimplvars, int* covered);
//SCIP_RETCODE selectCand(SCIP* scip, int totalRings, int custo, SCIP_VAR** pvar, SCIP_VAR** varlist, int n1, int nfracimplvars, int* covered);
//SCIP_VAR* selectCand(SCIP* scip, SCIP_VAR* tabu, int totalRings, int custo, int lucro);
//SCIP_RETCODE updateMasterAndPricingProblem(SCIP* scip, SCIP* pricer, SCIP_VAR* var);
int aleatoria(SCIP* scip, SCIP_SOL** sol, SCIP_HEUR* heur);

/** creates the aleatoria primal heuristic and includes it in SCIP */

SCIP_RETCODE SCIPincludeHeurAleatoria(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
