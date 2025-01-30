/**@file   problem.c
 * @brief  This file contains routines specific for the problem and the functions loadInstance(), freeInstance, 
 * printInstance, and loadProblem must be implemented 
 *
 **/ 
#include<stdio.h>
#include<math.h>
#include "scip/scip.h"
#include "problem.h"
#include "probdata_mochila.h"

void freeInstance(instanceT* I)
{
  if(I){
    free(I->item);
    free(I->C);
    free(I);
    I = NULL;
  }
}

// criando a instancia = alocando as estruturas
void createInstance(instanceT** I, int n, int m)
{
  *I = (instanceT*) malloc(sizeof(instanceT));
  (*I)->item = (itemType*) malloc(sizeof(itemType)*n);
  (*I)->C = (int*) malloc(sizeof(int)*m);
  (*I)->n = n;
  (*I)->m = m;
}
void printInstance(instanceT* I)
{
  int i;
  printf("\nInstance with n=%d items, m=%d knapsacks", I->n, I->m);
  for(i=0;i<I->n;i++){
     printf("\nKnapsack %d Capacity=%d", i+1, I->C[i]);
  }
  printf("\nItems= \n");
  for(i=0;i<I->n;i++){
     printf("%d value=%d weight=%d\n", I->item[i].label, I->item[i].value, I->item[i].weight);
  }
}
int loadInstance(char* filename, instanceT** I)
{
  FILE* fin;
  int n, m, i;
  fin = fopen(filename, "r");
  if(!fin){
    printf("\nProblem to open file %s\n", filename);
    return 0;
  }
  fscanf(fin,"%d %d\n", &n, &m);

  createInstance(I, n, m);  // criando a instancia

  for(i=0; i<m; i++){
    // lendo as capacidades de cada mochila
     fscanf(fin, "%d\n", &((*I)->C[i]));
  }

  for(i=0; i<n; i++){
    // lendo o nome, o valor e o peso de cada item
     fscanf(fin, "%d %d %d\n", &((*I)->item[i].label), &((*I)->item[i].weight), &((*I)->item[i].value));
  }
  fclose(fin);
  return 1;
}

// load instance problem into SCIP
int loadProblem(SCIP* scip, char* probname, instanceT* I)
{
  SCIP_RETCODE ret_code;

  ret_code = SCIPprobdataCreate(scip, probname, I);
  if(ret_code!=SCIP_OKAY)
    return 0;
  return 1;
}
