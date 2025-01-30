#ifndef __PROBLEM__
#define __PROBLEM__
#include<stdio.h>
#include "scip/scip.h"

/** structure for each item */
typedef struct{
  int label; //  rotulo do item
  int value;  // valor do item
  int weight;  // peso do item
}itemType;

typedef struct{
   int n;   // quant de itens
   int m;   // quant de mochilas
   int *C;   // capacidade de cada mochila. C[0] = capacidade da mochila 0
   itemType *item; /**< data for each item in 0..n-1 */
} instanceT;

void freeInstance(instanceT* I);
void createInstance(instanceT** I, int n, int m);
void printInstance(instanceT* I);
// load instance from a file
int loadInstance(char* filename, instanceT** I);
// load instance problem into SCIP
int loadProblem(SCIP* scip, char* probname, instanceT* in);
#endif
