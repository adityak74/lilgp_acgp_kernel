// Provides prototypes for acgp_mwh.c 

//#ifndef _ACGP_MWH_H
//#define _ACGP_MWH_H

// #define DEBUG_LVL1

#define MUTATION_SOURCE 0
#define REGROW_SOURCE 1
#define WRAPPER_SOURCE 2

#include "lilgp.h"

#define L2minweight 0.0000001

#define SMALL 0.000000001

// typedef struct typeHyperCube typelevel2info;

typedef struct _typeHyperCube
{
	int *count;				// hypercube of dimension Dim implemented as 1dim array
	double *prob;
	double *wheel;
	int wheel_calculated;
	int dim;					// dimension of our Hypercube
	int *offset;			// offset of the index to each array
	int *size;         // size of each dimension      
	int num_elems;			// Total elements in our hypercube
	int count_tot;	
	
	struct _typeHyperCube *L2;
	
}	typeHyperCube;

typedef struct
{
	typeHyperCube *f;	// array of hypercubes, of length fset
	int total;		// number of hypercubes

} typeLevel1_Cnt;

extern typeLevel1_Cnt Level1_Cnt_mwh;

extern typeLevel1_Cnt Level2_Cnt_mwh;

extern int child_list[10];

// How deep to store root information, default of 1 if using
// extra information operators
// int acgp_level;  

extern int *level1context_mwh;
extern int *level2context_mwh;
extern int size_context1_mwh;
extern int size_context2_mwh;
extern int depth_mwh;

int root_fset_index_to_counter_index(int trueIndex);

void inc_counter_with_context_mwh(int trueIndex,int expressed_count);
void update_context_mwh(int depth,int index);

void reset_counters_mwh(void);

void print_cnts_mwh(int curGen,const char *basefile,int flag);
void print_weights_mwh(int curGen,const char *basefile, int flag);
void adjust_weights_with_counters_mwh(double prctChnge, double thresholdPrct);
void gather_info_mwh(lnode *data,lnode *target,int *replace_depth,int *replace_arity,int *trueIndex,int *context);

/* Wrapper functions for generate_random_*_tree functions to take into account lvl 1 information
   Basically they generate the lvl1 information and then call the normal generate_random_*_tree functions
   to generate past lvl 1 */

void create_mwh(void);
			/* Will allocate and initialize our level 1 counters by accessing
				global fset function table */

// #endif
