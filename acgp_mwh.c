/* acgp_mwh.c
	File to keep track of complete depth information up to level 1
	and complete local information of depth 1
	
	Written by Mark Hauschild with extensive help by Dr. Cezary Janikow
	
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "lilgp.h"

#define MINWGHT 0.00100


int this_iteration = -1;

// sets number of times to try and regrow terminals when needing to end
// If we cannot find a set of terminals after this number of tries, exit
const int num_retries = 100;

typeLevel1_Cnt Level1_Cnt_mwh;  // global count datastructure to hold root level counters
typeLevel1_Cnt Level1_Cnt_f_mwh; // wheel for functions only
typeLevel1_Cnt Level1_Cnt_t_mwh;  // wheel for terminals only

typeLevel1_Cnt Level2_Cnt_mwh;

typeLevel1_Cnt Local1_mwh;      // global data structure to hold local counters
typeLevel1_Cnt Local1_t_mwh; // wheel for terminals only
typeLevel1_Cnt Local1_f_mwh; // wheel for functions only

int *level1context_mwh;	// list of the indices of lvl 1 context data
int *level2context_mwh;  // list of indices of lvl 2 context data
int size_context1_mwh;   // number of elements in level1context_mwh
int size_context2_mwh;   // number of elements in level2context_mwh

int child_list[10];  // parameter returned by spin_level1_information
int depth_mwh;

void pull_L2_probs_from_MS_czj(void);

void convert_int_to_indices(int *index,int i,int dim, int *offset) 
// Converts a onedimensional array index 'i' to a multidimensional array index 'index'
// Needs an array of offsets and the dimension of the offset array
{
	int remain = i;
	int d;
	for (d = dim - 1; d >=0; d--) {
		index[d] = remain / offset[d];
		remain = remain - index[d]*offset[d];
	}

}

int root_fset_index_to_counter_index(int trueIndex) 
// Takes a integer that represents a function fset, that is at the root
// and converts it to an integer representing this root in our level1counters
{
	int f = 0;
	int NumFT = fset[0].function_count;
	
	// convert trueIndex(in fset) to an index to our counters
	f = MS_czj[NumFT][0].convert_to_ms[trueIndex];
	return f;
}

void generate_wheels_lvl1_mwh(int f)
// Generates a regrow wheel at lvl1
{
	if (Level1_Cnt_mwh.f[f].wheel_calculated == 1)
		return;
		
	Level1_Cnt_mwh.f[f].wheel_calculated = 1;
	if (Level1_Cnt_mwh.f[f].wheel == NULL)
		Level1_Cnt_mwh.f[f].wheel 
		= (double *) calloc(Level1_Cnt_mwh.f[f].num_elems,sizeof(double));
	
	// if our prob is below min weight, then set the space on wheel to 0
	if (Level1_Cnt_mwh.f[f].prob[0] > L2minweight + SMALL)
		Level1_Cnt_mwh.f[f].wheel[0] = Level1_Cnt_mwh.f[f].prob[0];
	else
		Level1_Cnt_mwh.f[f].wheel[0] = 0;
	
	int c;
	for (c = 0; c < Level1_Cnt_mwh.f[f].num_elems-1; c++) {
		
		// Take into account that we might want to have 0 space on this wheel
		if (Level1_Cnt_mwh.f[f].prob[c+1] > L2minweight + SMALL) {
			Level1_Cnt_mwh.f[f].wheel[c+1] = Level1_Cnt_mwh.f[f].wheel[c] 
				+ Level1_Cnt_mwh.f[f].prob[c+1];
		}
		else {
			Level1_Cnt_mwh.f[f].wheel[c+1] = Level1_Cnt_mwh.f[f].wheel[c] 
				+ 0;
		}
	}
	
	// Terminals wheel code
	Level1_Cnt_t_mwh.f[f].wheel_calculated = 1;
	if (Level1_Cnt_t_mwh.f[f].wheel == NULL)
		Level1_Cnt_t_mwh.f[f].wheel = (double*) calloc(Level1_Cnt_t_mwh.f[f].num_elems,sizeof(double));
	
	 
	int NumF = fset[0].function_count;
	
	// First generate our probabilities for terminal wheel
	int *index = (int *) calloc(Level1_Cnt_t_mwh.f[f].dim,sizeof(int));
	for (c=0; c < Level1_Cnt_t_mwh.f[f].num_elems; c++) {
		// first convert this to index to our terminals array
		convert_int_to_indices(index,c,Level1_Cnt_t_mwh.f[f].dim,Level1_Cnt_t_mwh.f[f].offset);
		
		// now upsize this index to an index into the entire array
		int d;
		for (d = 0; d < Level1_Cnt_t_mwh.f[f].dim; d++) {
			index[d] += MS_czj[NumF][0].numF;
		}
		
		// Now convert 'index' to an index into our full array
		int index_to_full = 0;
		int i;
		for (i = 0; i < Level1_Cnt_t_mwh.f[f].dim; i++) {
			index_to_full += index[i]*Level1_Cnt_mwh.f[f].offset[i];
		}
		
		// Now set our weight in the terminal array
		Level1_Cnt_t_mwh.f[f].wheel[c] = Level1_Cnt_mwh.f[f].prob[index_to_full];	
	} 
	free(index);
	
	// Take care of having below minimum weight, so we do not want any chance to spin it
	if (Level1_Cnt_t_mwh.f[f].wheel[0] <= L2minweight + SMALL)
		Level1_Cnt_t_mwh.f[f].wheel[0] = 0;
	
	// Then turn the probabilities into a wheel
	for (c = 1; c < Level1_Cnt_t_mwh.f[f].num_elems; c++) {
		// Make sure to account for being below minweight, then allocate 0 space for it on our wheel
		if (Level1_Cnt_t_mwh.f[f].wheel[c] <= L2minweight + SMALL) {
			Level1_Cnt_t_mwh.f[f].wheel[c] = Level1_Cnt_t_mwh.f[f].wheel[c-1];
		}
		else {
			Level1_Cnt_t_mwh.f[f].wheel[c] = Level1_Cnt_t_mwh.f[f].wheel[c] + Level1_Cnt_t_mwh.f[f].wheel[c-1];
		}
	}
	
	
	// Functions wheel code
	Level1_Cnt_f_mwh.f[f].wheel_calculated = 1;
	if (Level1_Cnt_f_mwh.f[f].wheel == NULL)
		Level1_Cnt_f_mwh.f[f].wheel = (double*) calloc(Level1_Cnt_f_mwh.f[f].num_elems,sizeof(double));
	
	// First generate our probabilities for function wheel
	index = (int *) calloc(Level1_Cnt_f_mwh.f[f].dim,sizeof(int));
	for (c=0; c < Level1_Cnt_f_mwh.f[f].num_elems; c++) {
		// first convert this to index to our terminals array
		convert_int_to_indices(index,c,Level1_Cnt_f_mwh.f[f].dim,Level1_Cnt_f_mwh.f[f].offset);
				
		// Now convert 'index' to an index into our full array
		int index_to_full = 0;
		int i;
		for (i = 0; i < Level1_Cnt_f_mwh.f[f].dim; i++) {
			index_to_full += index[i]*Level1_Cnt_mwh.f[f].offset[i];
		}
		
		// Now set our weight in the terminal array
		Level1_Cnt_f_mwh.f[f].wheel[c] = Level1_Cnt_mwh.f[f].prob[index_to_full];	
	
	} 
	free(index);	
	
	// Take care of having below minimum weight, so we do not want any chance to spin it
	if (Level1_Cnt_f_mwh.f[f].wheel[0] <= L2minweight + SMALL)
		Level1_Cnt_f_mwh.f[f].wheel[0] = 0;
	
	// Then turn the probabilities into a wheel
	for (c = 1; c < Level1_Cnt_f_mwh.f[f].num_elems; c++) {
		// Make sure to account for being below minweight, then allocate 0 space for it on our wheel
		if (Level1_Cnt_f_mwh.f[f].wheel[c] <= L2minweight + SMALL) {
			Level1_Cnt_f_mwh.f[f].wheel[c] = Level1_Cnt_f_mwh.f[f].wheel[c-1];
		}
		else {
			Level1_Cnt_f_mwh.f[f].wheel[c] = Level1_Cnt_f_mwh.f[f].wheel[c] + Level1_Cnt_f_mwh.f[f].wheel[c-1];
		}
	}
}

void generate_regrow_wheel_lcl1_mwh(int f)
// Generates only one regrow wheel based on the index of the parent 'f'
{
	if (Local1_mwh.f[f].wheel_calculated == 1)
		return;
		
	Local1_mwh.f[f].wheel_calculated = 1;
	if (Local1_mwh.f[f].wheel == NULL)
		Local1_mwh.f[f].wheel 
		= (double *) calloc(Local1_mwh.f[f].num_elems,sizeof(double));
	
	// if our prob is below min weight, then set the space on wheel to 0	
	if ( Local1_mwh.f[f].prob[0] > L2minweight + SMALL)
		Local1_mwh.f[f].wheel[0] = Local1_mwh.f[f].prob[0];
	else
		Local1_mwh.f[f].wheel[0] = 0;
	
	int c;
	for (c = 0; c < Local1_mwh.f[f].num_elems-1; c++) {
		if (Local1_mwh.f[f].prob[c+1] > L2minweight + SMALL) {
			Local1_mwh.f[f].wheel[c+1] = Local1_mwh.f[f].wheel[c] 
				+ Local1_mwh.f[f].prob[c+1];
		}
		else {
			Local1_mwh.f[f].wheel[c+1] = Local1_mwh.f[f].wheel[c]
				+ 0;
		}
	}
	
	// Terminals wheel code
	Local1_t_mwh.f[f].wheel_calculated = 1;
	if (Local1_t_mwh.f[f].wheel == NULL)
		Local1_t_mwh.f[f].wheel
		= (double *) calloc(Local1_t_mwh.f[f].num_elems,sizeof(double));
	
	// We have wheel for our total, now lets generate one for terminals only
	
	// First generate our probabilities
	int *index = (int *) calloc(Local1_t_mwh.f[f].dim,sizeof(int));
	for (c=0; c < Local1_t_mwh.f[f].num_elems; c++) {
		// first convert this to index to our terminals array
		convert_int_to_indices(index,c,Local1_t_mwh.f[f].dim,Local1_t_mwh.f[f].offset);
		
		// now upsize this index to an index into the entire array
		int d;
		for (d = 0; d < Local1_t_mwh.f[f].dim; d++) {
			index[d] += MS_czj[f][d].numF;
		}
		
		// Now convert 'index' to an index into our full array
		int index_to_full = 0;
		int i;
		for (i = 0; i < Local1_t_mwh.f[f].dim; i++) {
			index_to_full += index[i]*Local1_mwh.f[f].offset[i];
		}
		
		// Now set our weight in the terminal array
		Local1_t_mwh.f[f].wheel[c] = Local1_mwh.f[f].prob[index_to_full];	
	} 
	
	// Take care of having below minimum weight, so we do not want any chance to spin it
	if (Local1_t_mwh.f[f].wheel[0] <= L2minweight + SMALL)
		Local1_t_mwh.f[f].wheel[0] = 0;
	
	// Then turn the probabilities into a wheel
	for (c = 1; c < Local1_t_mwh.f[f].num_elems; c++) {
		// Make sure to account for being below minweight, then allocate 0 space for it on our wheel
		if (Local1_t_mwh.f[f].wheel[c] <= L2minweight + SMALL) {
			Local1_t_mwh.f[f].wheel[c] = Local1_t_mwh.f[f].wheel[c-1];
		}
		else {
			Local1_t_mwh.f[f].wheel[c] = Local1_t_mwh.f[f].wheel[c] + Local1_t_mwh.f[f].wheel[c-1];
		}
	}
	
	
	
	free(index);
	// Now generate wheel for functions only 
	
	Local1_f_mwh.f[f].wheel_calculated = 1;
	if (Local1_f_mwh.f[f].wheel == NULL)
		Local1_f_mwh.f[f].wheel
		= (double *) calloc(Local1_f_mwh.f[f].num_elems,sizeof(double));
	
	// First generate our probabilities
	index = (int *) calloc(Local1_f_mwh.f[f].dim,sizeof(int));
	for (c=0; c < Local1_f_mwh.f[f].num_elems; c++) {
		// first convert this to index to our terminals array
		convert_int_to_indices(index,c,Local1_f_mwh.f[f].dim,Local1_f_mwh.f[f].offset);
		
		// Now convert 'index' to an index into our full array
		int index_to_full = 0;
		int i;
		for (i = 0; i < Local1_t_mwh.f[f].dim; i++) {
			index_to_full += index[i]*Local1_mwh.f[f].offset[i];
		}
		
		// Now set our weight in the terminal array
		Local1_f_mwh.f[f].wheel[c] = Local1_mwh.f[f].prob[index_to_full];	
	} 
	
	
	// Take care of having below minimum weight, so we do not want any chance to spin it
	if (Local1_f_mwh.f[f].wheel[0] <= L2minweight + SMALL)
		Local1_f_mwh.f[f].wheel[0] = 0;
	
	// Then turn the probabilities into a wheel
	for (c = 1; c < Local1_f_mwh.f[f].num_elems; c++) {
		// Make sure to account for being below minweight, then allocate 0 space for it on our wheel
		if (Local1_f_mwh.f[f].wheel[c] <= L2minweight + SMALL) {
			Local1_f_mwh.f[f].wheel[c] = Local1_f_mwh.f[f].wheel[c-1];
		}
		else {
			Local1_f_mwh.f[f].wheel[c] = Local1_f_mwh.f[f].wheel[c] + Local1_f_mwh.f[f].wheel[c-1];
		}
	}
	
	free(index);
}

void spin_regrow_lcl1_information_mwh(int f,int *child_lcl1_regrow)
{
		
	// now trueIndex should be the index to our cube for parIndex
	if (Local1_mwh.f[f].wheel_calculated == 0)
		generate_regrow_wheel_lcl1_mwh(f);
	
#ifdef DEBUG_LCL1
	int j;
	printf("Local regrow wheel inside nonterminal lcl1 spin for %d is\n",f);
	printf("[");
	for (j = 0; j < Local1_mwh.f[f].num_elems; j++) 
	{
		printf("%f,",Local1_mwh.f[f].wheel[j]);
	}
	printf("]\n");
	if (Local1_mwh.f[f].wheel[0] == 0)
	{
		printf("In spin regrow lcl1 nonterminals but wheel not calculated\n");
		exit(1);
	}
#endif	
	
	
	double spin = random_double()
	   * Local1_mwh.f[f].wheel[Local1_mwh.f[f].num_elems-1];
	int winner = 0;
	//printf("Actual local spin is %f,",spin);
	while (spin > Local1_mwh.f[f].wheel[winner])
		winner++;
	
	//printf(" internal winner index was %d\n",winner); fflush(stdout);
	
	// Now lets find out the index for this
	int i; 
	convert_int_to_indices(child_lcl1_regrow,winner,Local1_mwh.f[f].dim,Local1_mwh.f[f].offset);
	
	// Now convert it to fset indices
	int c;
	for (c =0; c < (fset->cset)[f].arity; c++) {
		child_lcl1_regrow[c] = MS_czj[f][c].members[child_lcl1_regrow[c]];
	}	
	
}

void spin_regrow_lcl1_information_terminals_mwh(int f,int *child_lcl1_regrow)
{
		
	// now trueIndex should be the index to our cube for parIndex
	if (Local1_mwh.f[f].wheel_calculated == 0) {
		//printf("Generating regrow wheel in spin-regrow_lcl1_terminals_mwhw\n");
		generate_regrow_wheel_lcl1_mwh(f);
	}
	
#ifdef DEBUG_LCL1
	int j;
	printf("Local regrow wheel inside terminal lcl1 spin for %d is\n",f);
	printf("[");
	for (j = 0; j < Local1_t_mwh.f[f].num_elems; j++) 
	{
		printf("%f,",Local1_t_mwh.f[f].wheel[j]);
	}
	printf("]\n");
	
#endif	
	// Check if we can even generate all terminals for this one
	if (Local1_t_mwh.f[f].num_elems <= 0 || Local1_t_mwh.f[f].wheel[ Local1_t_mwh.f[f].num_elems - 1] <= L2minweight + SMALL) {
		//// If not, lets default by calling the more general one to allow anything
		//spin_regrow_lcl1_information_mwh(f,child_lcl1_regrow);
		
		// Lets use lvl 1 cgp heuristics
		Function_czj = f;
		int d;
		for (d = 0; d < Local1_t_mwh.f[f].dim; d++) {
			Argument_czj = d;
			child_lcl1_regrow[d]  = random_T_czj(); 
		}
		
		return;
	}
	
	double spin = random_double() * (Local1_t_mwh.f[f].wheel[Local1_t_mwh.f[f].num_elems-1]);
	int winner = 0;
	//printf("Actual local spin is %f,",spin);
	while (spin > Local1_t_mwh.f[f].wheel[winner])
		winner++;
	
	//printf(" internal winner index was %d\n",winner); fflush(stdout);
	
	// Now lets find out the index for this
	int i; 
	convert_int_to_indices(child_lcl1_regrow,winner,Local1_t_mwh.f[f].dim,Local1_t_mwh.f[f].offset);
	
	// now this index is only index to terminal wheel, need to upsize it to the full one
	int d;
	for (d = 0; d < Local1_t_mwh.f[f].dim; d++) {
		child_lcl1_regrow[d] += MS_czj[f][d].numF;
	}
	
	// Now convert it to fset indices
	int c;
	for (c =0; c < (fset->cset)[f].arity; c++) {
		child_lcl1_regrow[c] = MS_czj[f][c].members[child_lcl1_regrow[c]];
	}	
	
}


void spin_regrow_lcl1_information_functions_mwh(int f,int *child_lcl1_regrow)
{
		
	// now trueIndex should be the index to our cube for parIndex
	if (Local1_mwh.f[f].wheel_calculated == 0)
		generate_regrow_wheel_lcl1_mwh(f);
	
#ifdef DEBUG_LCL1
	int j;
	printf("Local regrow wheel inside function lcl1 spin for %d is\n",f);
	printf("[");
	for (j = 0; j < Local1_f_mwh.f[f].num_elems; j++) 
	{
		printf("%f,",Local1_f_mwh.f[f].wheel[j]);
	}
	printf("]\n");
	
#endif	
	
	// Check if we can even generate all functions for this one
	if (Local1_f_mwh.f[f].num_elems <= 0 || Local1_f_mwh.f[f].wheel[ Local1_f_mwh.f[f].num_elems - 1] <= L2minweight + SMALL) {
		//// If not, lets default by calling the more general one to allow anything
		//spin_regrow_lcl1_information_mwh(f,child_lcl1_regrow);
		//return;
		
		// Lets use lvl 1 cgp heuristics
		Function_czj = f;
		int d;
		for (d = 0; d < Local1_t_mwh.f[f].dim; d++) {
			Argument_czj = d;
			child_lcl1_regrow[d]  = random_F_czj(); 
		}
		
		return;
	}
	
	
	double spin = random_double()
	   * Local1_f_mwh.f[f].wheel[Local1_f_mwh.f[f].num_elems-1];
	int winner = 0;

	while (spin > Local1_f_mwh.f[f].wheel[winner])
		winner++;
	
	// Now lets find out the index for this
	int i; 
	convert_int_to_indices(child_lcl1_regrow,winner,Local1_f_mwh.f[f].dim,Local1_f_mwh.f[f].offset);

	// Now convert it to fset indices
	int c;
	for (c =0; c < (fset->cset)[f].arity; c++) {
		child_lcl1_regrow[c] = MS_czj[f][c].members[child_lcl1_regrow[c]];
	}	
	
}



void spin_lvl1_information_mwh(int parIndex)
{
	// trueindex is index to our lvl1 counter, parIndex is raw fset index
	// need to convert
	// int trueIndex = MS_czj[NumFT][0].members[parIndex];
	int NumFT = fset[0].function_count;
	int trueIndex;
	
	trueIndex = MS_czj[NumFT][0].convert_to_ms[parIndex];
	
	// now trueIndex should be the index to our cube for parIndex
	if (Level1_Cnt_mwh.f[trueIndex].wheel_calculated == 0) {
		//printf("Generating wheels in spin_lvl1_information_mwh\n");
		generate_wheels_lvl1_mwh(trueIndex);
	}
	
	//printf("Showing wheel for %d\n[",trueIndex);
	int c;
	//for (c = 0; c < Level1_Cnt_mwh.f[trueIndex].num_elems; c++)
	//	printf("%f ",Level1_Cnt_mwh.f[trueIndex].wheel[c]);
	//printf("]\n"); 
	
	double spin = random_double()
	   * Level1_Cnt_mwh.f[trueIndex].wheel[Level1_Cnt_mwh.f[trueIndex].num_elems-1];
	int winner = 0;
	
	while (spin > Level1_Cnt_mwh.f[trueIndex].wheel[winner])
		winner++;
	
	// Now lets find out the index for this
	int i; 
	
	convert_int_to_indices(child_list,winner,Level1_Cnt_mwh.f[trueIndex].dim,Level1_Cnt_mwh.f[trueIndex].offset);
		
	// now our child_list is set, but still needs to be converted back to actual indices	
	// Now we have our child_list[c], but it is index to our cube, we need to convert it to indexes to fset
	//int c;
	for (c =0; c < (fset->cset)[parIndex].arity; c++) {
		child_list[c] = MS_czj[parIndex][c].members[child_list[c]];
	}
}

void spin_lvl1_information_terminals_mwh(int parIndex)
{
	// trueindex is index to our lvl1 counter, parIndex is raw fset index
	// need to convert
	// int trueIndex = MS_czj[NumFT][0].members[parIndex];
	int NumFT = fset[0].function_count;
	int trueIndex;
	
	trueIndex = MS_czj[NumFT][0].convert_to_ms[parIndex];
	
	//printf("---INSIDE spin_lvl1_information_terminals_mwh---\n");
	
	// now trueIndex should be the index to our cube for parIndex
	if (Level1_Cnt_mwh.f[trueIndex].wheel_calculated == 0) {
		//printf("Generating wheels in spin_lvl1_information_terminals_mwh\n");
		generate_wheels_lvl1_mwh(trueIndex);
	}
	
	int c;
	
	// Check to see if we have a valid wheel to spin, if not go with L1 heuristics to generate them all
	if (Level1_Cnt_t_mwh.f[trueIndex].num_elems <= 0 || Level1_Cnt_t_mwh.f[trueIndex].wheel[ Level1_Cnt_t_mwh.f[trueIndex].num_elems - 1] <= L2minweight + SMALL) {
		//printf("Cannot spin for a terminal, spinning for all\n");
		//spin_lvl1_information_mwh(parIndex);
		//return;
		
		// Lets use lvl 1 cgp heuristics
		Function_czj = parIndex;
		int d;
		for (d = 0; d < Local1_t_mwh.f[parIndex].dim; d++) {
			Argument_czj = d;
			child_list[d]  = random_T_czj(); 
		}
		
		return;
	}
	
	double spin = random_double()
	   * Level1_Cnt_t_mwh.f[trueIndex].wheel[Level1_Cnt_t_mwh.f[trueIndex].num_elems-1];
	int winner = 0;
	//printf("Actual spin is %f,",spin);
	while (spin > Level1_Cnt_t_mwh.f[trueIndex].wheel[winner])
		winner++;
	
	//printf(" internal winner index was %d\n",winner);
	
	convert_int_to_indices(child_list,winner,Level1_Cnt_t_mwh.f[trueIndex].dim,Level1_Cnt_t_mwh.f[trueIndex].offset);

	// now this index is only index to terminal wheel, need to upsize it to the full one
	int d;
	for (d = 0; d < Level1_Cnt_t_mwh.f[trueIndex].dim; d++) {
		// Should this be += MS_czj[NumFT][d].numF?
		child_list[d] += MS_czj[parIndex][d].numF;
	}
	
	// Now we have our child_list[c], but it is index to our cube, we need to convert it to indexes to fset
	//int c;
	for (c =0; c < (fset->cset)[parIndex].arity; c++) {
		child_list[c] = MS_czj[parIndex][c].members[child_list[c]];
	}
	
}

void spin_lvl1_information_functions_mwh(int parIndex)
{
	// trueindex is index to our lvl1 counter, parIndex is raw fset index
	// need to convert
	// int trueIndex = MS_czj[NumFT][0].members[parIndex];
	int NumFT = fset[0].function_count;
	int trueIndex;
	
	trueIndex = MS_czj[NumFT][0].convert_to_ms[parIndex];
	
	// now trueIndex should be the index to our cube for parIndex
	if (Level1_Cnt_mwh.f[trueIndex].wheel_calculated == 0)
		generate_wheels_lvl1_mwh(trueIndex);
	
	int c;
	
	// Check to see if we have a valid wheel to spin, if not go with L1 heuristics to generate them all
	if (Level1_Cnt_f_mwh.f[trueIndex].num_elems <= 0 || Level1_Cnt_f_mwh.f[trueIndex].wheel[ Level1_Cnt_f_mwh.f[trueIndex].num_elems - 1] <= L2minweight + SMALL) {
		//printf("Could not find a function, spinning all\n");
		//spin_lvl1_information_mwh(parIndex);
		//return;
		
		// Lets use lvl 1 cgp heuristics
		Function_czj = parIndex;
		int d;
		for (d = 0; d < Local1_t_mwh.f[parIndex].dim; d++) {
			Argument_czj = d;
			child_list[d]  = random_F_czj(); 
		}
		
		return;
	}
	
	double spin = random_double()
	   * Level1_Cnt_f_mwh.f[trueIndex].wheel[Level1_Cnt_f_mwh.f[trueIndex].num_elems-1];
	int winner = 0;
	while (spin > Level1_Cnt_f_mwh.f[trueIndex].wheel[winner])
		winner++;
	
	convert_int_to_indices(child_list,winner,Level1_Cnt_f_mwh.f[trueIndex].dim,Level1_Cnt_f_mwh.f[trueIndex].offset);
	
	// now our child_list is set, but still needs to be converted back to actual indices	
	// Now we have our child_list[c], but it is index to our cube, we need to convert it to indexes to fset
	//int c;
	for (c =0; c < (fset->cset)[parIndex].arity; c++) {
		child_list[c] = MS_czj[parIndex][c].members[child_list[c]];
	}
}



void inc_local_level1_mwh(int trueIndex, int *child_list,int arity, int expressed_count) {
	int i,j;
		
	// arity is 0 so no lcl1 info for this
	if (arity == 0)
	   return;
		
	#ifdef DEBUG_LCL1
	printf("Storing local (%d)(",trueIndex);
	for (i = 0; i < arity; i++) {
		printf("%d,",child_list[i]);
	}
	printf(")\n");
	fflush(stdout);
	#endif
	
	// Need to convert child_list values to indices to our counters
	// by using MST_czj
	for (i=0; i < arity; i++) {
		child_list[i] = MS_czj[trueIndex][i].convert_to_ms[ child_list[i] ];
	}
	
	#ifdef DEBUG_LCL1
	printf("After conversion we have local (%d)(",trueIndex);
	for (i = 0; i < arity; i++) {
		printf("%d,",child_list[i]);
		if (child_list[i] < 0) {
			printf("Invalid conversion, have a negative index when trying to incr. counters\n");
			printf("in inc_local_level1_mwh. Probably means an operator is generating an invalid\n");
			printf("tree that does not match constraints.\n");
			exit(1);
		}
	}
	printf(")\n");
	fflush(stdout);
	#endif
	
	// Generate our final index by using offset array
	int num_indexes = Local1_mwh.f[trueIndex].dim;
	int index = 0;
	for (i = 0; i < num_indexes; i++) {
		index += child_list[i]*Local1_mwh.f[trueIndex].offset[i];
	}	
	
	// now increment Parent's hypercube at index by 1
	Local1_mwh.f[trueIndex].count[index] += expressed_count;
	Local1_mwh.f[trueIndex].count_tot += expressed_count;
}

void inc_counter_with_context1_mwh(int trueIndex, int expressed_count) 
	/* Using context information from a tree
		increment the appropriate counter 
		Should be a global variable Parent as well as a
		global list of offsets representing the indexes
		parIndex is index to fset
		we need to convert it to trueIndex, index to our counters first
					*/
{
		
	int i,j;
	// Take care if we are called with a terminal at root
	if (size_context1_mwh == 0)
		return;
	
	// printf("Received %d (index to fset) from acgp_count, need to convert to lvl1counter\n",parIndex);
	// fflush(stdout);
	
	int NumFT = fset[0].function_count;
		
	// convert trueIndex to an index to our mutation set	
	int f = MS_czj[NumFT][0].convert_to_ms[trueIndex];
		
#ifdef DEBUG_LVL1
	printf("Storing lvl1 (%d)(",trueIndex);
	for (i = 0; i < size_context1_mwh; i++) {
		printf("%d,",level1context_mwh[i]);
	}
	printf(")\n");
	fflush(stdout);

	printf("Inside inc_c_with_context, storing Lvl1 of: %s (",fset[0].cset[parIndex].string);
	for (i = 0; i < size_context_mwh; i++)
		printf("%s ",fset[0].cset[level1context_mwh[i]].string);
	printf("),size_context_mwh=%d\nraw lvl1contextinfo=[",size_context_mwh);
	for (i=0; i < size_context_mwh; i++) {
		printf("%d ",level1context_mwh[i]);
	}
	printf("\n");
	fflush(stdout);
#endif
	
	// Need to convert level1context_mwh values to indices to our counters
	// by using MST_czj
	for (i=0; i < size_context1_mwh; i++) {
		level1context_mwh[i] = MS_czj[trueIndex][i].convert_to_ms[level1context_mwh[i]];
	}
	
	// Generate our final index by using offset array
	int num_indexes = Level1_Cnt_mwh.f[f].dim;
	int index = 0;
	for (i = 0; i < num_indexes; i++) {
		index += level1context_mwh[i]*Level1_Cnt_mwh.f[f].offset[i];
	}

#ifdef DEBUG_LVL1
	if (index < 0 || index >= Level1_Cnt_mwh.f[f].num_elems){
		printf("Trying to index a hypercube out of bounds!\n");
		printf("Index was %d, max elements was %d\n",index,Level1_Cnt_mwh.f[f].num_elems);
		
		printf("Printing out level1context_mwh,[");
		for (i = 0; i < num_indexes; i++) {
			printf("%d ",level1context_mwh[i]);
		}
		printf("]\n");
		
		printf("Printing out offsets,[");
		for (i = 0; i < num_indexes; i++) {
			printf("%d ",Level1_Cnt_mwh.f[f].offset[i]);
		}
		printf("]\n");
		
		exit(1);
	}
#endif	
	
	// now increment Parent's hypercube at index by 1
	Level1_Cnt_mwh.f[f].count[index] += expressed_count;
	Level1_Cnt_mwh.f[f].count_tot += expressed_count;

}



void inc_counter_with_context_mwh(int trueIndex, int expressed_count)
{
	inc_counter_with_context1_mwh(trueIndex, expressed_count);
	
	size_context2_mwh = 0;
	size_context1_mwh = 0;	
}

double acgp_new_weight_lvl1_mwh(double oldW, double prctChnge, double statW,int mutSetSize, double thresholdPrct) 
{ double newW;
  newW = oldW*(1-prctChnge)+statW*prctChnge;
  /*if (newW<thresholdPrct/mutSetSize || newW<MINWGHT)
  	return MINWGHT;
	return MINWGHT/mutSetSize;
else
	return newW;*/
   if (newW < 0.000001)
   	return 0.000001;
   return newW;	
	
	
}


void adjust_lvl1_weights_with_counters_mwh(double prctChnge,double thresholdPrct) 
/*
Update all probabilities based on the counts for lvl1 functions
and then update all wheels. Right now always calculates all wheels immediately
*/
{
	int f,c;
	for (f = 0; f < Level1_Cnt_mwh.total; f++) {
		Level1_Cnt_mwh.f[f].wheel_calculated = 0;
		Level1_Cnt_t_mwh.f[f].wheel_calculated = 0;
		Level1_Cnt_f_mwh.f[f].wheel_calculated = 0;
		
	        if (Level1_Cnt_mwh.f[f].count_tot) {
			for (c = 0; c < Level1_Cnt_mwh.f[f].num_elems;c++) {
				Level1_Cnt_mwh.f[f].prob[c] =
					 
			   	acgp_new_weight_lvl1_mwh(
			        	Level1_Cnt_mwh.f[f].prob[c],
			   		prctChnge,
					(double)((double)Level1_Cnt_mwh.f[f].count[c]/(double)Level1_Cnt_mwh.f[f].count_tot),
					/*Level1_Cnt_mwh.f[f].count_tot, */
					Level1_Cnt_mwh.f[f].num_elems,
					thresholdPrct);
	
	
			}
		}
		
		// Generate our new wheels
		generate_wheels_lvl1_mwh(f);	
	}
}

void adjust_local1_weights_with_counters_mwh(double prctChnge,double thresholdPrct) 
/*
Update all probabilities based on the counts for lcl1 functions
and then update all wheels. Right now always calculates all wheels immediately
*/
{
	int f,c;
	for (f = 0; f < Local1_mwh.total; f++) {
		Local1_mwh.f[f].wheel_calculated = 0;
		Local1_t_mwh.f[f].wheel_calculated = 0;
		Local1_f_mwh.f[f].wheel_calculated = 0;
	        if (Local1_mwh.f[f].count_tot) {
			for (c = 0; c < Local1_mwh.f[f].num_elems;c++) {
				Local1_mwh.f[f].prob[c] =
					 
			   	acgp_new_weight_lvl1_mwh(
			        	Local1_mwh.f[f].prob[c],
			   		prctChnge,
					(double)((double)Local1_mwh.f[f].count[c]/(double)Local1_mwh.f[f].count_tot),
					/*Level1_Cnt_mwh.f[f].count_tot, */
					Local1_mwh.f[f].num_elems,
					thresholdPrct);
	
	
			}
		}
		
		// Generate our new wheels
		generate_regrow_wheel_lcl1_mwh(f);
	}

}

void adjust_weights_with_MS_czj_mwh() {
	
	// Pull our probabilities from MS_czj
	pull_L2_probs_from_MS_czj();

	// now recompute wheels
	int f;
	for (f = 0; f < Level1_Cnt_mwh.total; f++) {
		// Generate our new wheels
		generate_wheels_lvl1_mwh(f);	
	}
	
	for (f=0; f < Local1_mwh.total; f++) {
		// Generate our new wheels
		generate_regrow_wheel_lcl1_mwh(f);
	}

}

void adjust_weights_with_counters_mwh(double prctChnge, double thresholdPrct)
/*
Updates all probabilities based on counts for all functions
and then updates all wheels
*/
{
	// pull from MS_czj for validation only
	// adjust_weights_with_MS_czj_mwh();
	
	//Change this back after testing
	adjust_lvl1_weights_with_counters_mwh(prctChnge,thresholdPrct);
	adjust_local1_weights_with_counters_mwh(prctChnge,thresholdPrct);

}

void print_header_lvl1_mwh(FILE *fp)
/*
   Prints header information into a .gcnt or .gwgt file
*/
{
	fprintf(fp,"Allowable Roots = %d\n",Level1_Cnt_mwh.total);
	int f;
	int NumFT = fset[0].function_count;
	// Print the list of allowable roots
	for (f = 0; f < Level1_Cnt_mwh.total; f++) {
		int fsetrootindex = MS_czj[NumFT][0].members[f];
		fprintf(fp,"%s,",fset[0].cset[fsetrootindex].string);
	}

	fprintf(fp,"\n\n");

	for (f=0; f < Level1_Cnt_mwh.total; f++) {
		int fsetrootindex = MS_czj[NumFT][0].members[f];
		fprintf(fp,"'%s' has arity = %d\n",fset[0].cset[fsetrootindex].string,Level1_Cnt_mwh.f[f].dim);
		
		// Loop through the arity's of our 'f' function
		int i;
		for (i = 0; i < Level1_Cnt_mwh.f[f].dim; i++) {
			fprintf(fp,"Child_%d,total = %d\n",i,Level1_Cnt_mwh.f[f].size[i]);
			int c;
			// Loop through all possible children for this arity of the root
			for (c = 0; c < Level1_Cnt_mwh.f[f].size[i]; c++) {
				int fset_index = MS_czj[fsetrootindex][i].members[c];
				fprintf(fp,"%s,",fset[0].cset[fset_index].string);
			}
			fprintf(fp,"\n");
		}
		fprintf(fp,"\n");
	}
}

void print_header_local1_mwh(FILE *fp)
/*
   Prints header information into a .lcnt or .lwgt file
*/
{
	fprintf(fp,"Allowable Parents = %d\n",Level1_Cnt_mwh.total);
	int f;
	int NumFT = fset[0].function_count;
	// Print the list of allowable parents
	for (f = 0; f < Local1_mwh.total; f++) {
		fprintf(fp,"%s,",fset[0].cset[f].string);
	}

	fprintf(fp,"\n\n");

	for (f=0; f < Local1_mwh.total; f++) {
		fprintf(fp,"'%s' has arity = %d\n",fset[0].cset[f].string,Local1_mwh.f[f].dim);
		
		// Loop through the arity's of our 'f' function
		int i;
		for (i = 0; i < Local1_mwh.f[f].dim; i++) {
			fprintf(fp,"Child_%d,total = %d\n",i,Local1_mwh.f[f].size[i]);
			int c;
			// Loop through all possible children for this arity of the root
			for (c = 0; c < Local1_mwh.f[f].size[i]; c++) {
				int fset_index = MS_czj[f][i].members[c];
				fprintf(fp,"%s,",fset[0].cset[fset_index].string);
			}
			fprintf(fp,"\n");
		}
		fprintf(fp,"\n");
	}

}

void print_lvl1_cnts_mwh(int curGen,const char *basefile,int flag)
	/* Display for debugging all level 1 cnts
		Note: Large for anything but extremely trivial problems */
{
	FILE *fp;
	char fname[BUFSIZ];
	
	strcpy(fname,basefile);
	strcat(fname,".gcnt");
	if ((fp=fopen(fname,"a+"))==0)
	{
		fprintf(stderr,"Could not open %s to write counts\n",fname);
		exit(1);
	}
	
	int NumFT = fset[0].function_count;
	
	if (curGen == 0)
		print_header_lvl1_mwh(fp);
	
	fprintf(fp,"\nGen %d",curGen);
	int f,i;
	for (f = 0; f < Level1_Cnt_mwh.total; f++) {
		int fsetrootindex = MS_czj[NumFT][0].members[f];
		
		
		fprintf(fp,"\n\tRoot %d ( %s ): arity %d,total elements %d,trees with this root=%d\n",f,
				fset[0].cset[fsetrootindex].string,		
							Level1_Cnt_mwh.f[f].dim,
							Level1_Cnt_mwh.f[f].num_elems,
							Level1_Cnt_mwh.f[f].count_tot);
		fprintf(fp,"Offsets=[");
		for (i = 0; i < Level1_Cnt_mwh.f[f].dim; i++)
			fprintf(fp,"%d,",Level1_Cnt_mwh.f[f].offset[i]);
		fprintf(fp,"]\n");
		fprintf(fp,"Sizes=[");
		for (i = 0; i < Level1_Cnt_mwh.f[f].dim; i++)
			fprintf(fp,"%d,",Level1_Cnt_mwh.f[f].size[i]);
		fprintf(fp,"]\n");
		int c;
		int i = 0; 
		int *this_index = (int *) calloc(Level1_Cnt_mwh.f[f].dim,sizeof(int));
		if (this_index==NULL) {
			printf("Couldnt allocate this_index array in acgp_mwh\n");
			exit(1);
		}
	
		
		for (c = 0; c < Level1_Cnt_mwh.f[f].num_elems;c++) {
			// first print generation
			fprintf(fp,"%d %d ",curGen,this_iteration);
			
			// next print root
			fprintf(fp,"%s ",fset[0].cset[fsetrootindex].string);
			
			//now we want to find the indexes of this 'c'
			convert_int_to_indices(this_index,c,Level1_Cnt_mwh.f[f].dim,Level1_Cnt_mwh.f[f].offset);
					
			// now this_index[i] is the i'th child of f
			for (i = 0; i < Level1_Cnt_mwh.f[f].dim; i++) {
				fprintf(fp,"%s ",fset[0].cset[MS_czj[f][i].members[this_index[i]]].string);
			}
			
			fprintf(fp,"%d\n",Level1_Cnt_mwh.f[f].count[c]);
		
		}
						
		free(this_index);
	}
	fclose(fp);
}



void print_local1_cnts_mwh(int curGen,const char *basefile,int flag)
	/* Display for debugging all level 1 cnts
		Note: Large for anything but extremely trivial problems */
{
	FILE *fp;
	char fname[BUFSIZ];
	
	strcpy(fname,basefile);
	strcat(fname,".lcnt");
	if ((fp=fopen(fname,"a+"))==0)
	{
		fprintf(stderr,"Could not open %s to write counts\n",fname);
		exit(1);
	}
	
	int NumFT = fset[0].function_count;
	
	if (curGen == 0)
		print_header_local1_mwh(fp);
	
	fprintf(fp,"\n\nGen %d",curGen);
	int f,i;
	for (f = 0; f < Local1_mwh.total; f++) {
		
		// not needed for local counters
		//int fsetrootindex = MS_czj[NumFT][0].members[f];
		
		
		fprintf(fp,"\n\tParent %d ( %s ): arity %d,total elements %d,trees with this as parent=%d\n",f,
				fset[0].cset[/*fsetrootindex*/ f].string,		
							Local1_mwh.f[f].dim,
							Local1_mwh.f[f].num_elems,
							Local1_mwh.f[f].count_tot);
		fprintf(fp,"Offsets=[");
		for (i = 0; i < Local1_mwh.f[f].dim; i++)
			fprintf(fp,"%d,",Local1_mwh.f[f].offset[i]);
		fprintf(fp,"]\n");
		fprintf(fp,"Sizes=[");
		for (i = 0; i < Local1_mwh.f[f].dim; i++)
			fprintf(fp,"%d,",Local1_mwh.f[f].size[i]);
		fprintf(fp,"]\n");
		int c;
		int i = 0; 
		int *this_index = (int *) calloc(Local1_mwh.f[f].dim,sizeof(int));
		if (this_index==NULL) {
			printf("Couldnt allocate this_index array in acgp_mwh\n");
			exit(1);
		}
		
		
		for (c = 0; c < Local1_mwh.f[f].num_elems;c++) {
			// first print generation
			fprintf(fp,"%d %d ",curGen,this_iteration);
			
			// next print root
			fprintf(fp,"%s ",fset[0].cset[f].string);
			
			//now we want to find the indexes of this 'c'
			convert_int_to_indices(this_index,c,Local1_mwh.f[f].dim,Local1_mwh.f[f].offset);
					
			// now this_index[i] is the i'th child of f
			for (i = 0; i < Local1_mwh.f[f].dim; i++) {
				fprintf(fp,"%s ",fset[0].cset[MS_czj[f][i].members[this_index[i]]].string);
			}
			
			fprintf(fp,"%d\n",Local1_mwh.f[f].count[c]);
		
		}
		
				
		free(this_index);
	}
	fclose(fp);
}

void print_local1_weights_mwh(int curGen,const char *basefile,int flag)
	/* Display for debugging all level 1 wgts
		Note: Large for anything but extremely trivial problems */
{
	FILE *fp;
	char fname[BUFSIZ];
	
	strcpy(fname,basefile);
	strcat(fname,".lwgt");
	if ((fp=fopen(fname,"a+"))==0)
	{
		fprintf(stderr,"Could not open %s to write counts\n",fname);
		exit(1);
	}
	
	int NumFT = fset[0].function_count;
	
	if (curGen == 0)
		print_header_local1_mwh(fp);
	
	fprintf(fp,"\n\nGen %d",curGen);
	int f,i;
	for (f = 0; f < Local1_mwh.total; f++) {
		
		// not needed for local counters
		//int fsetrootindex = MS_czj[NumFT][0].members[f];
		
		
		fprintf(fp,"\n\tParent %d ( %s ): arity %d,total elements %d,trees with this as parent=%d\n",f,
				fset[0].cset[/*fsetrootindex*/ f].string,		
							Local1_mwh.f[f].dim,
							Local1_mwh.f[f].num_elems,
							Local1_mwh.f[f].count_tot);
		fprintf(fp,"Offsets=[");
		for (i = 0; i < Local1_mwh.f[f].dim; i++)
			fprintf(fp,"%d,",Local1_mwh.f[f].offset[i]);
		fprintf(fp,"]\n");
		fprintf(fp,"Sizes=[");
		for (i = 0; i < Local1_mwh.f[f].dim; i++)
			fprintf(fp,"%d,",Local1_mwh.f[f].size[i]);
		fprintf(fp,"]\n");
		int c;
		int i = 0; 
		int *this_index = (int *) calloc(Local1_mwh.f[f].dim,sizeof(int));
		if (this_index==NULL) {
			printf("Couldnt allocate this_index array in acgp_mwh\n");
			exit(1);
		}
		
		for (c = 0; c < Local1_mwh.f[f].num_elems;c++) {
			// first print generation
			fprintf(fp,"%d %d ",curGen,this_iteration);
			
			// next print root
			fprintf(fp,"%s ",fset[0].cset[f].string);
			
			//now we want to find the indexes of this 'c'
			convert_int_to_indices(this_index,c,Local1_mwh.f[f].dim,Local1_mwh.f[f].offset);
					
			// now this_index[i] is the i'th child of f
			for (i = 0; i < Local1_mwh.f[f].dim; i++) {
				fprintf(fp,"%s ",fset[0].cset[MS_czj[f][i].members[this_index[i]]].string);
			}
			
			fprintf(fp,"%f\n",Local1_mwh.f[f].prob[c]);
		
		}		
		free(this_index);
	}
	fclose(fp);
}

void print_lvl1_weights_mwh(int curGen,const char *basefile,int flag)
	/* Display for debugging all level 1 cnts
		Note: Large for anything but extremely trivial problems */
{
	FILE *fp;
	char fname[BUFSIZ];
	
	strcpy(fname,basefile);
	strcat(fname,".gwgt");
	if ((fp=fopen(fname,"a+"))==0)
	{
		fprintf(stderr,"Could not open %s to write counts\n",fname);
		exit(1);
	}
	
	if (curGen == 0)
		print_header_lvl1_mwh(fp);
	
	int NumFT = fset[0].function_count;
	fprintf(fp,"\n\nGen %d",curGen);
	int f,i;
	for (f = 0; f < Level1_Cnt_mwh.total; f++) {
		int fsetrootindex = MS_czj[NumFT][0].members[f];
		fprintf(fp,"\n\tRoot %d (%s): arity %d,total elements %d,trees with this root=%d\n",f,
							fset[0].cset[fsetrootindex].string,
							Level1_Cnt_mwh.f[f].dim,
							Level1_Cnt_mwh.f[f].num_elems,
							Level1_Cnt_mwh.f[f].count_tot);
		fprintf(fp,"Offsets=[");
		for (i = 0; i < Level1_Cnt_mwh.f[f].dim; i++)
			fprintf(fp,"%d ",Level1_Cnt_mwh.f[f].offset[i]);
		fprintf(fp,"]\n");
		fprintf(fp,"Sizes=[");
		for (i = 0; i < Level1_Cnt_mwh.f[f].dim; i++)
			fprintf(fp,"%d ",Level1_Cnt_mwh.f[f].size[i]);
		fprintf(fp,"]\n");
		int c;
		int i = 0;
		
		int *this_index = (int *) calloc(Level1_Cnt_mwh.f[f].dim,sizeof(int));
		if (this_index==NULL) {
			printf("Couldnt allocate this_index array in acgp_mwh\n");
			exit(1);
		}
		
		
		for (c = 0; c < Level1_Cnt_mwh.f[f].num_elems;c++) {
			// first print generation
			fprintf(fp,"%d %d ",curGen,this_iteration);
			
			// next print root
			fprintf(fp,"%s ",fset[0].cset[fsetrootindex].string);
			
			//now we want to find the indexes of this 'c'
			convert_int_to_indices(this_index,c,Level1_Cnt_mwh.f[f].dim,Level1_Cnt_mwh.f[f].offset);
					
			// now this_index[i] is the i'th child of f
			for (i = 0; i < Level1_Cnt_mwh.f[f].dim; i++) {
				fprintf(fp,"%s ",fset[0].cset[MS_czj[f][i].members[this_index[i]]].string);
			}
			
			fprintf(fp,"%f\n",Level1_Cnt_mwh.f[f].prob[c]);
		
		}
		free(this_index);
	}
	fclose(fp);
}


void print_cnts_mwh(int curGen,const char *basefile,int flag)
/* Prints counts to count file
*/
{
	if (flag)
		this_iteration++;
		
	print_lvl1_cnts_mwh(curGen,basefile,flag);
	print_local1_cnts_mwh(curGen,basefile,flag);
	
}

void print_weights_mwh(int curGen,const char *basefile,int flag)
/* prints weights to weight file
*/
{
	print_lvl1_weights_mwh(curGen,basefile,flag);
	print_local1_weights_mwh(curGen,basefile,flag);

}

void list_counted_level1_mwh(int curGen, const char *baseline, int flag) 
	/* Prints out a list of all counted level 1 trees for debugging purposes 
							*/
{
	int NumFT = fset[0].function_count;
	FILE *fp;
	char fname[BUFSIZ];
	strcpy(fname,baseline);
	strcat(fname,".gcnt");
	if ((fp=fopen(fname,"a+"))==0)
	{
		fprintf(stderr,"Could not open %s to write counts\n",fname);
		exit(1);
	}

        fprintf(fp,"======================================\n");
	fprintf(fp,"MWH_Displaying level 1 trees stored at generation %d\n",curGen);
	int f,i;
	
	// Loop through all our functions at root
	for (f=0; f < Level1_Cnt_mwh.total; f++) {
		int c,remain;
		i = 0;
		int *this_index = (int *) calloc(Level1_Cnt_mwh.f[f].dim,sizeof(int));
		if (this_index==NULL) {
			printf("Couldnt allocate this_index array in acgp_mwh\n");
			exit(1);
		}
		// Go through counter and pick out trees, then display them
		for (c = 0; c < Level1_Cnt_mwh.f[f].num_elems; c++) {
			if (Level1_Cnt_mwh.f[f].count[c] == 0)
				continue;
				
			remain = c;
			// extract from our counter the tree into this_index array
			for (i = Level1_Cnt_mwh.f[f].dim-1; i >=0; i--) {
				this_index[i] = remain / Level1_Cnt_mwh.f[f].offset[i];
			//	printf("Remain = %d, doing mod %d with it, get %d",remain,Level1_Cnt_mwh.f[f].offset[i],this_index[i]);
				
				remain = remain - this_index[i]*Level1_Cnt_mwh.f[f].offset[i];
				
			}
#ifdef DEBUG_LVL1
			printf(", Raw indices are [");
			for (i = 0; i < Level1_Cnt_mwh.f[f].dim;i++)
				printf("%d ",this_index[i]);
			// Use MS_czj to convert these raw to actual fset ones
			printf("], Converted indices are [");
			for (i = 0; i < Level1_Cnt_mwh.f[f].dim; i++)
				printf("%d ",MS_czj[f][i].members[this_index[i]]);
			printf("]\n");
			
#endif
			int index_to_fset = MS_czj[NumFT][0].members[f];   // was just 'f'
			fprintf(fp,"%d of : %s(",Level1_Cnt_mwh.f[f].count[c],fset[0].cset[index_to_fset].string);
			for (i = 0; i < Level1_Cnt_mwh.f[f].dim; i++)
				fprintf(fp,"%s ",fset[0].cset[MS_czj[f][i].members[this_index[i]]].string);
			fprintf(fp,")\n");
		
		}
		free(this_index);
	}
	
	fprintf(fp,"======================================\n");	
	fclose(fp);
}

void update_context_mwh(int depth,int index)
/*
   Gets as input the depth(starting at 1 for the root) and current function and then determines if
   we want to add that to one of our context information arrays 
*/
{
        if (depth == 2) {
		level1context_mwh[size_context1_mwh] = index;
		size_context1_mwh++;
	}
	
}

void reset_counters_level1_mwh(void) 
	/* Reset all counters to 0 */
{	
	int i;
	// Loop through all our hypercubes
	for (i = 0; i < Level1_Cnt_mwh.total; i++) {
		
		Level1_Cnt_mwh.f[i].count_tot = 0;
		// Loop through this functions hypercube
		int j;
		for (j= 0; j < Level1_Cnt_mwh.f[i].num_elems; j++) {
			Level1_Cnt_mwh.f[i].count[j] = 0;
		}
	}
}

void reset_counters_local1_mwh(void)
{	
	int i;
	// Loop through all our hypercubes
	for (i = 0; i < Local1_mwh.total; i++) {
		
		Local1_mwh.f[i].count_tot = 0;
		// Loop through this functions hypercube
		int j;
		for (j= 0; j < Local1_mwh.f[i].num_elems; j++) {
			Local1_mwh.f[i].count[j] = 0;
		}
	}
}


void reset_counters_mwh(void) 
	/* Reset all counters to 0 */
{	
	reset_counters_level1_mwh();
	reset_counters_local1_mwh();
}

void pull_L2_probs_from_MS_czj(void) 
/*
Using the probabilities from MS_czj, set the l2 probabilities
Should be used to set initial probabilties, as well as possibly in testing
*/
{

	int parIndex;
	int NumF = fset[0].function_count;
	
	// First do it for Global probabilities
	for (parIndex = 0; parIndex < Level1_Cnt_mwh.total ; parIndex++) {
		Level1_Cnt_mwh.f[parIndex].wheel_calculated = 0;
		Level1_Cnt_t_mwh.f[parIndex].wheel_calculated = 0;
		Level1_Cnt_f_mwh.f[parIndex].wheel_calculated = 0;
		
		int *this_index = (int *) calloc(Level1_Cnt_mwh.f[parIndex].dim,sizeof(int)); /* index to this element */
		
		// indicate wheel is not updated
		Level1_Cnt_mwh.f[parIndex].wheel_calculated = 0;
	
		// convert parIndex into index to MS_czj
		int trueIndex = MS_czj[NumF][0].members[parIndex];
		int c,remain;
		for (c = 0; c < Level1_Cnt_mwh.f[parIndex].num_elems; c++) {
			
			// convert c into an index into this cube
			
			convert_int_to_indices(this_index,c,Level1_Cnt_mwh.f[parIndex].dim,Level1_Cnt_mwh.f[parIndex].offset);
			
			// Now we want to multiply together the probabilities p(index[i]|parent)
			double this_prob = 1.0;
			double this_weight = 0;
			
			int isminweight = 0;
						
			// Now lets get the children probabilities given parent
			int child;
			for (child = 0; child < Level1_Cnt_mwh.f[parIndex].dim; child++) {
											
				// get the weight of this child given the parent
				int thisindexchild;
				thisindexchild = this_index[child];
				
				/*printf("itm=%d ",thisindexchild);
				fflush(stdout); */
				
				int afterindex;
				afterindex = MS_czj[trueIndex][child].members[thisindexchild] ;
				
				/* printf("itw=%d ",afterindex);
				fflush(stdout); */
								
				this_weight =MS_czj[trueIndex][child].weights[
					afterindex  ];
				
				// printf("weight=%f",this_weight);
				
				// is this a minimal weight, such that this element should never happen?
				if (this_weight <= MINWGHT + SMALL) {
					isminweight = 1;
					break;
				}
				
				this_prob = this_prob * this_weight; /* / weight_total; */
									
			}  /* loop over children, this_index[child] */
			
			if (!isminweight) {
				Level1_Cnt_mwh.f[parIndex].prob[c] = this_prob;
			}
			else {
				Level1_Cnt_mwh.f[parIndex].prob[c] = L2minweight;
			}
			
		}  /* loop over all elements in counter cube */
		free(this_index);
	}  /* loop over possible roots */	
	
	
	// Now do it for Local probabilities
	for (parIndex = 0; parIndex < Local1_mwh.total ; parIndex++) {
		Local1_mwh.f[parIndex].wheel_calculated = 0;
		Local1_t_mwh.f[parIndex].wheel_calculated = 0;
		Local1_f_mwh.f[parIndex].wheel_calculated = 0;
		int *this_index = (int *) calloc(Local1_mwh.f[parIndex].dim,sizeof(int)); /* index to this element */
		
		// indicate wheel is not updated
		Local1_mwh.f[parIndex].wheel_calculated = 0;
	
		// convert parIndex into index to MS_czj
		int trueIndex = parIndex;   // not needed I dont think atm for local counters
		int c,remain;
		for (c = 0; c < Local1_mwh.f[parIndex].num_elems; c++) {
			
			// convert c into an index into this cube
			
			convert_int_to_indices(this_index,c,Local1_mwh.f[parIndex].dim,Local1_mwh.f[parIndex].offset);
			
			// Now we want to multiply together the probabilities p(index[i]|parent)
			double this_prob = 1.0;
			double this_weight = 0;
			
			/*
			printf("c=%d,parentIndex=%d,indexes are [",c,parIndex);
			for (i = 0; i < Level1_Cnt_mwh.f[parIndex].dim;i++) {
				printf("%d ",this_index[i]);
			}
			printf("\n"); fflush(stdout); */
			
			int isminweight = 0;
						
			// Now lets get the children probabilities given parent
			int child;
			for (child = 0; child < Local1_mwh.f[parIndex].dim; child++) {
											
				// get the weight of this child given the parent
				int thisindexchild;
				thisindexchild = this_index[child];
				
				/*printf("itm=%d ",thisindexchild);
				fflush(stdout); */
				
				int afterindex;
				afterindex = MS_czj[trueIndex][child].members[thisindexchild] ;
				
				/* printf("itw=%d ",afterindex);
				fflush(stdout); */
								
				this_weight =MS_czj[trueIndex][child].weights[
					afterindex  ];
				
				//printf("weight=%f",this_weight);
				
				// is this a minimal weight, such that this element should never happen?
				if (this_weight <= MINWGHT + SMALL) {
					isminweight = 1;
					break;
				}
				
				this_prob = this_prob * this_weight; /* / weight_total; */
									
			}  /* loop over children, this_index[child] */
			
			Local1_mwh.f[parIndex].prob[c] = this_prob;
			
			if (!isminweight) {
				Local1_mwh.f[parIndex].prob[c] = this_prob;
			}
			else {
				Local1_mwh.f[parIndex].prob[c] = L2minweight;
			}
			
		}  /* loop over all elements in counter cube */
		free(this_index);
	}  /* loop over possible roots */
}


void create_mwh_level1(void) 
	/* Allocate and initialize complete level 1 counters(Level1_Cnt)
	 from the root */
	
{
	// Need to remove our old counts file
	FILE *fp;
	const char *basefile;
	char fname[BUFSIZ];
	basefile = get_parameter("output.basename");
	
	char *p;
	if (p=get_parameter("acgp.level")) {
		Acgp_level = atof(p);
	}
	else {
		Acgp_level = 1;
	}
	
	if (Acgp_level < 2) {
		printf("Initializing level 2 heuristics but acgp_level not greater than 1!\n");
		printf("Please set acgp_level to greater than 1 and try again.\n");
		exit(1);
	}
	
	if (Acgp_level < 3) {
		strcpy(fname,basefile);
		strcat(fname,".gcnt");
	
		if ((fp=fopen(fname,"w"))==0) {
			printf("Could not open lvl1 counter file in create_mwh\n");
			exit(1);
		}
		fclose(fp);
	
		strcpy(fname,basefile);
		strcat(fname,".gwgt");
		if ((fp=fopen(fname,"w"))==0) {
			printf("Could not open lvl1 weights file in create_mwh\n");
		}
		fclose(fp);
		
		strcpy(fname,basefile);
		strcat(fname,".lcnt");
		if ((fp=fopen(fname,"w"))==0) {
			printf("Could not open local lvl1 counter file in create_mwh\n");
		}
		fclose(fp);
	
		strcpy(fname,basefile);
		strcat(fname,".lwgt");
		if ((fp=fopen(fname,"w"))==0) {
			printf("Could not open local lvl1 weights file in create_mwh\n");
		}
		fclose(fp);
	}

	int NumF = fset[0].function_count;
	//printf("Creating Lvl1 Counters\n");
	Level1_Cnt_t_mwh.total = Level1_Cnt_f_mwh.total = Level1_Cnt_mwh.total = MS_czj[NumF][0].numF;
	
	//printf("Total hypercubes at level 1 will be %d\n",Level1_Cnt_mwh.total);
	
	Level1_Cnt_mwh.f = (typeHyperCube *)calloc(Level1_Cnt_mwh.total,sizeof(typeHyperCube));
	Level1_Cnt_f_mwh.f = (typeHyperCube *)calloc(Level1_Cnt_f_mwh.total,sizeof(typeHyperCube));
	Level1_Cnt_t_mwh.f = (typeHyperCube *)calloc(Level1_Cnt_t_mwh.total,sizeof(typeHyperCube));
	if (Level1_Cnt_mwh.f == NULL) {
		printf("Could not allocate our initial array of hypercube counters in acgp_mwh.c\n");
		exit(1);
	}
	
	int f;   // index into our array of hypercubes
	for (f = 0; f < Level1_Cnt_mwh.total; f++) {
		int parIndex = MS_czj[NumF][0].members[f];     /* Find f's index into fset and MS_czj*/
		Level1_Cnt_t_mwh.f[f].dim = Level1_Cnt_f_mwh.f[f].dim = Level1_Cnt_mwh.f[f].dim = fset[0].cset[parIndex].arity;  /* Find f's arity */
		
		
		// Allocate our size array
		Level1_Cnt_mwh.f[f].size = (int *) calloc(Level1_Cnt_mwh.f[f].dim,sizeof(int));
		Level1_Cnt_f_mwh.f[f].size = (int *) calloc(Level1_Cnt_f_mwh.f[f].dim,sizeof(int));
		Level1_Cnt_t_mwh.f[f].size = (int *) calloc(Level1_Cnt_t_mwh.f[f].dim,sizeof(int));
		if (Level1_Cnt_mwh.f[f].size == NULL) {
			printf("Could not allocate size array for hypercube %d in acgp_mwh.c\n",f);
			exit(1);
		}
		
		Level1_Cnt_f_mwh.f[f].wheel_calculated = Level1_Cnt_t_mwh.f[f].wheel_calculated = Level1_Cnt_mwh.f[f].wheel_calculated = 0;
		Level1_Cnt_f_mwh.f[f].wheel = Level1_Cnt_t_mwh.f[f].wheel = Level1_Cnt_mwh.f[f].wheel = NULL;
		
		// Now find the sizes of each of our dimensions
		int i;
		Level1_Cnt_mwh.f[f].num_elems = 1;
		Level1_Cnt_f_mwh.f[f].num_elems = 1;
		Level1_Cnt_t_mwh.f[f].num_elems = 1;
		for (i = 0; i < Level1_Cnt_mwh.f[f].dim; i++) {
			Level1_Cnt_mwh.f[f].size[i] = MS_czj[parIndex][i].numFT;   // size of this arity
			Level1_Cnt_f_mwh.f[f].size[i] = MS_czj[parIndex][i].numF;
			Level1_Cnt_t_mwh.f[f].size[i] = MS_czj[parIndex][i].numT;
			Level1_Cnt_mwh.f[f].num_elems *= Level1_Cnt_mwh.f[f].size[i];
			Level1_Cnt_f_mwh.f[f].num_elems *= Level1_Cnt_f_mwh.f[f].size[i];
			Level1_Cnt_t_mwh.f[f].num_elems *= Level1_Cnt_t_mwh.f[f].size[i];
		}
		
		// Allocate adn intialize one hypercubes counters
		// Do not need for terminal and function wheels
		Level1_Cnt_mwh.f[f].count 
			= (int *)calloc(Level1_Cnt_mwh.f[f].num_elems,sizeof(int));
		Level1_Cnt_mwh.f[f].prob
			= (double *)calloc(Level1_Cnt_mwh.f[f].num_elems,sizeof(double));
		if (Level1_Cnt_mwh.f[f].count == NULL || Level1_Cnt_mwh.f[f].prob == NULL) {
			printf("Could not allocate hypercube %d in acgp_mwh.c\n",f);
			exit(1);
		} 
		
				
		// Allocate and initialize offset array
		Level1_Cnt_mwh.f[f].offset
			= (int *) calloc(Level1_Cnt_mwh.f[f].dim,sizeof(int));
		Level1_Cnt_f_mwh.f[f].offset
			= (int *) calloc(Level1_Cnt_f_mwh.f[f].dim,sizeof(int));
		Level1_Cnt_t_mwh.f[f].offset
			= (int *) calloc(Level1_Cnt_t_mwh.f[f].dim,sizeof(int));
		if (Level1_Cnt_mwh.f[f].offset == NULL) {
			printf("Could not allocate offset array for %d in acgp_mwh.c\n",f);
			exit(1);
		}
		int this_index;
		Level1_Cnt_mwh.f[f].offset[0] = 1;
		Level1_Cnt_f_mwh.f[f].offset[0] = 1;
		Level1_Cnt_t_mwh.f[f].offset[0] = 1;
		for (this_index = 1; this_index < Level1_Cnt_mwh.f[f].dim; this_index++) {
			Level1_Cnt_mwh.f[f].offset[this_index] = 
			    Level1_Cnt_mwh.f[f].offset[this_index-1] * 
			    	Level1_Cnt_mwh.f[f].size[this_index-1];
			Level1_Cnt_f_mwh.f[f].offset[this_index] = 
			    Level1_Cnt_f_mwh.f[f].offset[this_index-1] * 
			    	Level1_Cnt_f_mwh.f[f].size[this_index-1];
			Level1_Cnt_t_mwh.f[f].offset[this_index] = 
			    Level1_Cnt_t_mwh.f[f].offset[this_index-1] * 
			    	Level1_Cnt_t_mwh.f[f].size[this_index-1];
			    /* MS_czj[parIndex][this_index].numFT; */
		}
		/*
		printf("Offsets =[");
		for (i = 0; i < Level1_Cnt_mwh.f[f].dim; i++)
			printf("%d ",Level1_Cnt_mwh.f[f].offset[i]);
		printf("]\n");
		*/
	
	}
	
			
	//printf("Allocating our context information\n");
	// TODO: Fix to max arity of any function
	level1context_mwh = (int *) calloc(10,sizeof(int));

	// =================================================================================
	// Now allocate our global local counts Local1_mwh
	
	
	//int NumF = fset[0].function_count;
	//printf("Creating Local lcl1 Counters\n");
	Local1_f_mwh.total = Local1_t_mwh.total = Local1_mwh.total = fset[0].function_count;
	
	
	//printf("Total hypercubes for local level 1 will be %d\n",Local1_mwh.total);
	
	Local1_mwh.f = (typeHyperCube *)calloc(Local1_mwh.total,sizeof(typeHyperCube));
	Local1_t_mwh.f = (typeHyperCube *)calloc(Local1_t_mwh.total,sizeof(typeHyperCube));
	Local1_f_mwh.f = (typeHyperCube *)calloc(Local1_f_mwh.total,sizeof(typeHyperCube));
	
	if (Local1_mwh.f == NULL) {
		printf("Could not allocate our initial array of local hypercube counters in acgp_mwh.c\n");
		exit(1);
	}
	
	for (f = 0; f < Local1_mwh.total; f++) {
		int parIndex = f;
		Local1_f_mwh.f[f].dim = Local1_t_mwh.f[f].dim = Local1_mwh.f[f].dim = fset[0].cset[parIndex].arity;  /* Find f's arity */
		
		//printf("Dim of %d=%d\n",f,Local1_mwh.f[f].dim);
	
		// Allocate our size array
		Local1_mwh.f[f].size = (int *) calloc(Local1_mwh.f[f].dim,sizeof(int));
		Local1_t_mwh.f[f].size = (int *) calloc(Local1_t_mwh.f[f].dim,sizeof(int));
		Local1_f_mwh.f[f].size = (int *) calloc(Local1_f_mwh.f[f].dim,sizeof(int));
		if (Local1_mwh.f[f].size == NULL) {
			printf("Could not allocate size array for local hypercube %d in acgp_mwh.c\n",f);
			exit(1);
		}
		
		Local1_t_mwh.f[f].wheel_calculated = Local1_mwh.f[f].wheel_calculated = 0;
		Local1_f_mwh.f[f].wheel = Local1_t_mwh.f[f].wheel = Local1_mwh.f[f].wheel = NULL;
		
		// Now find the sizes of each of our dimensions
		int i;
		Local1_mwh.f[f].num_elems = 1;
		Local1_t_mwh.f[f].num_elems = 1;
		Local1_f_mwh.f[f].num_elems = 1;
		for (i = 0; i < Local1_mwh.f[f].dim; i++) {
			Local1_mwh.f[f].size[i] = MS_czj[parIndex][i].numFT;   // size of this arity
			Local1_t_mwh.f[f].size[i] = MS_czj[parIndex][i].numT;
			Local1_f_mwh.f[f].size[i] = MS_czj[parIndex][i].numF;
			
			Local1_mwh.f[f].num_elems *= Local1_mwh.f[f].size[i];
			Local1_t_mwh.f[f].num_elems *= Local1_t_mwh.f[f].size[i];
			Local1_f_mwh.f[f].num_elems *= Local1_f_mwh.f[f].size[i];
		}
		
		// Allocate adn intialize one hypercubes counters
		// Note: Terminal and function wheels do not need this
		Local1_mwh.f[f].count 
			= (int *)calloc(Local1_mwh.f[f].num_elems,sizeof(int));
		Local1_mwh.f[f].prob
			= (double *)calloc(Local1_mwh.f[f].num_elems,sizeof(double));
		if (Local1_mwh.f[f].count == NULL || Local1_mwh.f[f].prob == NULL) {
			printf("Could not allocate local hypercube %d in acgp_mwh.c\n",f);
			exit(1);
		} 
		
				
		// Allocate and initialize offset array
		Local1_mwh.f[f].offset
			= (int *) calloc(Local1_mwh.f[f].dim,sizeof(int));
		Local1_t_mwh.f[f].offset
			= (int *) calloc(Local1_t_mwh.f[f].dim,sizeof(int));
		Local1_f_mwh.f[f].offset
			= (int *) calloc(Local1_f_mwh.f[f].dim,sizeof(int));
		if (Local1_mwh.f[f].offset == NULL) {
			printf("Could not allocate local offset array for %d in acgp_mwh.c\n",f);
			exit(1);
		}
		
		int this_index;
		Local1_mwh.f[f].offset[0] = 1;
		Local1_t_mwh.f[f].offset[0] = 1;
		Local1_f_mwh.f[f].offset[0] = 1;
		for (this_index = 1; this_index < Local1_mwh.f[f].dim; this_index++) {
			Local1_mwh.f[f].offset[this_index] = 
			    Local1_mwh.f[f].offset[this_index-1] * 
			    	Local1_mwh.f[f].size[this_index-1];
			Local1_t_mwh.f[f].offset[this_index] =
				Local1_t_mwh.f[f].offset[this_index-1] *
					Local1_t_mwh.f[f].size[this_index-1];
			Local1_f_mwh.f[f].offset[this_index] =
				Local1_f_mwh.f[f].offset[this_index-1] *
					Local1_f_mwh.f[f].size[this_index-1];
			    /* MS_czj[parIndex][this_index].numFT; */
		}
		/*
		printf("Offsets =[");
		for (i = 0; i < Local1_mwh.f[f].dim; i++)
			printf("%d ",Local1_mwh.f[f].offset[i]);
		printf("]\n");
	
		printf("Sizes =[");
		for (i = 0; i < Local1_mwh.f[f].dim; i++)
			printf("%d ",Local1_mwh.f[f].size[i]);
		printf("]\n");
		*/
	}
	
	// now initialize our starting probabilities for our global and local counters
	pull_L2_probs_from_MS_czj();

	// Create cache to convert from fset to mutation set
	
	int numFT = fset[0].function_count + fset[0].terminal_count;
	int numT = fset[0].terminal_count;
	// First take care of all but the root
	int i;
	for (f = 0; f < NumF; f++) {
		int a;
		for (a = 0; a < fset[0].cset[f].arity; a++) {
			MS_czj[f][a].convert_to_ms = (int *)MALLOC(numFT * sizeof(int));
		
			// first initialize all to -1
			for (i = 0; i < numFT; i++) {
				MS_czj[f][a].convert_to_ms[i] = -1;
			}
			
			for (i = 0; i < MS_czj[f][a].numFT; i++) {
				MS_czj[f][a].convert_to_ms[ MS_czj[f][a].members[i] ] = i;
			}
		}
	}
	// Now the root
	MS_czj[NumF][0].convert_to_ms = (int *)MALLOC(numFT * sizeof(int));
	for (i = 0; i < numFT; i++) {
		MS_czj[NumF][0].convert_to_ms[i] = -1;
	}
	
	for (i = 0; i < MS_czj[NumF][0].numFT; i++) {
		MS_czj[NumF][0].convert_to_ms[ MS_czj[NumF][0].members[i] ] = i;
	}
}



void create_mwh(void) {
	
	// create level one counters
	//printf("Creating level1 counters\n");
	create_mwh_level1();
	
	/*// create level two counters
	if (acgp_level == 2) {
		printf("Current level 2 counts are not implemented\n");
		exit(1);
		//create_mwh_level2();
	}
	*/
	size_context1_mwh = 0;
	size_context2_mwh = 0;
}


int gather_info_recurse_mwh ( lnode **l, lnode *sub,int *trueIndex,int *replace_depth ,int *replace_arity, int thisArity,int *context)
{
     function *parent;
     
     parent = (**l).f;
     int i, j;

     
     int found_flag = 0;

     if ( *l == sub ) {
          //return depth;
	//  printf("Found our target node in gather_info\n");
	  *replace_depth = depth_mwh;
	  *replace_arity = thisArity;
	  found_flag = 1;
     }

     int child_is_target = 0;;	  

     // mwh: Check if we are at the right depth and if so
     // add this ones index to our context information
     depth_mwh++;
   //  update_context_mwh(depth_mwh,(**l).f->index); /* update context arrays depending on depth */
    
    int *child_list;
    child_list = (int *)MALLOC(parent->arity * sizeof(int));
    
     ++*l;                      /* l is now first child */
     switch ( parent->type )
     {  case FUNC_DATA:
        case EVAL_DATA:
          for ( i = 0; i < parent->arity; ++i )
          {    
	       child_list[i] = (**l).f->index;   // save child[i] to our parent
               if (gather_info_recurse_mwh(l,sub,trueIndex,replace_depth,replace_arity,i,context) == 1)
	          child_is_target = 1; 
		depth_mwh--;   /* going back up */  
          }
          break;
        case FUNC_EXPR:
        case EVAL_EXPR:
          for ( i = 0; i < parent->arity; ++i )
          {    ++*l;           /* l was skipnode, now its the child */
              
	      child_list[i] = (**l).f->index;  // save child[i] to our parent
	      
               if (gather_info_recurse_mwh(l,sub,trueIndex,replace_depth,replace_arity,i,context) == 1)
	          child_is_target = 1;
		depth_mwh--;  /* going back up*/
          }
          break;
     }
     
     if (child_is_target == 1) {
        *trueIndex = parent->index;
//	printf("Child was found below %d, childlist is [",*trueIndex);
     	// the target of mutate or crossover is one of this child_list[]'s children
	// so lets save this context information
     	for (i=0; i < parent->arity; i++)
	{
		context[i] = child_list[i];
//		printf("%d ",child_list[i]);
	}
//	printf("]\n");
     }
    
    FREE(child_list);
    return found_flag;
}	  
    

void gather_info_mwh(lnode *data,lnode *target,int *replace_depth,int *replace_Arity,int *trueIndex,int *context)
/*
  Goes through the tree 'data' gathering context information and when it finds the correct target
  , stores its conext in 'level1context_mwh' and
  also stores its parent in trueIndex
  
*/
{
	lnode *l = data;
	function *parent;
	parent = (*l).f;
	
	depth_mwh = 0;
	size_context1_mwh = 0;
	int found = 0;
	if (gather_info_recurse_mwh(&l,target,trueIndex,replace_depth,replace_Arity,0,context) == 1)
	{
		// target was at root, in this case, *context is undefined, make sure not to use this in operators then
		//printf("Gather info found target at root\n");
		*trueIndex = parent->index;	
	}
}
