/*  lil-gp Genetic Programming System, version 1.0, 11 July 1995
 *  Copyright (C) 1995  Michigan State University
 * 
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of version 2 of the GNU General Public License as
 *  published by the Free Software Foundation.
 * 
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 * 
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 *  
 *  Douglas Zongker       (zongker@isl.cps.msu.edu)
 *  Dr. Bill Punch        (punch@isl.cps.msu.edu)
 *
 *  Computer Science Department
 *  A-714 Wells Hall
 *  Michigan State University
 *  East Lansing, Michigan  48824
 *  USA
 *  
 */
 /********* czj 
   mutation modified to always grow from the root
***********/

#include <lilgp.h>

// Our global Level1_Cnt
extern typeLevel1_Cnt Level1_Cnt_mwh;  // global count datastructure
extern typeLevel1_Cnt Local1_mwh;    // global local datastructure
extern int *level1context_mwh;	// list of the indices of lvl 1 context data
extern int *level2context_mwh;  // list of indices of lvl 2 context data
extern int size_context1_mwh;   // number of elements in level1context_mwh
extern int size_context2_mwh;   // number of elements in level2context_mwh
extern int acgp_level;

/********************************************************************/
/********************************************************************
   Mutationlvl12Local operator
**********************************************************************/
/********************************************************************/

typedef struct
{
     int keep_trying;
     double internal;
     double external;
     double *tree;
     double treetotal;
     char *sname;
     sel_context *sc;
     int method;
     int mindepth, maxdepth;
     int absDepth_czj;
} mutate_data;

/* operator_mutate_init()
 *
 * initialize a breedphase record for mutation.  much of this is cut and
 * paste from the crossover operator; see that for (a few) more comments.
 */

int operator_mutate2Local_init ( char *options, breedphase *bp )
{
     int errors = 0;
     mutate_data *md;
     int i, j, k, m;
     double r;
     char **argv, **targv;
     int internalset = 0, externalset = 0;
     char *cp;

     // mwh : set to not use level 2 heuristics unless requested
     char *p;
     if (p=get_parameter( "acgp.level")) {
         printf("We are using level 2 heuristics operator mutation2Local.\n");
         Acgp_level = atof(p);
     }
     else {
         Acgp_level = 1;
     }

     if (Acgp_level < 2) {
     	printf("Trying to use level 2 heuristic operator(mutation2Local) when acgp_level is not greater than 1!\n");
	printf("Please change to a different operator or set acgp_level=2 in parameter file\n");
	exit(1);
     }

     md = (mutate_data *)MALLOC ( sizeof ( mutate_data ) );

     /* fill in the breedphase record. */
     bp->operator = OPERATOR_MUTATION2LOCAL;
     bp->data = (void *)md;
     bp->operator_free = operator_mutate2Local_free;
     bp->operator_start = operator_mutate2Local_start;
     bp->operator_end = operator_mutate2Local_end;
     bp->operator_operate = operator_mutate2Local;

     /* default values for the mutation-specific data structure. */
     md->keep_trying = 0;
     md->internal = 0.9;
     md->external = 0.1;
     md->tree = (double *)MALLOC ( tree_count * sizeof ( double ) );
     for ( j = 0; j < tree_count; ++j )
          md->tree[j] = 0.0;
     md->treetotal = 0.0;
     md->sname = NULL;
     md->method = GENERATE_HALF_AND_HALF;
     md->mindepth = -1;
     md->maxdepth = -1;
     md->absDepth_czj = 0;

     j = parse_o_rama ( options, &argv );

     for ( i = 0; i < j; ++i )
     {
          if ( strcmp ( "keep_trying", argv[i] ) == 0 )
          {
               md->keep_trying = translate_binary ( argv[++i] );
               if ( md->keep_trying == -1 )
               {
                    ++errors;
                    error ( E_ERROR, "mutation: \"%s\" is not a valid setting for \"keep_trying\".",
                           argv[i] );
               }
          }
          else if ( strcmp ( "internal", argv[i] ) == 0 )
          {
               internalset = 1;
               md->internal = strtod ( argv[++i], NULL );
               if ( md->internal < 0.0 )
               {
                    ++errors;
                    error ( E_ERROR, "mutation: \"internal\" must be nonnegative." );
               }
          }
          else if ( strcmp ( "external", argv[i] ) == 0 )
          {
               externalset = 1;
               md->external = strtod ( argv[++i], NULL );
               if ( md->external < 0.0 )
               {
                    ++errors;
                    error ( E_ERROR, "mutation: \"external\" must be nonnegative." );
               }
          }
          else if ( strcmp ( "select", argv[i] ) == 0 )
          {
               if ( !exists_select_method ( argv[++i] ) )
               {
                    ++errors;
                    error ( E_ERROR, "mutation: \"%s\" is not a known selection method.",
                           argv[i] );
               }
               FREE ( md->sname );
               md->sname = (char *)MALLOC ( (strlen(argv[i])+1) * sizeof ( char ) );
               strcpy ( md->sname, argv[i] );
          }
          else if ( strcmp ( "method", argv[i] ) == 0 )
          {
               ++i;
               if ( strcmp ( argv[i], "half_and_half" ) == 0 )
                    md->method = GENERATE_HALF_AND_HALF;
               else if ( strcmp ( argv[i], "grow" ) == 0 )
                    md->method = GENERATE_GROW;
               else if ( strcmp ( argv[i], "full" ) == 0 )
                    md->method = GENERATE_FULL;
               else
               {
                    ++errors;
                    error ( E_ERROR, "mutation: \"%s\" is not a known generation method.",
                           argv[i] );
               }
          }
          else if ( strcmp ( "depth", argv[i] ) == 0 )
          {
               md->mindepth = strtol ( argv[++i], &cp, 10 );
               if ( *cp == 0 )
                    md->maxdepth = md->mindepth;
               else if ( *cp == '-' )
               {
                    md->maxdepth = strtol ( cp+1, &cp, 10 );
                    if ( *cp )
                    {
                         ++errors;
                         error ( E_ERROR, "mutation: malformed depth string \"%s\".",
                                argv[i] );
                    }
               }
               else
               {
                    ++errors;
                    error ( E_ERROR, "mutation: malformed depth string \"%s\".",
                           argv[i] );
               }
          }
          else if ( strcmp ( "tree", argv[i] ) == 0 )
          {
               k = parse_o_rama ( argv[++i], &targv );
               if ( k != tree_count )
               {
                    ++errors;
                    error ( E_ERROR, "mutation: wrong number of tree fields: \"%s\".",
                           argv[i] );
               }
               else
               {
                    for ( m = 0; m < k; ++m )
                    {
                         md->tree[m] = strtod ( targv[m], &cp );
                         if ( *cp )
                         {
                              ++errors;
                              error ( E_ERROR, "mutation: \"%s\" is not a number.",
                                     targv[m] );
                         }
                    }
               }
               
               free_o_rama ( k, &targv );
          }
          else if ( strncmp ( "tree", argv[i], 4 ) == 0 )
          {
               k = strtol ( argv[i]+4, &cp, 10 );
               if ( *cp )
               {
                    ++errors;
                    error ( E_ERROR, "mutation: unknown option \"%s\".",
                           argv[i] );
               }
               if ( k < 0 || k >= tree_count )
               {
                    ++errors;
                    error ( E_ERROR, "mutation: \"%s\" is out of range.",
                           argv[i] );
               }
               else
               {
                    md->tree[k] = strtod ( argv[++i], &cp );
                    if ( *cp )
                    {
                         ++errors;
                         error ( E_ERROR, "mutation: \"%s\" is not a number.",
                                argv[i] );
                    }
               }
          }

/************* begin czj changes */
          else if ( strcmp ( "depth_abs", argv[i] ) == 0 )
          {
               md->absDepth_czj = translate_binary ( argv[++i] );
               if ( md->absDepth_czj == -1 )
               { ++errors;
                 error(E_ERROR,
         "mutation: \"%s\" is not a valid setting for \"depth_abs\".",argv[i]);
               }
          }
/************* end czj changes */


          else
          {
               ++errors;
               error ( E_ERROR, "mutation: unknown option \"%s\".",
                      argv[i] );
          }
     }
     
     free_o_rama ( j, &argv );
     
     if ( internalset && !externalset )
          md->external = 0.0;
     else if ( !internalset && externalset )
          md->internal = 0.0;
     
     if ( md->sname == NULL )
     {
          ++errors;
          error ( E_ERROR, "mutation: no selection method specified." );
     }

     if ( md->mindepth == -1 && md->maxdepth == -1 )
     {
          md->mindepth = 0;
          md->maxdepth = 4;
     }
     
     if ( md->mindepth < 0 || md->maxdepth < 0 ||
         md->maxdepth < md->mindepth )
     {
          ++errors;
          error ( E_ERROR, "mutation: bad depth range.\n" );
     }
     
     for ( j = 0; j < tree_count; ++j )
          md->treetotal += md->tree[j];
     if ( md->treetotal == 0.0 )
     {
          for ( j = 0; j < tree_count; ++j )
               md->tree[j] = 1.0;
          md->treetotal = tree_count;
     }
          
     r = 0.0;
     for ( j = 0; j < tree_count; ++j )
          r = (md->tree[j] += r);

#ifdef DEBUG
     if ( !errors )
     {
          printf ( "mutation options:\n" );
          printf ( "   internal: %lf  external: %lf\n", md->internal, md->external );
          printf ( "   keep_trying: %d\n", md->keep_trying );
          printf ( "   selection: %s\n", md->sname==NULL?"NULL":md->sname );
          printf ( "   method: %d    mindepth: %d   maxdepth: %d\n",
                  md->method, md->mindepth, md->maxdepth );
          printf ( "   tree total: %lf\n", md->treetotal );
          for ( j = 0; j < tree_count; ++j )
               printf ( "   tree %d: %lf\n", j, md->tree[j] );
     }
#endif
     
     return errors;
}

/* operator_mutate_free()
 *
 * free mutation stuff.
 */

void operator_mutate2Local_free ( void *data )
{
     mutate_data * md;

     md = (mutate_data *)data;
     FREE ( md->sname );
     FREE ( md->tree );
     FREE ( md );
}

/* operator_mutate_start()
 *
 * get selection context for mutation operator.
 */

double *mutation_wheel;

int spin_mutation_wheel_lcl1_for_lvl1(int trueIndex,int *context,int fixed)
/*
   'context' is an array of children of the root 'trueIndex'
   'fixed' is the part that we are mutating from
*/

{
	int f;
	int NumFT = fset[0].function_count;
	
	// Convert our true index to index to our counters
	f = root_fset_index_to_counter_index(trueIndex);

#ifdef DEBUG_MUTATE
	printf("%d converted to lvl1 counter index of %d\n",trueIndex,f);
	if (Level1_Cnt_mwh.f == NULL)
		printf("Level 1 counters not allocated!\n");
#endif	

	int farity = Level1_Cnt_mwh.f[f].dim;

	// number of functions and terminals that will be on the wheel
	int size_fixed = Level1_Cnt_mwh.f[f].size[fixed];
	
	// Now we need to convert 'context' to our indices
	int i,j;
	for (i = 0; i < farity; i++) {
		if (i != fixed) {
			j = 0;
			while (context[i] != MS_czj[trueIndex][i].members[j])
				j++;
			// DEBUG
			if ( j > MS_czj[trueIndex][i].numFT) {
				printf("Couldnt convert fset context info to counter info\n");
				printf("in spin_mutation_wheel\n");
				exit(1);
			}
			context[i] = j;
		}
	}
	
	
	// Find the base offset not including our 'fixed' index
	int base_offset = 0;
	for (i = 0; i < farity; i++) {
		if (i != fixed)
			base_offset += context[i]*Level1_Cnt_mwh.f[f].offset[i];
	}

	mutation_wheel = (double *)calloc(size_fixed,sizeof(double));
	int this_offset = 0;
	for (i = 0; i < size_fixed; i++) {
		this_offset = base_offset +
			Level1_Cnt_mwh.f[f].offset[fixed]*i;
			
		mutation_wheel[i] = Level1_Cnt_mwh.f[f].prob[this_offset];
	}
	
	// Now turn this list of probabilities into a wheel we can spin
	for (i = 1; i < size_fixed; i++) {
		mutation_wheel[i] += mutation_wheel[i-1];
	}
	
	double spin = random_double() * mutation_wheel[size_fixed - 1];
	
	int winner = 0;
	while (spin > mutation_wheel[winner])
		winner++;
			
	// Convert 'winner', which is the index to fixed, to 
	// fset index
	int true_winner;
	true_winner = MS_czj[trueIndex][fixed].members[winner];
	
#ifdef DEBUG_MUTATE
	printf("Spinning wheel for arity %d\n",fixed);
	printf("Wheel is [");
	for (i = 0; i < size_fixed; i++) {
		printf("%f ",mutation_wheel[i]);
	}
	printf("]\n");
	printf("Wheel spun, spin=%f, winner=%d, true fset winner=%d\n",spin,winner,true_winner);
	
#endif
	free(mutation_wheel);
	return true_winner;
}


int spin_mutation_wheel_lcl1_for_lcl1(int trueIndex,int *context,int fixed)
/*
   'context' is an array of children of the root 'trueIndex'
   'fixed' is the child(a arity) that we are mutating from
*/

{
	int f;
	
	// NOT NEEDED FOR LOCAL
	//int NumFT = fset[0].function_count;
	//// Convert our true index to index to our counters
	//f = root_fset_index_to_counter_index(trueIndex);
	f = trueIndex;
	
#ifdef DEBUG_MUTATE
	printf("%d converted to lvl1 counter index of %d\n",trueIndex,f);
	if (Level1_Cnt_mwh.f == NULL)
		printf("Level 1 counters not allocated!\n");
#endif	

	int farity = Local1_mwh.f[f].dim;

	// number of functions and terminals that will be on the wheel
	int size_fixed = Local1_mwh.f[f].size[fixed];
	
	// Now we need to convert 'context' to our indices
	int i;
	for (i = 0; i < farity; i++) {
		if (i != fixed) {
			context[i] = MS_czj[f][i].convert_to_ms[ context[i] ];
		}
	}
	
	
	// Find the base offset not including our 'fixed' index
	int base_offset = 0;
	for (i = 0; i < farity; i++) {
		if (i != fixed)
			base_offset += context[i]*Local1_mwh.f[f].offset[i];
	}

	mutation_wheel = (double *)calloc(size_fixed,sizeof(double));
	int this_offset = 0;
	for (i = 0; i < size_fixed; i++) {
		this_offset = base_offset +
			Local1_mwh.f[f].offset[fixed]*i;
			
		mutation_wheel[i] = Local1_mwh.f[f].prob[this_offset];
	}
	
	// Now turn this list of probabilities into a wheel we can spin
	for (i = 1; i < size_fixed; i++) {
		mutation_wheel[i] += mutation_wheel[i-1];
	}
	
	double spin = random_double() * mutation_wheel[size_fixed - 1];
	
	int winner = 0;
	while (spin > mutation_wheel[winner])
		winner++;
			
	// Convert 'winner', which is the index to fixed, to 
	// fset index
	int true_winner;
	true_winner = MS_czj[f][fixed].members[winner];
	
#ifdef DEBUG_MUTATE
	printf("Spinning local info mutation wheel for arity %d\n",fixed);
	printf("Wheel is [");
	for (i = 0; i < size_fixed; i++) {
		printf("%f ",mutation_wheel[i]);
	}
	printf("]\n");
	printf("Wheel spun, spin=%f, winner=%d, true fset winner=%d\n",spin,winner,true_winner);
	
#endif
	free(mutation_wheel);
	return true_winner;
}

int spin_mutation_wheel_lcl1_for_lcl1_terminal(int trueIndex,int *context,int fixed)
/*
   'context' is an array of children of the root 'trueIndex'
   'fixed' is the child(a arity) that we are mutating from
*/

{
	int f;
	
	// NOT NEEDED FOR LOCAL
	//int NumFT = fset[0].function_count;
	//// Convert our true index to index to our counters
	//f = root_fset_index_to_counter_index(trueIndex);
	f = trueIndex;
	
#ifdef DEBUG_MUTATE
	printf("%d converted to lvl1 counter index of %d\n",trueIndex,f);
	if (Level1_Cnt_mwh.f == NULL)
		printf("Level 1 counters not allocated!\n");
#endif	

	int farity = Local1_mwh.f[f].dim;

	// number of functions and terminals that will be on the wheel
	int size_fixed = MS_czj[trueIndex][fixed].numT;
	// int size_fixed = Local1_mwh.f[f].size[fixed];
	int start_terminals = MS_czj[trueIndex][fixed].numF;
	
	// Now we need to convert 'context' to our indices
	int i;
	for (i = 0; i < farity; i++) {
		if (i != fixed) {
			context[i] = MS_czj[f][i].convert_to_ms[ context[i] ];
		}
	}
	
	
	// Find the base offset not including our 'fixed' index
	int base_offset = 0;
	for (i = 0; i < farity; i++) {
		if (i != fixed)
			base_offset += context[i]*Local1_mwh.f[f].offset[i];
	}

	mutation_wheel = (double *)calloc(size_fixed,sizeof(double));
	int this_offset = 0;
	for (i = 0; i < size_fixed; i++) {
		this_offset = base_offset +
			Local1_mwh.f[f].offset[fixed]*(i+start_terminals);
			
		mutation_wheel[i] = Local1_mwh.f[f].prob[this_offset];
	}
	
	// Now turn this list of probabilities into a wheel we can spin
	for (i = 1; i < size_fixed; i++) {
		mutation_wheel[i] += mutation_wheel[i-1];
	}
	
	double spin = random_double() * mutation_wheel[size_fixed - 1];
	
	int winner = 0;
	
	// now this wheel starts at 0 but represents terminals past that
	while (spin > mutation_wheel[winner])
		winner++;
	
	// so now increment winner to the correct one
	winner += start_terminals;
			
	// Convert 'winner', which is the index to fixed, to 
	// fset index
	int true_winner;
	true_winner = MS_czj[f][fixed].members[winner];
	
#ifdef DEBUG_MUTATE
	printf("Spinning local info mutation wheel for arity %d\n",fixed);
	printf("Only want terminal here in spin_mutation_wheel_lcl1_for_lcl1_terminal\n");
	printf("Wheel is [");
	for (i = 0; i < size_fixed; i++) {
		printf("%f ",mutation_wheel[i]);
	}
	printf("]\n");
	printf("Wheel spun, spin=%f, winner=%d, true fset winner=%d\n",spin,winner,true_winner);
	
#endif
	free(mutation_wheel);
	return true_winner;
}

int spin_mutation_wheel_lcl1_for_lcl1_function(int trueIndex,int *context,int fixed)
/*
   'context' is an array of children of the root 'trueIndex'
   'fixed' is the child(a arity) that we are mutating from
*/

{
	int f;
	
	// NOT NEEDED FOR LOCAL
	//int NumFT = fset[0].function_count;
	//// Convert our true index to index to our counters
	//f = root_fset_index_to_counter_index(trueIndex);
	f = trueIndex;
	
#ifdef DEBUG_MUTATE
	printf("%d converted to lvl1 counter index of %d\n",trueIndex,f);
	if (Level1_Cnt_mwh.f == NULL)
		printf("Level 1 counters not allocated!\n");
#endif	

	int farity = Local1_mwh.f[f].dim;

	// number of functions and terminals that will be on the wheel
	int size_fixed = MS_czj[trueIndex][fixed].numF;
	
	// Now we need to convert 'context' to our indices
	int i;
	for (i = 0; i < farity; i++) {
		if (i != fixed) {
			context[i] = MS_czj[f][i].convert_to_ms[ context[i] ];
		}
	}
	
	
	// Find the base offset not including our 'fixed' index
	int base_offset = 0;
	for (i = 0; i < farity; i++) {
		if (i != fixed)
			base_offset += context[i]*Local1_mwh.f[f].offset[i];
	}

	mutation_wheel = (double *)calloc(size_fixed,sizeof(double));
	int this_offset = 0;
	for (i = 0; i < size_fixed; i++) {
		this_offset = base_offset +
			Local1_mwh.f[f].offset[fixed]*i;
			
		mutation_wheel[i] = Local1_mwh.f[f].prob[this_offset];
	}
	
	// Now turn this list of probabilities into a wheel we can spin
	for (i = 1; i < size_fixed; i++) {
		mutation_wheel[i] += mutation_wheel[i-1];
	}
	
	double spin = random_double() * mutation_wheel[size_fixed - 1];
	
	int winner = 0;
	
	while (spin > mutation_wheel[winner])
		winner++;
			
	// Convert 'winner', which is the index to fixed, to 
	// fset index
	int true_winner;
	true_winner = MS_czj[f][fixed].members[winner];
	
#ifdef DEBUG_MUTATE
	printf("Spinning local info mutation wheel for arity %d\n",fixed);
	printf("Only want function here in spin_mutation_wheel_lcl1_for_lcl1_function\n");
	printf("Wheel is [");
	for (i = 0; i < size_fixed; i++) {
		printf("%f ",mutation_wheel[i]);
	}
	printf("]\n");
	printf("Wheel spun, spin=%f, winner=%d, true fset winner=%d\n",spin,winner,true_winner);
	
#endif
	free(mutation_wheel);
	return true_winner;
}

void operator_mutate2Local_start ( population *oldpop, void *data )
{
     mutate_data * md;
     select_context_func_ptr select_con;

     md = (mutate_data *)data;
     select_con = get_select_context ( md->sname );
     md->sc = select_con ( SELECT_INIT, NULL, oldpop, md->sname );
}

/* operator_mutate_end()
 *
 * free selection context for mutation operator.
 */

void operator_mutate2Local_end ( void *data )
{
     mutate_data * md;

     md = (mutate_data *)data;
     md->sc->context_method ( SELECT_CLEAN, md->sc, NULL, NULL );
}

/* operator_mutate()
 *
 * do the mutation.
 */

void operator_mutate2Local ( population *oldpop, population *newpop,
                      void *data )
{
     int i,j;
     int ps;
     lnode *replace[2];
     int l, ns;
     int badtree;
     int repcount;
     mutate_data * md;
     int t;
     double r;
     int depth;
     int totalnodes;
     int p;
     int forceany;
     double total;
     
     
     md = (mutate_data *)data;
     total = md->internal + md->external;

     /* choose a tree to mutate. */
     r = random_double() * md->treetotal;
     for ( t = 0; r >= md->tree[t]; ++t );

     /* select an individual to mutate. */
     p = md->sc->select_method ( md->sc ); 
     ps = tree_nodes ( oldpop->ind[p].tr[t].data );
     forceany = (ps==1||total==0.0);

#ifdef DEBUG_MUTATE
	fprintf(stderr,"ENTERING MUTATION\n");
	fflush(stdout);
	fflush(stderr);
     fprintf ( stderr, "the parent size is %d\n", ps );
     fprintf ( stderr, "    parent %4d: ", p );
     print_tree ( oldpop->ind[p].tr[t].data, stderr );
#endif          

     while(1)
     {

	  if ( forceany )
	  {
	       /* choose any point. */
	       l = random_int ( ps );
	       replace[0] = get_subtree ( oldpop->ind[p].tr[t].data, l );
	  }
	  else if ( total*random_double() < md->internal )
	  {
	       /* choose an internal point. */
	       l = random_int ( tree_nodes_internal ( oldpop->ind[p].tr[t].data ) );
	       replace[0] = get_subtree_internal ( oldpop->ind[p].tr[t].data, l );
	  }
	  else
	  {
	       /* choose an external point. */
	       l = random_int ( tree_nodes_external ( oldpop->ind[p].tr[t].data ) );
	       replace[0] = get_subtree_external ( oldpop->ind[p].tr[t].data, l );
	  }
	  
	  // We need to traverse the tree to gather context information
	  int replace_depth = 0;
	  
	  int trueIndex;
	  int newNode;
	  int thisArity;
	  gather_info_mwh(oldpop->ind[p].tr[t].data,replace[0],&replace_depth,&thisArity,&trueIndex,level1context_mwh);
	 
	 	  
#ifdef DEBUG_MUTATE
          fprintf(stderr,"subtree number: %d\n",l);         /* added by czj */
          fprintf ( stderr, "selected for replacement: " );
          print_tree ( replace[0], stderr );
	  printf("Context information for mutation\n");
	  printf("Depth of our mutation is %d\n",replace_depth);
	  if (replace_depth > 0)
	  {
	  printf("Context information is (%s)(",fset[0].cset[trueIndex].string);
	  for (i = 0; i < fset[0].cset[trueIndex].arity; i++) {
	  	printf("%s ",fset[0].cset[level1context_mwh[i]].string);
	  }
	  printf(")\n");
	  }
	  else 
	  {
	  	printf("No context information since target found at root\n");
	  }
#endif
	  // if mutating at the root
	  if (replace_depth == 0) {
	  	// In this case, we want to regrow from the root
		// so we want to spin a wheel for entire first level,
		// like our regrow operator, so lets just call our regrow operator and do it
	  	
		gensp_reset ( 1 );
	 	 /* pick a value from the depth ramp. */
          	depth = md->mindepth + random_int ( md->maxdepth - md->mindepth + 1 );
		if (md->absDepth_czj)
			MinDepth_czj = depth-md->mindepth;
	  	/* grow the tree. */
          	switch ( md->method )
          	{
             	case GENERATE_GROW:
               		generate_random_grow_tree_wrapper_lcl1_mwh ( 1, depth, fset+tree_map[t].fset,REGROW_SOURCE );
               		break;
                case GENERATE_FULL:
               		generate_random_full_tree_wrapper_lcl1_mwh ( 1, depth, fset+tree_map[t].fset,REGROW_SOURCE );
               		break;
                case GENERATE_HALF_AND_HALF:
               		if ( random_double() < 0.5 )
                    		generate_random_grow_tree_wrapper_lcl1_mwh ( 1, depth, fset+tree_map[t].fset,REGROW_SOURCE );  // calls mwh version for lvl 1 information
               		else
                    		generate_random_full_tree_wrapper_lcl1_mwh ( 1, depth, fset+tree_map[t].fset,REGROW_SOURCE );  // calls wrapper for lvl 1 information
               		break;
          	}
		
	  
	  }
	  
	  // if mutating a child at lvl 1
	  // we want to just do it at lvl 1 like it is normal mutate, then
	  // regrow from that point on.
	  if (replace_depth == 1) {
#ifdef DEBUG_MUTATE
	  	printf("Doing level one mutation on arity %d of our root\n",thisArity);
	  	printf("So mutating at '%s'\n",fset[0].cset[level1context_mwh[thisArity]].string);
#endif
		gensp_reset ( 1 );
		
		depth = md->mindepth + random_int ( md->maxdepth - md->mindepth + 1);
		if (md->absDepth_czj)
			MinDepth_czj = depth-md->mindepth;
	  	
		switch(md->method) 
		{
			case GENERATE_GROW:
				if (depth <= 0)
					i = random_T_czj();
				else if (depth > MinDepth_czj)
					i = random_F_czj();
				else
					i = random_FT_czj();
			break;
			case GENERATE_FULL:
				if (depth <= 0)
					i = random_T_czj();
				else
					i = random_F_czj();
			break;
			case GENERATE_HALF_AND_HALF:
				if (random_double() < 0.5) {
					if (depth <= 0)
						i = random_T_czj();
					else if (depth > MinDepth_czj)
						i = random_F_czj();
					else
						i = random_FT_czj();
				}
				else {
					if (depth <= 0)
						i = random_T_czj();
					else
						i = random_F_czj();
				}
			break;
		}
		
		// printf("The winning lvl1 mutate is %d\n",i);
		int save;
		int space = 1;
		
		// Lets generate our space for the root
		gensp_next(space)->f = (fset->cset)+i;
		
	        // Now, if this is a function, I want to spin a wheel to generate its children
		int *child_list_lcl;
		
		if (i < fset->function_count)
                {		
     	             // now spin local information to find all children of 'i'
     	             child_list_lcl = (int *)MALLOC( (fset->cset)[i].arity * sizeof(int));
		     switch(md->method)
		     {
		     case GENERATE_GROW:
	             	if (depth <= 1)
		          	spin_regrow_lcl1_information_terminals_mwh(i,child_list_lcl);
	             	else if (depth > MinDepth_czj + 1)
     		          	spin_regrow_lcl1_information_functions_mwh(i,child_list_lcl);
			  else
			  	spin_regrow_lcl1_information_mwh(i,child_list_lcl);
			break;
		     case GENERATE_FULL:
		        if (depth <= 1)
				spin_regrow_lcl1_information_terminals_mwh(i,child_list_lcl);
			else
				spin_regrow_lcl1_information_functions_mwh(i,child_list_lcl);
			break;
			
		     case GENERATE_HALF_AND_HALF:
		        if (random_double() < 0.5) {
				if (depth <= 1)
		          		spin_regrow_lcl1_information_terminals_mwh(i,child_list_lcl);
	             		else if (depth > MinDepth_czj + 1)
     		          		spin_regrow_lcl1_information_functions_mwh(i,child_list_lcl);
			  	else
			  		spin_regrow_lcl1_information_mwh(i,child_list_lcl);
				}
			else {
				if (depth <= 1)
					spin_regrow_lcl1_information_terminals_mwh(i,child_list_lcl);
				else
					spin_regrow_lcl1_information_functions_mwh(i,child_list_lcl);
			}
			
			break;	
		    }		      
                }
		
		
		if (i>= fset->function_count) {
			if ((fset->cset)[i].ephem_gen)
				gensp_next(space)->d = new_ephemeral_const((fset->cset)+i);
				
		}
		else {
			switch( (fset->cset)[i].type)
				{ case FUNC_DATA: case EVAL_DATA:
					for (j=0; j < (fset->cset)[i].arity; ++j)
					{ Function_czj = i;
					  Argument_czj = j;
					  switch ( md->method )
          				  {
             	                          case GENERATE_GROW:
					       generate_random_grow_tree_lcl1_mwh(space,depth-1,fset,child_list_lcl[j]);
					       break;
					  case GENERATE_FULL:
					       generate_random_full_tree_lcl1_mwh(space,depth-1,fset,child_list_lcl[j]);
					       break;     
					  case GENERATE_HALF_AND_HALF:
               					if ( random_double() < 0.5 )   
						     generate_random_grow_tree_lcl1_mwh(space,depth-1,fset,child_list_lcl[j]);
						else
						     generate_random_full_tree_lcl1_mwh(space,depth-1,fset,child_list_lcl[j]);
						break;
					  }  
					}
					break;
				   case FUNC_EXPR: case EVAL_EXPR:
		   			for (j = 0; j < (fset->cset)[i].arity; ++j)
					{	Function_czj = i;
					Argument_czj = j;
					save = gensp_next_int(space);
					
					switch ( md->method )
          				  {
             	                          case GENERATE_GROW:
					       generate_random_grow_tree_lcl1_mwh(space,depth-1,fset,child_list_lcl[j]);
					       break;
					  case GENERATE_FULL:
					       generate_random_full_tree_lcl1_mwh(space,depth-1,fset,child_list_lcl[j]);
					       break;     
					  case GENERATE_HALF_AND_HALF:
               					if ( random_double() < 0.5 )   
						     generate_random_grow_tree_lcl1_mwh(space,depth-1,fset,child_list_lcl[j]);
						else
						     generate_random_full_tree_lcl1_mwh(space,depth-1,fset,child_list_lcl[j]);
						break;
					  }  
										
					gensp[space].data[save].s = gensp[space].used-save-1;
					}
					break;
		
		
				}
		
			}
	  }

	  // if this is going to be a local mutate
	  if (replace_depth > 1) 
	  {
	      	  gensp_reset ( 1 );
		
		  /* pick a value from the depth ramp. */
		    // This picks how deep we want to grow, this depth IS NOT the depth
		
	  	  depth = md->mindepth + random_int ( md->maxdepth - md->mindepth + 1);
		  if (md->absDepth_czj)
		  	MinDepth_czj = depth-md->mindepth;
	  
	  	  // Spin our local mutation wheel
		  // depending if we want a terminal or not
		  
		  switch( md->method)
		  {
		  	case GENERATE_GROW:
		  		if (depth <= 0)
		  			i = spin_mutation_wheel_lcl1_for_lcl1_terminal(trueIndex,level1context_mwh,thisArity);
		  		else if (depth>MinDepth_czj)
		  	  	        i = spin_mutation_wheel_lcl1_for_lcl1(trueIndex,level1context_mwh,thisArity);  
				else
					i = spin_mutation_wheel_lcl1_for_lcl1_function(trueIndex,level1context_mwh,thisArity);
		  		break;
				
			case GENERATE_FULL:
				if (depth <= 0)
					i = spin_mutation_wheel_lcl1_for_lcl1_terminal(trueIndex,level1context_mwh,thisArity);
				else
					i = spin_mutation_wheel_lcl1_for_lcl1_function(trueIndex,level1context_mwh,thisArity);
				break;
			
			case GENERATE_HALF_AND_HALF:
				if (random_double() < 0.5) {
					if (depth <= 0)
		  				i = spin_mutation_wheel_lcl1_for_lcl1_terminal(trueIndex,level1context_mwh,thisArity);
		  			else if (depth>MinDepth_czj)
		  	  		        i = spin_mutation_wheel_lcl1_for_lcl1(trueIndex,level1context_mwh,thisArity);  
					else
						i = spin_mutation_wheel_lcl1_for_lcl1_function(trueIndex,level1context_mwh,thisArity);
				}
				else {
					if (depth <= 0)
						i = spin_mutation_wheel_lcl1_for_lcl1_terminal(trueIndex,level1context_mwh,thisArity);
					else
						i = spin_mutation_wheel_lcl1_for_lcl1_function(trueIndex,level1context_mwh,thisArity);
				}
				break;
		  }
		  
		  //printf("The winning local mutate is %d\n",i);
		  int save;
		  int space = 1;
		  
		  int *child_list_lcl;
		  // Now lets generate children for it if we have to
		  if (i < fset->function_count)
		  {
		  	
						
		  	child_list_lcl = (int *)MALLOC( (fset->cset)[i].arity * sizeof(int));
			// now spin local information to find all children of 'i'
			
			switch (md->method)
			{
				case GENERATE_GROW:
					if (depth <= 1)
						spin_regrow_lcl1_information_terminals_mwh(i,child_list_lcl);
					else if (depth > MinDepth_czj + 1)
						spin_regrow_lcl1_information_mwh(i,child_list_lcl);
					else
						spin_regrow_lcl1_information_functions_mwh(i,child_list_lcl);
				break;
				case GENERATE_FULL:
					if (depth <= 1)
						spin_regrow_lcl1_information_terminals_mwh(i,child_list_lcl);
					else
						spin_regrow_lcl1_information_functions_mwh(i,child_list_lcl);
				break;
				case GENERATE_HALF_AND_HALF:
					if (random_double() < 0.5) {
						if (depth <= 1)
							spin_regrow_lcl1_information_terminals_mwh(i,child_list_lcl);
						else if (depth > MinDepth_czj + 1)
							spin_regrow_lcl1_information_mwh(i,child_list_lcl);
						else
							spin_regrow_lcl1_information_functions_mwh(i,child_list_lcl);
					}
					else {
						if (depth <= 1)
							spin_regrow_lcl1_information_terminals_mwh(i,child_list_lcl);
						else
							spin_regrow_lcl1_information_functions_mwh(i,child_list_lcl);
					}
			}
		  
		  }
		  
		  // save space for our new item
                  gensp_next(space)->f = (fset->cset)+i;
		  		  
		  
		  if (i>= fset->function_count) 
		  {
			if ((fset->cset)[i].ephem_gen)
		 		gensp_next(space)->d = new_ephemeral_const((fset->cset)+i);
				
		  }
		  else 
		  {
			
			switch( (fset->cset)[i].type)
				{ case FUNC_DATA: case EVAL_DATA:
					for (j=0; j < (fset->cset)[i].arity; ++j)
					{ Function_czj = i;
					  Argument_czj = j;
					  // Need to go depending on type of trees
					  // we grow(grow/full/halfnhalf
					  switch ( md->method )
	          				{
       		      			        case GENERATE_GROW:
               				          generate_random_grow_tree_lcl1_mwh (space, depth-1, fset+tree_map[t].fset,child_list_lcl[j] );
               		              		  break;
             	 				case GENERATE_FULL:
                                                  generate_random_full_tree_lcl1_mwh ( space, depth-1, fset+tree_map[t].fset,child_list_lcl[j] );
               		                          break;
                                                case GENERATE_HALF_AND_HALF:
               		                          if ( random_double() < 0.5 )
                    		                      generate_random_grow_tree_lcl1_mwh(space,depth-1,fset+tree_map[t].fset,child_list_lcl[j] );
               		                          else
                    		                      generate_random_full_tree_lcl1_mwh(space,depth-1,fset+tree_map[t].fset,child_list_lcl[j] );
               		                          break;
          	                            }
					  
					  
					    //generate_random_full_tree(space,depth-1,fset+tree_map[t].fset);
					}
					break;
				   case FUNC_EXPR: case EVAL_EXPR:
		   			for (j = 0; j < (fset->cset)[i].arity; ++j)
					{	Function_czj = i;
					Argument_czj = j;
					save = gensp_next_int(space);
					// Need to go depending on type of trees
					// we grow(grow/full/halfnhalf
					switch ( md->method )
	          			{
       		      			        case GENERATE_GROW:
               				          generate_random_grow_tree_lcl1_mwh (space, depth-1, fset+tree_map[t].fset,child_list_lcl[j] );
               		              		  break;
             	 				case GENERATE_FULL:
                                                  generate_random_full_tree_lcl1_mwh ( space, depth-1, fset+tree_map[t].fset,child_list_lcl[j] );
               		                          break;
                                                case GENERATE_HALF_AND_HALF:
               		                          if ( random_double() < 0.5 )
                    		                      generate_random_grow_tree_lcl1_mwh ( space, depth-1,
						      fset+tree_map[t].fset,child_list_lcl[j] );
               		                          else
                    		                      generate_random_full_tree_lcl1_mwh ( space, depth-1,
						      fset+tree_map[t].fset,child_list_lcl[j] );
               		                          break;
          	                        }
					
					// generate_random_full_tree(space,depth-1,fset+tree_map[t].fset);
					gensp[space].data[save].s = gensp[space].used-save-1;
					}
					break;
		
		
				}
		
		}
		
		MinDepth_czj = 0;
		if (i < fset->function_count)
			FREE(child_list_lcl);
	}
	  
	  
#ifdef DEBUG_MUTATE
          fprintf ( stderr, "the new subtree is: " );
          print_tree ( gensp[1].data, stderr );
#endif

	  /* count the nodes in the new tree. */
          ns = ps - tree_nodes ( replace[0] ) + tree_nodes ( gensp[1].data );
          totalnodes = ns;

	  /* check the mutated tree against node count and/or size limits. */
          badtree = 0;
          if ( tree_map[t].nodelimit > -1 && ns > tree_map[t].nodelimit )
               badtree = 1;
          else if ( tree_map[t].depthlimit > -1 )
          {
               ns = tree_depth_to_subtree ( oldpop->ind[p].tr[t].data,
                                           replace[0] ) +
                    tree_depth ( gensp[1].data );
               if ( ns > tree_map[t].depthlimit )
                    badtree = 1;
          }

/************ begin czj changes */
          else
            if (md->absDepth_czj && tree_depth( gensp[1].data ) < md->mindepth)
              badtree = 1;
/************ end czj changes */

	  /* if tree is too big and keep_trying is set, then skip to the
	     stop and choose a new mutation point/mutant subtree. */
          if ( md->keep_trying && badtree )
               continue;

	  /* check mutated tree against whole-individual node limits. */
          if ( ind_nodelimit > -1 )
          {
               for ( i = 0; i < tree_count; ++i )
                    if ( i != t )
                         totalnodes += oldpop->ind[p].tr[i].nodes;
               badtree |= (totalnodes > ind_nodelimit);
          }

	  /* if tree is too big and keep_trying is set, then skip to the
	     stop and choose a new mutation point/mutant subtree. */
          if ( md->keep_trying && badtree )
               continue;

          if ( badtree )
          {
#ifdef DEBUG_MUTATE
               fprintf ( stderr,
                        "new tree is too big; reproducing parent.\n" );
#endif
	       /* tree too big but keep_trying not set, just reproduce
		  parent tree. */
               duplicate_individual ( (newpop->ind)+newpop->next,
                                     (oldpop->ind)+p );
          }
          else
          {
#ifdef DEBUG_MUTATE
               fprintf ( stderr, "new tree is permissible.\n" );
#endif
	       /* copy the parent tree to the offspring position. */
               duplicate_individual ( (newpop->ind)+newpop->next,
                                     (oldpop->ind)+p );
	       /* free the tree selected for mutation. */
               free_tree ( newpop->ind[newpop->next].tr+t );
               
               /* copy the selected tree, replacing the subtree at the
		  mutation point with the randomly generated tree. */
               replace[1] = gensp[1].data;
               copy_tree_replace_many ( 0, oldpop->ind[p].tr[t].data,
                                       replace, replace+1, 1, &repcount );
               if ( repcount != 1 )
               {
                    error ( E_FATAL_ERROR,
                           "botched mutation:  this can't happen." );
               }
#if VERIFY_MUTATION
  if (verify_tree_czj(oldpop->[p].tr[t].data))
  { oprintf(OUT_SYS, 10, "INVALID TREE in mutation\n");
    exit(1);
  }
#endif
	       /* copy the tree to the new individual. */
               gensp_dup_tree ( 0, newpop->ind[newpop->next].tr+t );
               newpop->ind[newpop->next].evald = EVAL_CACHE_INVALID;
               newpop->ind[newpop->next].flags = FLAG_NONE;
               
#ifdef DEBUG_MUTATE
               fprintf ( stderr, "    the mutated individual is: " );
               print_individual ( newpop->ind+newpop->next, stderr );
#endif
          }

          ++newpop->next;
          break;
     }

#ifdef DEBUG_MUTATE
     printf ( "MUTATION COMPLETE.\n\n\n" );
#endif
}
