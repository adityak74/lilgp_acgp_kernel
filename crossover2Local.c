
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
extern typeLevel1_Cnt Local1_mwh;   // global local counter information
extern int *level1context_mwh;	// list of the indices of lvl 1 context data
extern int *level2context_mwh;  // list of indices of lvl 2 context data
extern int size_context1_mwh;   // number of elements in level1context_mwh
extern int size_context2_mwh;   // number of elements in level2context_mwh
extern int acgp_level;

/********************************************************************/
/********************************************************************
   crossover2Local operator
*********************************************************************/
/********************************************************************/



typedef struct
{
     int keep_trying;
     double internal;
     double external;
     double *tree;       /* probability that a given tree
			    will be selected for crossover. */
     double *treecumul;  /* running sum of "tree" field. */
     double treetotal;   /* total of all tree fields. */
     double *func;       /* probability that a given function
			    set will be selected for crossover. */
     char *sname;
     sel_context *sc;
     char *sname2;
     sel_context *sc2;
} crossover_data;


double *crossover_weights_mwh;

// void generate_crossover_weights_lcl1_for_lvl1_mwh(int trueIndex,int *context,int fixed)
/*
   'context' is an array of children of the root 'trueIndex'
   'fixed' is the part that we are crossing at, basically fixed in that
   we do not consider it as part of our context
*/
/*
{
	int f;
	int NumF = fset[0].function_count;
	
	// Convert our true index to index to our counters
	f = root_fset_index_to_counter_index(trueIndex);

#ifdef DEBUG_CROSSOVER
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
				printf("in generate_crossover_weights_mwh\n");
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
	double *weights_for_this_MS;
	weights_for_this_MS = (double *)calloc(size_fixed,sizeof(double));
	int this_offset = 0;
	for (i = 0; i < size_fixed; i++) {
		this_offset = base_offset +
			Level1_Cnt_mwh.f[f].offset[fixed]*i;
			
		weights_for_this_MS[i] = Level1_Cnt_mwh.f[f].prob[this_offset];
	}
	
	// Now weights_for_this_MS is the weights for each of the functions in this mutation set BUT
	// we want a weights list for ALL possible functions/terminals in fset, so lets give the others
	// weights of 0
	int numFT = NumF + fset[0].terminal_count;
	crossover_weights_mwh = (double *)calloc(numFT,sizeof(double));
	
	int c = 0; // counter of what number function we are in the MS for this context
	for (i = 0; i < numFT; i++) {
		// check to see if this label is part of our MS, ie: does it have a weight
		if (i == MS_czj[trueIndex][fixed].members[c] ) {
			// set it to our weight
			crossover_weights_mwh[i] = weights_for_this_MS[c];
			c++;  // move to next member of our mutation set
		}
		else {
			// Wasnt in our possible mutation set, so its weight is 0
			crossover_weights_mwh[i] = 0;
		}
	
	}
	
	free(weights_for_this_MS);
}
*/	

void generate_crossover_weights_lcl1_for_lcl1_mwh(int trueIndex,int *context,int fixed)
/*
   'context' is an array of children of the root 'trueIndex'
   'fixed' is the part that we are crossing at, basically fixed in that
   we do not consider it as part of our context
*/

{
	int f;
	// int NumF = fset[0].function_count;
	
	// Convert our true index to index to our counters
	// f = root_fset_index_to_counter_index(trueIndex);
	f = trueIndex;
	
#ifdef DEBUG_CROSSOVER
	printf("Parent of our context is %d,fixed index is %d\n",f,fixed);
	if (Local1_mwh.f == NULL)
		printf("Level 1 local counters not allocated!\n");
	fflush(stdout);
#endif	

	int farity = Local1_mwh.f[f].dim;

	// number of functions and terminals that will be on the wheel
	int size_fixed = Local1_mwh.f[f].size[fixed];
	
	// Now we need to convert 'context' to our indices
	int i,j;
	for (i = 0; i < farity; i++) {
		if (i != fixed) {
			j = 0;
			while (context[i] != MS_czj[f][i].members[j])
				j++;
			// DEBUG
			if ( j > MS_czj[f][i].numFT) {
				printf("Couldnt convert fset context info to counter info\n");
				printf("in generate_crossover_weights_mwh\n");
				exit(1);
			}
			context[i] = j;
		}
	}
	
	
	// Find the base offset not including our 'fixed' index
	int base_offset = 0;
	for (i = 0; i < farity; i++) {
		if (i != fixed)
			base_offset += context[i]*Local1_mwh.f[f].offset[i];
	}
	double *weights_for_this_MS;
	weights_for_this_MS = (double *)calloc(size_fixed,sizeof(double));
	int this_offset = 0;
	for (i = 0; i < size_fixed; i++) {
		this_offset = base_offset +
			Local1_mwh.f[f].offset[fixed]*i;
			
		weights_for_this_MS[i] = Local1_mwh.f[f].prob[this_offset];
	}
	
	// Now weights_for_this_MS is the weights for each of the functions in this mutation set BUT
	// we want a weights list for ALL possible functions/terminals in fset, so lets give the others
	// weights of 0
	int numFT = fset[0].function_count + fset[0].terminal_count;
	crossover_weights_mwh = (double *)calloc(numFT,sizeof(double));
	
/*	// BUG IN THIS CODE
	printf("Making weights for crossover\n");
	int c = -1; // counter of what number function we are in the MS for this context
	printf("Doing it for parent of %d,arity of %d\n",f,fixed);
	for (i = 0; i < numFT; i++) {
		// check to see if this label is part of our MS, ie: does it have a weight
		printf("i=%d,c=%d\n",i,c);
		if (i == MS_czj[f][fixed].members[c+1] ) {
			// set it to our weight
			c++;  
			crossover_weights_mwh[i] = weights_for_this_MS[c];
		}
		else {
			// Wasnt in our possible mutation set, so its weight is 0
			crossover_weights_mwh[i] = 0;
		}
	
	}
*/
	
	// Lets try a new way
	// This time start by setting them all to zero
	// Then set the nonzero ones
	for (i = 0; i < numFT; i++) {
		crossover_weights_mwh[i] = 0;
	}
	
	for (i = 0; i < MS_czj[f][fixed].numFT; i++) {
		crossover_weights_mwh[ MS_czj[f][fixed].members[i] ] = weights_for_this_MS[i];
	}
	
	free(weights_for_this_MS);
}	



int spin_after_generating_crossover_weights_lcl1_mwh(lnode *l, int size_tree,int trueIndex, int fixed) 
/*
    Using the global crossover_weights_mwh array, traverses the source tree and generates a wheel to
    replace the child of the root in position "fixed"
    
    trueIndex is index of the root
    
    Then spins the wheel and returns the winner
*/
{
	int f;
	int NumFT = fset[0].function_count;
	
	// Convert our true index to index to our counters
	//f = root_fset_index_to_counter_index(trueIndex);
	
	// Local counters so no reason to convert
	// since all possible parents will have a hypercube
	f = trueIndex;
	
        int i;
        double *crossover_wheel;
	crossover_wheel = (double *)calloc(size_tree,sizeof(double));
	
	// generate our crossover wheel
	int cur_node=0;
        generate_crossover_wheel_recurse_lcl1(&l,crossover_wheel,&cur_node,trueIndex,fixed);

	// Now turn this list of probabilities into a wheel we can spin
	for (i = 1; i < size_tree; i++) {
		crossover_wheel[i] += crossover_wheel[i-1];
	}
	
	double spin = random_double() * crossover_wheel[size_tree - 1];
	
	int winner = 0;
	while (spin > crossover_wheel[winner])
		winner++;
			
	// Convert 'winner', which is the index to fixed, to 
	// fset index
	//int true_winner;
	// true_winner = MS_czj[trueIndex][fixed].members[winner];
	
#ifdef DEBUG_CROSSOVER
	printf("Spinning local wheel for arity %d\n",fixed);
	printf("Wheel is [");
	for (i = 0; i < size_tree; i++) {
		printf("%f ",crossover_wheel[i]);
	}
	printf("]\n");
	printf("Wheel spun, spin=%f, subtree winner of source is=%d\n",spin,winner);
	
#endif
	free(crossover_weights_mwh);
	free(crossover_wheel);
	return winner;
}

int generate_crossover_wheel_recurse_lcl1( lnode **l,double *crossover_wheel, int *cur_node , int trueIndex, int fixed)
{
     function *parent;
     
     parent = (**l).f;
     int parentIndex = (**l).f->index;
     int i, j;
    
     // based on 'parent', lets update crossover_wheel[cur_node] based on our weights array crossover_weights_mwh
    // printf("Parent index is %d\n",parentIndex);
     double weight = crossover_weights_mwh[parentIndex]; 
     //printf("cur_node=%d\n",*cur_node);
     crossover_wheel[*cur_node] = weight;
	
#ifdef DEBUG_CROSSOVER
	//printf("Node=%d fsetIndex=%d, weight=%f",*cur_node,parentIndex,weight);
#endif	

        // Finished with this node, move to find next wheel weight
	*cur_node = *cur_node + 1;
		         
     ++*l;                      /* l is now first child */
     switch ( parent->type )
     {  case FUNC_DATA:
        case EVAL_DATA:
          for ( i = 0; i < parent->arity; ++i )
          {    
               generate_crossover_wheel_recurse_lcl1(l,crossover_wheel,cur_node,trueIndex,fixed);
		/* going back up */  
          }
          break;
        case FUNC_EXPR:
        case EVAL_EXPR:
          for ( i = 0; i < parent->arity; ++i )
          {    ++*l;           /* l was skipnode, now its the child */
              
               generate_crossover_wheel_recurse_lcl1(l,crossover_wheel,cur_node,trueIndex,fixed);
		/* going back up*/
          }
          break;
     }
}	



/* operator_crossver_init()
 *
 * called to parse crossover options and initialize one record
 * of a breedphase table appropriately.
 */

int operator_crossover2Local_init ( char *options, breedphase *bp )
{
     int errors = 0;
     crossover_data *cd;
     int i, j, k, m;
     double r;
     char **argv, **targv;
     int internalset = 0, externalset = 0;
     char *cp;


    // mwh : set to not use level 2 heuristics unless requested
    char *p;
    if (p=get_parameter( "acgp.level")) {
         printf("We are using level 2 heuristics operator crossover2Local.\n");
         Acgp_level = atof(p);
    }
    else {
         Acgp_level = 1;
    }

     if (Acgp_level < 2) {
     	printf("Trying to use level 2 heuristic operator(crossover2Local) when acgp_level is not greater than 1!\n");
	printf("Please change to a different operator or set acgp_level=2 in parameter file\n");
	exit(1);
     }

     cd = (crossover_data *)MALLOC ( sizeof ( crossover_data ) );

     /* place values into the breedphase table record. */
     bp->operator = OPERATOR_CROSSOVER2LOCAL;
     bp->data = (void *)cd;
     bp->operator_free = operator_crossover2Local_free;
     bp->operator_start = operator_crossover2Local_start;
     bp->operator_end = operator_crossover2Local_end;
     bp->operator_operate = operator_crossover2Local;

     /* default values for all the crossover options. */
     cd->keep_trying = 0;
     cd->internal = 0.9;
     cd->external = 0.1;
     cd->tree = (double *)MALLOC ( tree_count * sizeof ( double ) );
     cd->treecumul = (double *)MALLOC ( tree_count * sizeof ( double ) );
     for ( j = 0; j < tree_count; ++j )
          cd->tree[j] = 0.0;
     cd->treetotal = 0.0;
     cd->func = (double *)MALLOC ( fset_count * sizeof ( double ) );
     cd->sname = NULL;
     cd->sname2 = NULL;

     /* break the options string into an argv-style array of strings. */
     j = parse_o_rama ( options, &argv );

     for ( i = 0; i < j; ++i )
     {
	  /* parse "keep_trying" option. */
          if ( strcmp ( "keep_trying", argv[i] ) == 0 )
          {
	       /* translate a string into a binary value.  returns -1 if
		  the string is not one of the valid strings meaning
		  yes or no. */
               cd->keep_trying = translate_binary ( argv[++i] );
               if ( cd->keep_trying == -1 )
               {
                    ++errors;
                    error ( E_ERROR, "crossover: \"%s\" is not a valid setting for \"keep_trying\".",
                           argv[i] );
               }
          }
	  /* parse "internal" option. */
          else if ( strcmp ( "internal", argv[i] ) == 0 )
          {
               internalset = 1;
               cd->internal = strtod ( argv[++i], NULL );
               if ( cd->internal < 0.0 )
               {
                    ++errors;
                    error ( E_ERROR, "crossover: \"internal\" must be nonnegative." );
               }
          }
	  /* parse "external" option. */
          else if ( strcmp ( "external", argv[i] ) == 0 )
          {
               externalset = 1;
               cd->external = strtod ( argv[++i], NULL );
               if ( cd->external < 0.0 )
               {
                    ++errors;
                    error ( E_ERROR, "crossover: \"external\" must be nonnegative." );
               }
          }
	  /* parse "select" option. */
          else if ( strcmp ( "select", argv[i] ) == 0 )
          {
               if ( !exists_select_method ( argv[++i] ) )
               {
                    ++errors;
                    error ( E_ERROR, "crossover: \"%s\" is not a known selection method.",
                           argv[i] );
               }
               FREE ( cd->sname );
               cd->sname = (char *)MALLOC ( (strlen(argv[i])+1) * sizeof ( char ) );
               strcpy ( cd->sname, argv[i] );
               if ( cd->sname2 == NULL )
                    cd->sname2 = cd->sname;
          }
	  /* parse "select2" option. */
          else if ( strcmp ( "select2", argv[i] ) == 0 )
          {
               if ( !exists_select_method ( argv[++i] ) )
               {
                    ++errors;
                    error ( E_ERROR, "crossover: \"%s\" is not a known selection method.",
                           argv[i] );
               }
               if ( cd->sname2 && cd->sname != cd->sname2 )
                    FREE ( cd->sname2 );
               cd->sname2 = (char *)MALLOC ( (strlen(argv[i])+1) * sizeof ( char ) );
               strcpy ( cd->sname2, argv[i] );
          }
	  /* parse "tree" option. */
          else if ( strcmp ( "tree", argv[i] ) == 0 )
          {
               k = parse_o_rama ( argv[++i], &targv );
               if ( k != tree_count )
               {
                    ++errors;
                    error ( E_ERROR, "crossover: wrong number of tree fields: \"%s\".",
                           argv[i] );
               }
               else
               {
                    for ( m = 0; m < k; ++m )
                    {
                         cd->tree[m] = strtod ( targv[m], &cp );
                         if ( *cp )
                         {
                              ++errors;
                              error ( E_ERROR, "crossover: \"%s\" is not a number.",
                                     targv[m] );
                         }
                    }
               }
               
               free_o_rama ( k, &targv );
          }
	  /* parse "tree#" option. */
          else if ( strncmp ( "tree", argv[i], 4 ) == 0 )
          {
               k = strtol ( argv[i]+4, &cp, 10 );
               if ( *cp )
               {
                    ++errors;
                    error ( E_ERROR, "crossover: unknown option \"%s\".",
                           argv[i] );
               }
               if ( k < 0 || k >= tree_count )
               {
                    ++errors;
                    error ( E_ERROR, "crossover: \"%s\" is out of range.",
                           argv[i] );
               }
               else
               {
                    cd->tree[k] = strtod ( argv[++i], &cp );
                    if ( *cp )
                    {
                         ++errors;
                         error ( E_ERROR, "crossover: \"%s\" is not a number.",
                                argv[i] );
                    }
               }
          }
          else
          {
               ++errors;
               error ( E_ERROR, "crossover: unknown option \"%s\".",
                      argv[i] );
          }
     }
     
     free_o_rama ( j, &argv );
     
     if ( internalset && !externalset )
          cd->external = 0.0;
     else if ( !internalset && externalset )
          cd->internal = 0.0;
     
     if ( cd->sname == NULL )
     {
          ++errors;
          error ( E_ERROR, "crossover: no selection method specified." );
     }

     /** compute "func" array from the "tree" array. **/
     
     for ( j = 0; j < tree_count; ++j )
          cd->treetotal += cd->tree[j];
     if ( cd->treetotal == 0.0 )
     {
          for ( j = 0; j < tree_count; ++j )
               cd->tree[j] = 1.0;
          cd->treetotal = tree_count;
     }
          
     for ( j = 0; j < fset_count; ++j )
          cd->func[j] = 0.0;
     for ( j = 0; j < tree_count; ++j )
          cd->func[tree_map[j].fset] += cd->tree[j];
     
     r = 0.0;
     for ( j = 0; j < fset_count; ++j )
          r = (cd->func[j] += r);

#ifdef DEBUG
     if ( !errors )
     {
          printf ( "crossover options:\n" );
          printf ( "   internal: %lf  external: %lf\n", cd->internal, cd->external );
          printf ( "   keep_trying: %d\n", cd->keep_trying );
          printf ( "   primary selection: %s\n", cd->sname==NULL?"NULL":cd->sname );
          printf ( "   second selection: %s\n", cd->sname2==NULL?"NULL":(cd->sname2==cd->sname?"same as primary":cd->sname2) );
          printf ( "   tree total: %lf\n", cd->treetotal );
          for ( j = 0; j < tree_count; ++j )
               printf ( "   tree %d: %lf\n", j, cd->tree[j] );
          for ( j = 0; j < fset_count; ++j )
               printf ( "   fset %d: %lf\n", j, cd->func[j] );
     }
#endif
     
     return errors;
}

/* operator_crossover_free()
 *
 * free the crossover-specific data structure.
 */

void operator_crossover2Local_free ( void *data )
{
     crossover_data * cd;

     cd = (crossover_data *)data;

     FREE ( cd->sname );
     if ( cd->sname != cd->sname2 )
          FREE ( cd->sname2 );
     FREE ( cd->tree );
     FREE ( cd->treecumul );
     FREE ( cd->func );
     FREE ( cd );
}

/* operator_crossover_start()
 *
 * called at the start of the breeding process each generation.
 * initializes the selection contexts for this phase.
 */

void operator_crossover2Local_start ( population *oldpop, void *data )
{
     crossover_data * cd;
     select_context_func_ptr select_con;

     cd = (crossover_data *)data;
     
     select_con = get_select_context ( cd->sname );
     cd->sc = select_con ( SELECT_INIT, NULL, oldpop, cd->sname );

     /* if there is a separate selection method specified for the
	second parent... */
     if ( cd->sname2 != cd->sname )
     {
	  /* ...then initialize it too. */
          select_con = get_select_context ( cd->sname2 );
          cd->sc2 = select_con ( SELECT_INIT, NULL, oldpop, cd->sname2 );
     }
		 else
	  /* ...otherwise use the first context. */
          cd->sc2 = cd->sc;
}

/* operator_crossover_end()
 *
 * called when breeding is finished each generation.  frees up selection
 * contexts for this phase.
 */

void operator_crossover2Local_end ( void *data )
{
     crossover_data * cd;

     cd = (crossover_data *)data;

     cd->sc->context_method ( SELECT_CLEAN, cd->sc, NULL, NULL );
     if ( cd->sname != cd->sname2 )
          cd->sc2->context_method ( SELECT_CLEAN, cd->sc2, NULL, NULL );
}

/* operator_crossover()
 *
 * performs the crossover, inserting one or both offspring into the
 * new population.
 */
/* czj: modified for dealing with constraints */

void operator_crossover2Local ( population *oldpop, population *newpop,
                         void *data )
{ crossover_data * cd;
  int p1, p2;
  int ps1, ps2;
  int l1,l2;
  lnode *st[3];
  int sts1, sts2;
  int ns1, ns2;
  int badtree=0;
  double total;
  int forceany1, forceany2;
  int repcount;
  int repSrc_czj=0;            /* counts repeats if no feasible source found */
  int repBad_czj=0;             /* counts repeats if offspring violates size */
  int f, t1, t2, j;
  double r, r2;
  int totalnodes;
  int i;
#ifdef DEBUG_CROSSOVER 
  static int numCalls_czj=0;                           /* czj for debugging */
  numCalls_czj++; /* czj debug */
#endif
   
  /* get the crossover-specific data structure. */
  cd = (crossover_data *)data;
  total = cd->internal + cd->external;

  /* choose a function set. */
  r = random_double() * cd->treetotal;

  for ( f = 0; r >= cd->func[f]; ++f );

  /* fill in the "treecumul" array, zeroing all trees which
	don't use the selected function set. */
  r = 0.0;
  t1 = 0;
  for ( j = 0; j < tree_count; ++j )
  { if ( tree_map[j].fset == f )
      r = (cd->treecumul[j] = r + cd->tree[j]);
    else
      cd->treecumul[j] = r;
  }

  /* select the first and second trees. */ 
  r2 = random_double() * r;
  for ( t1 = 0; r2 >= cd->treecumul[t1]; ++t1 );

  r2 = random_double() * r;
  for ( t2 = 0; r2 >= cd->treecumul[t2]; ++t2 );

#ifdef DEBUG_CROSSOVER 
  printf ("crossover # %d\n",numCalls_czj); /* czj debug */
  printf ( "selected function set %d --> t1: %d; t2: %d\n", f, t1, t2 );
#endif 
     
  /* choose two parents */
  p1 = cd->sc->select_method ( cd->sc );
  ps1 = oldpop->ind[p1].tr[t1].nodes;

               /* if the tree only has one node, we obviously can't do
	                     fucntionpoint crossover.  use anypoint instead. */
  forceany1 = (ps1==1||total==0.0);

  p2 = cd->sc2->select_method ( cd->sc2 );
  ps2 = oldpop->ind[p2].tr[t2].nodes;
  forceany2 = (ps2==1||total==0.0);

#ifdef DEBUG_CROSSOVER
  printf ( "Parent 1(%d) is:\n" ,p1);
  print_individual ( oldpop->ind+p1, stderr ); 
  printf ( "Parent 2(%d) is:\n",p2 );
  print_individual ( oldpop->ind+p2, stderr ); 
printf("Parent1 has %d nodes, parent2 has %d nodes\n",ps1,ps2); /* czj debug */
#endif 

  repSrc_czj=0;
  while (repSrc_czj++<RepeatsSrc_czj)
               /* repeat when no feasible source found, up to RepeatsSrc_czj */
  { 
  
    // Selects a destination, that is where we are starting from
    if (forceany1)                                            /* choose dest */
    { l1 = random_int(ps1);
      st[1] = get_subtree(oldpop->ind[p1].tr[t1].data, l1);
    }
    else if (total*random_double() < cd->internal)   /* choose internal dest */
         { l1=random_int(tree_nodes_internal(oldpop->ind[p1].tr[t1].data)); 
           st[1] = get_subtree_internal(oldpop->ind[p1].tr[t1].data, l1);
        }
         else                                     /* choose an external dest */
         { l1=random_int(tree_nodes_external(oldpop->ind[p1].tr[t1].data)); 
           st[1] = get_subtree_external(oldpop->ind[p1].tr[t1].data,l1);
         }
	
    // We need to traverse the tree to gather context information
    // Either do myself or traverse internal to get_subtree and then
    // go the rest of the way. Also we need to do this to find out the depth
    // of our counter, so we know if we need to do this
      int replace_depth = 0;
	  
      int trueIndex;
      int newNode;
      int thisArity;
      // gather context info, putting it in global array level1context_mwh
      // Only deal with it if it has more than just a root
      if (ps1 > 1) {
      	   gather_info_mwh(oldpop->ind[p1].tr[t1].data,st[1],&replace_depth,&thisArity,&trueIndex,level1context_mwh);
       }
       else {
       		// Just has a root
		replace_depth = 0;
       } 
	 	  
#ifdef DEBUG_CROSSOVER
       fprintf(stderr,"subtree number: %d\n",l1/*l*/);         /* added by mwh */
       fprintf ( stderr, "selected for swap by crossover: " );
       print_tree ( /*replace[0]*/ st[1], stderr );
       printf("Context information for crossover\n");
       printf("DEPTH SOURCE=%d\n",replace_depth);
       if (replace_depth > 1) {
       	printf("Cntxt info is (%s)(",fset[0].cset[trueIndex].string);
      	 for (i = 0; i < fset[0].cset[trueIndex].arity; i++) {
  		    printf("%s ",fset[0].cset[level1context_mwh[i]].string);
      	 }
       printf(")\n");
       }
#endif
    
    // DO NOT NEED THIS CODE FOR ONLY LCL
       // Now generate our weight array to be used when we traverse source
    //   if (replace_depth == 1) 
    //   {
    //   	    // Sets crossover_weights_mwh[], so element 3 of this array would be
	    // fset's function 3's weight given the context information
//       	    generate_crossover_weights_lvl1_lcl1_for_lvl1_mwh(trueIndex,level1context_mwh,thisArity);
 //      } else
  
       if (replace_depth > 1) 
      {
            // Sets crossover_weights_mwh[], so element 3 of this array would be
	    // fset's function 3's weight given the context information
       	    generate_crossover_weights_lcl1_for_lcl1_mwh(trueIndex,level1context_mwh,thisArity);
       }
       
            
   	 //TODO: REMOVE OLD CODE BELOW HERE   
 	 /* czj: note that Function_/Argumnet_czj are set by calls to get_subtree* */
    repBad_czj=0;                /* reset the limit counter for the new dest */
 	if (markXNodes_czj(oldpop->ind[p2].tr[t2].data)==0)
 	    continue;                        /* no sources found, try another dest */
    while (repBad_czj++<RepeatsBad_czj)  /* limit on alt. srcs for this dest */
    { 
    
    
    
    if (replace_depth > 1) {
       // Now select our source from our own wheel, so spin it again
       // Might want to also, if it ends up being bad, to set the destination probability to
       // zero for the next time
       //TODO: Perhaps only let it take all nodes or just internal nodes
       int num_in_source = spin_after_generating_crossover_weights_lcl1_mwh(oldpop->ind[p2].tr[t2].data,ps2,trueIndex,thisArity);
       
       // now need to traverse this tree to our num_in_source
       st[2] = get_subtree(oldpop->ind[p2].tr[t2].data, num_in_source);
       
       
    }
    else {
      // Lets do normal czj picking of source
      if (forceany2)                 /* here only if some sources were found */
       st[2]=getSubtreeMarked_czj(oldpop->ind[p2].tr[t2].data,0);
     else if (total*random_double() < cd->internal)
           st[2]=getSubtreeMarked_czj(oldpop->ind[p2].tr[t2].data,1);
          else 
            st[2]=getSubtreeMarked_czj(oldpop->ind[p2].tr[t2].data,2);
	    
      }	     
	     
#ifdef DEBUG_CROSSOVER 
      printf("p1 (dest1) selected %dth node\n",l1);    
      printf ( "Dest subtree is: " );
      print_tree ( st[1], stdout );
      printf ( "Source subtree is: " );
      print_tree ( st[2], stdout );
     // if (replace_depth == 1) {
      //	printf("Exiting out until we figure out problem with lvl 1 duplicate trees\n");
       //  exit(0);
	//}
#endif 
      badtree = 0;
      sts1 = tree_nodes(st[1]);      /* count nodes in the selected subtrees */
      sts2 = tree_nodes(st[2]);
      ns1 = ps1 - sts1 + sts2;       /* calculate the sizes of the offspring */
      totalnodes = ns1;
#ifdef DEBUG_CROSSOVER
      printf("Dest subtr has size %d, source subtree has size %d\n",sts1,sts2);
      printf("Newtree 1 has size %d; limit is %d\n",ns1,tree_map[t1].nodelimit);
#endif
      /* validate the first offspring against the tree node and depth limits */ 
      if ( tree_map[t1].nodelimit > -1 && ns1 > tree_map[t1].nodelimit )
        badtree = 1;
      else if ( tree_map[t1].depthlimit > -1 )
           { ns1=tree_depth_to_subtree(oldpop->ind[p1].tr[t1].data,st[1])
                 + tree_depth(st[2]);
#ifdef DEBUG_CROSSOVER
             printf("Newtree 1 has depth %d; limit is %d\n",ns1,
                    tree_map[t1].depthlimit);
#endif
             if ( ns1 > tree_map[t1].depthlimit )
               badtree = 1;
           }
                                /* check node limits on the whole individual */
      if ( ind_nodelimit > -1 )
      { for ( i = 0; i < tree_count; ++i )
          if ( i != t1 )
            totalnodes += oldpop->ind[p1].tr[i].nodes;
        badtree |= (totalnodes > ind_nodelimit);
#ifdef DEBUG_CROSSOVER
        printf("Newind 1 has %d nodes; limit is %d\n",totalnodes,ind_nodelimit);
#endif
      }
      if (!(cd->keep_trying && badtree)) 
        break;    /* from inner loop: either not bad offspring or over limit */
    }
    duplicate_individual(newpop->ind+newpop->next,oldpop->ind+p1);
    if (!badtree)   /* if tree not bad replace the parent with the offspring */
    { 
#ifdef DEBUG_CROSSOVER
      fprintf ( stderr, "Offspring 1 is allowable.\n" );
#endif
                         /* make a copy of the crossover tree, replacing the */
         	          /* selected subtree with the crossed-over subtree. */
      copy_tree_replace_many(0,oldpop->ind[p1].tr[t1].data,st+1,st+2,1,
                             &repcount);
      if (repcount != 1)             /* this can't happen, but check anyway. */
        error ( E_FATAL_ERROR,"botched crossover:  this can't happen" );
                          /* free the appropriate tree of the new individual */
      free_tree ( newpop->ind[newpop->next].tr+t1 );
                             /* copy the crossovered tree to the freed space */
      gensp_dup_tree (0,newpop->ind[newpop->next].tr+t1);
               /* the new individual's fitness fields are of course invalid. */
      newpop->ind[newpop->next].evald = EVAL_CACHE_INVALID;
      newpop->ind[newpop->next].flags = FLAG_NONE;
    }
#ifdef DEBUG_CROSSOVER
    else
      fprintf ( stderr, "Offspring 1 not generated; copying parent 1.\n" );
#endif
    ++newpop->next;               /* an offspring created (or parent copied) */
#ifdef DEBUG_CROSSOVER
    fprintf ( stderr, "Offspring 1:" );
    if ( newpop->ind[newpop->next-1].evald == EVAL_CACHE_VALID )
      fprintf ( stderr, "  (valid)\n" );
    else
      fprintf ( stderr, "  (invalid)\n" );
    print_individual ( newpop->ind+(newpop->next-1), stderr );
#endif          
#if VERIFY_CROSSOVER_czj /* czj */
    if (verify_tree_czj ( newpop->ind[newpop->next-1].tr[t1].data))
    { oprintf(OUT_SYS, 10, "tree 1");
      /* print_tree ( newpop->ind[newpop->next-1].tr[t1].data, stderr ); */
      oprintf (OUT_SYS, 10, "INVALID TREE in crossover: \n");
      fprintf(stderr,"INVALID TREE in crossover \n");
      exit(1);
    }
#endif
    break;     /* will break from the outer since offspring has been created */
  }                                                 /* end of the outer loop */

             /* if the new population needs another member (it's not full) */
             /* czj: this is new code - do all over with p2 as destination */
  if ( newpop->next < newpop->size )   
  { 
        
       
    repSrc_czj=0;
    while (repSrc_czj++<RepeatsSrc_czj)
    { if (forceany2)                  
      { l2 = random_int(ps2);
        st[2] = get_subtree(oldpop->ind[p2].tr[t2].data,l2);
      }
      else if (total*random_double() < cd->internal)  
           { l2=random_int(tree_nodes_internal(oldpop->ind[p2].tr[t2].data)); 
             st[2]=get_subtree_internal( oldpop->ind[p2].tr[t2].data,l2);
           }
           else                                   
           { l2=random_int(tree_nodes_external(oldpop->ind[p2].tr[t2].data)); 
             st[2]=get_subtree_external(oldpop->ind[p2].tr[t2].data,l2);
           }                        /* Function_czj/Argument_czj will be set */
	   
	   int replace_depth = 0;
	   int trueIndex;
   	   int newNode;
     	 int thisArity;
      	// gather context info, putting it in global array level1context_mwh
	// Only gather info if it has more than just a root, that is has depth 0
	 if (ps2 > 1) {
     	 	gather_info_mwh(oldpop->ind[p2].tr[t2].data,st[2],&replace_depth,&thisArity,&trueIndex,level1context_mwh); 
		
		// DO NOT NEED THIS CODE IN LCL
		// If necessary, generate our weights
		//if (replace_depth == 1) 
		//{
		//	// Sets crossover_weights_mwh[], so element 3 of this array would be
	    	//	// fset's function 3's weight given the context information
       	    	//	generate_crossover_weights_lvl1_lcl1_for_lvl1_mwh(trueIndex,level1context_mwh,thisArity);
		//}
		//else 
		
		if (replace_depth > 1) 
       		{
            		// Sets crossover_weights_mwh[], so element 3 of this array would be
			// fset's function 3's weight given the context information
       	    		generate_crossover_weights_lvl1_lcl1_for_lcl1_mwh(trueIndex,level1context_mwh,thisArity);
       		}  
	 }
	 else {
	 	replace_depth = 0;
	 }
	  
	  	   
      repBad_czj=0;
      if (markXNodes_czj(oldpop->ind[p1].tr[t1].data)==0)
        continue;                            /* no sources, try another dest */
      while (repBad_czj++<RepeatsBad_czj)
      { 
        if (replace_depth > 1) {
		// is at lvl 1, so spin our wheel
	 // Now select our dest from our own wheel, so spin it again
      	 // Might want to also, if it ends up being bad, to set the destination probability to
      	 // zero for the next time
     	  //TODO: Perhaps only let it take all nodes or just internal nodes
     	  int num_in_source = spin_after_generating_crossover_weights_lcl1_mwh(oldpop->ind[p1].tr[t1].data,ps1,trueIndex,thisArity);
       
     	  // now need to traverse this tree to our num_in_source
     	  st[1] = get_subtree(oldpop->ind[p1].tr[t1].data, num_in_source);
	
	}
	else {
           // Use old method
     	   if (forceany1)
    	      st[1]=getSubtreeMarked_czj(oldpop->ind[p1].tr[t1].data,0);
    	    else if (total*random_double() < cd->internal)
     	          st[1]=getSubtreeMarked_czj(oldpop->ind[p1].tr[t1].data,1);
          	   else
        	       st[1]=getSubtreeMarked_czj(oldpop->ind[p1].tr[t1].data,2);
	}       
	       
	       
#ifdef DEBUG_CROSSOVER                                       /* czj modified */
        printf("p2 (dest) selected %dth node\n",l2); /* czj */
        printf ( "Dest subtree is: " );
        print_tree ( st[2], stdout );
        printf("p1 (source) selected %dth node\n",l1);
        printf("Source subtree is: ");
        print_tree(st[1],stdout);
	// mwh context information debugging
	 printf("Context information for crossover\n");
     	  printf("DEPTH SOURCE=%d\n",replace_depth);
      	 if (replace_depth > 1) {
      	 printf("Cntxt info is (%s)(",fset[0].cset[trueIndex].string);
      	 for (i = 0; i < fset[0].cset[trueIndex].arity; i++) {
  		    printf("%s ",fset[0].cset[level1context_mwh[i]].string);
       	}
       	printf(")\n");
       }
#endif 
        badtree = 0;
        sts2=tree_nodes(st[2]); 
        sts1=tree_nodes(st[1]);
        ns2=ps2-sts2+sts1; 
        ns1=ps1-sts1+sts2;   
        totalnodes = ns2;
#ifdef DEBUG_CROSSOVER
        printf("Dest subtree has size %d, source subtree2 has size %d\n",
               sts2,sts1);
        printf("Newtree 2 has size %d; limit is %d\n",ns1,
               tree_map[t2].nodelimit);
#endif
      /* validate the first offspring against the tree node and depth limits */ 
        if ( tree_map[t2].nodelimit > -1 && ns2 > tree_map[t2].nodelimit )
          badtree = 1;
        else if (tree_map[t2].depthlimit > -1)
             { ns2=tree_depth_to_subtree(oldpop->ind[p2].tr[t2].data,st[2])
                   + tree_depth(st[1]);
#ifdef DEBUG_CROSSOVER
               printf("Newtree 2 has depth %d; limit is %d\n",ns2,
                      tree_map[t2].depthlimit);
#endif
               if (ns2 > tree_map[t2].depthlimit)
                 badtree = 1;
             }
                                /* check node limits on the whole individual */
        if (ind_nodelimit > -1)
        { for (i = 0; i < tree_count; ++i)
            if (i != t1)
              totalnodes += oldpop->ind[p2].tr[i].nodes;
          badtree |= (totalnodes > ind_nodelimit);
#ifdef DEBUG_CROSSOVER
          printf("Newind 2 has %d nodes; limit is %d\n",totalnodes,
                 ind_nodelimit);
#endif
        }
        if (!(cd->keep_trying && badtree)) 
          break;         
      }
      duplicate_individual(newpop->ind+newpop->next,oldpop->ind+p2);
      if (!badtree)
      { 
#ifdef DEBUG_CROSSOVER
        fprintf(stderr,"Offspring 2 is allowable.\n");
#endif
        copy_tree_replace_many(0,oldpop->ind[p2].tr[t2].data,
                               st+2, st+1, 1, &repcount );
        if (repcount != 1)
          error(E_FATAL_ERROR,"botched crossover:  this can't happen");
        free_tree(newpop->ind[newpop->next].tr+t2);
        gensp_dup_tree(0, newpop->ind[newpop->next].tr+t2);
        newpop->ind[newpop->next].evald = EVAL_CACHE_INVALID;
        newpop->ind[newpop->next].flags = FLAG_NONE;
      }
#ifdef DEBUG_CROSSOVER
      else
        fprintf(stderr,"Offspring 2 not generated; copying parent 2\n");
#endif
      ++newpop->next;
        
#ifdef DEBUG_CROSSOVER
      fprintf ( stderr, "Offspring 2:" );
      if ( newpop->ind[newpop->next-1].evald == EVAL_CACHE_VALID )
        fprintf ( stderr, "  (valid)\n" );
      else
        fprintf ( stderr, "  (invalid)\n" );
      print_individual ( newpop->ind+(newpop->next-1), stderr );
#endif          
#if VERIFY_CROSSOVER_czj 
      if (verify_tree_czj ( newpop->ind[newpop->next-1].tr[t2].data))
      { oprintf(OUT_SYS, 10, "tree 2");
        /* print_tree ( newpop->ind[newpop->next-1].tr[t2].data, stdout ); */
        oprintf (OUT_SYS, 10, "INVALID TREE in crossover\n");
        fprintf(stderr,"INVALID TREE in crossover\n");
        exit(1);
      }
#endif
      break;
    }                                              /* end of the outer while */
  }                                                             /* end of if */ 
#ifdef DEBUG_CROSSOVER
  printf ( "CROSSOVER COMPLETE.\n\n\n" );
#endif
}
