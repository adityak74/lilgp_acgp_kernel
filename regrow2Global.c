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
extern int *level1context_mwh;	// list of the indices of lvl 1 context data
extern int *level2context_mwh;  // list of indices of lvl 2 context data
extern int size_context1_mwh;   // number of elements in level1context_mwh
extern int size_context2_mwh;   // number of elements in level2context_mwh
extern int acgp_level;


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
} regrow_data;

/* operator_regrow_init()
 *
 * initialize a breedphase record for mutation.  much of this is cut and
 * paste from the crossover operator; see that for (a few) more comments.
 */

int operator_regrow2Global_init ( char *options, breedphase *bp )
{
     int errors = 0;
     regrow_data *md;
     int i, j, k, m;
     double r;
     char **argv, **targv;
     int internalset = 0, externalset = 0;
     char *cp;


    // mwh : set to not use level 2 heuristics unless requested
    char *p;
    if (p=get_parameter( "acgp.level")) {
         printf("We are using level 2 heuristics operator regrow2Global.\n");
         Acgp_level = atof(p);
    }
    else {
         Acgp_level = 1;
    }

     if (Acgp_level < 2) {
     	printf("Trying to use level 2 heuristic operator(regrow2Global) when acgp_level is not greater than 1!\n");
	printf("Please change to a different operator or set acgp_level=2 in parameter file\n");
	exit(1);
     }

     md = (regrow_data *)MALLOC ( sizeof ( regrow_data ) );

     /* fill in the breedphase record. */
     bp->operator = OPERATOR_REGROW2GLOBAL;
     bp->data = (void *)md;
     bp->operator_free = operator_regrow2Global_free;
     bp->operator_start = operator_regrow2Global_start;
     bp->operator_end = operator_regrow2Global_end;
     bp->operator_operate = operator_regrow2Global;

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

/* operator_regrow_free()
 *
 * free mutation stuff.
 */

void operator_regrow2Global_free ( void *data )
{
     regrow_data * md;

     md = (regrow_data *)data;
     FREE ( md->sname );
     FREE ( md->tree );
     FREE ( md );
}

/* operator_regrow_start()
 *
 * get selection context for mutation operator.
 */

void operator_regrow2Global_start ( population *oldpop, void *data )
{
     regrow_data * md;
     select_context_func_ptr select_con;

     md = (regrow_data *)data;
     select_con = get_select_context ( md->sname );
     md->sc = select_con ( SELECT_INIT, NULL, oldpop, md->sname );
}

/* operator_regrow_end()
 *
 * free selection context for mutation operator.
 */

void operator_regrow2Global_end ( void *data )
{
     regrow_data * md;

     md = (regrow_data *)data;
     md->sc->context_method ( SELECT_CLEAN, md->sc, NULL, NULL );
}

/* operator_regrow()
 *
 * do the mutation.
 */


void generate_random_full_tree_wrapper_lvl1_mwh(int space, int depth, function_set *fset,const int source) 
/*
 Depth is the depth left until the bottom
*/
{
		
	
	int i,j;
	int save;
	
	if (source == REGROW_SOURCE) {
	
		if (depth <= 0)
			i=random_T_czj();
		else 
			i=random_F_czj();
		
		gensp_next(space)->f = (fset->cset)+i;          /* store whatever it is */
     		
		if (i>=fset->function_count)                   /* if this is a terminal */
     		{ if ( (fset->cset)[i].ephem_gen )      /* if ephem term generate value */
         		gensp_next(space)->d = new_ephemeral_const ( (fset->cset)+i );
       			return;
     		}
		
				
		// Spin the wheel, this gives us the children for 'i'
		if (depth  <= 1)
			// force all terminals
			spin_lvl1_information_terminals_mwh(i);
		else if (depth - 1 > MinDepth_czj)
			spin_lvl1_information_functions_mwh(i);
		else
			spin_lvl1_information_mwh(i);
		
		// Now we have our child_list[c], but it is index to our cube, we need to convert it to indexes to fset
		int c;
		for (c =0; c < (fset->cset)[i].arity; c++) {
			child_list[c] = MS_czj[i][c].members[child_list[c]];
		}
		
		/* printf("Wheel for %d spun, external actual fset winner was [",i);
		for (c = 0; c < (fset->cset)[i].arity; c++)
			printf("%d ",child_list[c]);
		printf("]\n"); */
				
		switch( (fset->cset)[i].type)
		{ case FUNC_DATA: case EVAL_DATA:
			for (j=0; j < (fset->cset)[i].arity; ++j)
			{ Function_czj = i;
			  Argument_czj = j;
			  generate_random_full_tree_wrapper_lvl1_mwh(space,depth-1,fset,WRAPPER_SOURCE);
			}
			break;
		   case FUNC_EXPR: case EVAL_EXPR:
		   	for (j = 0; j < (fset->cset)[i].arity; ++j)
			{	Function_czj = i;
				Argument_czj = j;
				save = gensp_next_int(space);
				generate_random_full_tree_wrapper_lvl1_mwh(space,depth-1,fset,WRAPPER_SOURCE);
				gensp[space].data[save].s = gensp[space].used-save-1;
			}
			break;
		
		
		}
		return;
	}
	
	
	if ( source == WRAPPER_SOURCE) {
		// Get our child from child_list[Argument_czj]
		i = child_list[Argument_czj];
		
		// Lets generate our space for the root
		gensp_next(space)->f = (fset->cset)+i;
		
		if (i>= fset->function_count) {
			if ((fset->cset)[i].ephem_gen)
				gensp_next(space)->d = new_ephemeral_const((fset->cset)+i);
				return;
		}
		
		switch( (fset->cset)[i].type)
		{ case FUNC_DATA: case EVAL_DATA:
			for (j=0; j < (fset->cset)[i].arity; ++j)
			{ Function_czj = i;
			  Argument_czj = j;
			  generate_random_full_tree(space,depth-1,fset);
			}
			break;
		   case FUNC_EXPR: case EVAL_EXPR:
		   	for (j = 0; j < (fset->cset)[i].arity; ++j)
			{	Function_czj = i;
				Argument_czj = j;
				save = gensp_next_int(space);
				generate_random_full_tree(space,depth-1,fset);
				gensp[space].data[save].s = gensp[space].used-save-1;
			}
			break;
		
		
		}
		return;
	}
	
	
	printf("Should never get here!\n");
	printf("Being here means that this function was called below lvl 1 and 2 by regrow operator\n");
	exit(1);
	

}


void generate_random_grow_tree_wrapper_lvl1_mwh(int space, int depth, function_set *fset,const int source) 
/*
 Depth is the depth left until the bottom
*/
{
		
	int i,j;
	int save;
	
	if (source == REGROW_SOURCE) {
	
		if (depth <= 0)
			i=random_T_czj();
		else if (depth>MinDepth_czj)
			i=random_F_czj();
		else
			i=random_FT_czj();
		
		// Lets generate our space for the root
		gensp_next(space)->f = (fset->cset)+i;
		if (i>=fset->function_count)                   /* if this is a terminal */
     		{ if ( (fset->cset)[i].ephem_gen )      /* if ephem term generate value */
         		gensp_next(space)->d = new_ephemeral_const ( (fset->cset)+i );
       			return;
     		}	
		
		// Spin the wheel, this gives us the children for 'i'
		if (depth  <= 1)
			// force all terminals
			spin_lvl1_information_terminals_mwh(i);
		else if (depth -1 > MinDepth_czj)
			// force at least one function
			spin_lvl1_information_functions_mwh(i);
		else
			// allow anything
			spin_lvl1_information_mwh(i);
			
		// Now we have our child_list[c], but it is index to our cube, we need to convert it to indexes to fset
		int c;
		for (c =0; c < (fset->cset)[i].arity; c++) {
			child_list[c] = MS_czj[i][c].members[child_list[c]];
		}
		
				
		/* printf("Wheel for %d spun, external actual fset winner was [",i);
		for (c = 0; c < (fset->cset)[i].arity; c++)
			printf("%d ",child_list[c]);
		printf("]\n"); */
		
		
		
		switch( (fset->cset)[i].type)
		{ case FUNC_DATA: case EVAL_DATA:
			for (j=0; j < (fset->cset)[i].arity; ++j)
			{ Function_czj = i;
			  Argument_czj = j;
			  generate_random_grow_tree_wrapper_lvl1_mwh(space,depth-1,fset,WRAPPER_SOURCE);
			}
			break;
		   case FUNC_EXPR: case EVAL_EXPR:
		   	for (j = 0; j < (fset->cset)[i].arity; ++j)
			{	Function_czj = i;
				Argument_czj = j;
				save = gensp_next_int(space);
				generate_random_grow_tree_wrapper_lvl1_mwh(space,depth-1,fset,WRAPPER_SOURCE);
				gensp[space].data[save].s = gensp[space].used-save-1;
			}
			break;
		
		
		}
		return;
	}
	
	
	if ( source == WRAPPER_SOURCE) {
		// Get our child from child_list[Argument_czj]
		i = child_list[Argument_czj];
		
		// Lets generate our space for the root
		gensp_next(space)->f = (fset->cset)+i;
		
		if (i>= fset->function_count) {
			if ((fset->cset)[i].ephem_gen)
				gensp_next(space)->d = new_ephemeral_const((fset->cset)+i);
				return;
		}
		
		switch( (fset->cset)[i].type)
		{ case FUNC_DATA: case EVAL_DATA:
			for (j=0; j < (fset->cset)[i].arity; ++j)
			{ Function_czj = i;
			  Argument_czj = j;
			  generate_random_grow_tree(space,depth-1,fset);
			}
			break;
		   case FUNC_EXPR: case EVAL_EXPR:
		   	for (j = 0; j < (fset->cset)[i].arity; ++j)
			{	Function_czj = i;
				Argument_czj = j;
				save = gensp_next_int(space);
				generate_random_grow_tree(space,depth-1,fset);
				gensp[space].data[save].s = gensp[space].used-save-1;
			}
			break;
		
		
		}
		return;
	}
	
	
	printf("Should never get here!\n");
	printf("Being here means that this function was called below lvl 1 and 2 by regrow operator\n");
	exit(1);
	

}

void operator_regrow2Global ( population *oldpop, population *newpop,
                      void *data )
{
     int i;
     int ps;
     lnode *replace[2];
     int l, ns;
     int badtree;
     int repcount;
     regrow_data * md;
     int t;
     double r;
     int depth;
     int totalnodes;
     int p;
     int forceany;
     double total;
     
     md = (regrow_data *)data;
     total = md->internal + md->external;

     /* choose a tree to regrow. */
     /* czj Note: will need to regrow all trees for the future */
     r = random_double() * md->treetotal;
     for ( t = 0; r >= md->tree[t]; ++t );

     /* select an individual to regrow. */
     p = md->sc->select_method ( md->sc ); 
     ps = tree_nodes ( oldpop->ind[p].tr[t].data );
     forceany = (ps==1||total==0.0);

#ifdef DEBUG_REGROW
     fprintf ( stderr, "the parent size is %d\n", ps );
     fprintf ( stderr, "    parent %4d: ", p );
     print_tree ( oldpop->ind[p].tr[t].data, stderr );
#endif          

     while(1)
     {

          /* czj change: just take the first element (0) from tr[t].data */
          l=0;
          replace[0] = get_subtree ( oldpop->ind[p].tr[t].data, l );
	  
#ifdef DEBUG_REGROW
          fprintf(stderr,"subtree number: %d\n",l);         /* added by czj */
          fprintf ( stderr, "selected for replacement: " );
          print_tree ( replace[0], stderr );
#endif

          gensp_reset ( 1 );
	  /* pick a value from the depth ramp. */
          depth = md->mindepth + random_int ( md->maxdepth - md->mindepth + 1 );
	  	if (md->absDepth_czj)
			MinDepth_czj=depth-md->mindepth;
	  /* grow the tree. */
          switch ( md->method )
          {
             case GENERATE_GROW:
               generate_random_grow_tree_wrapper_lvl1_mwh ( 1, depth, fset+tree_map[t].fset,REGROW_SOURCE );
               break;
             case GENERATE_FULL:
               generate_random_full_tree_wrapper_lvl1_mwh ( 1, depth, fset+tree_map[t].fset,REGROW_SOURCE );
               break;
             case GENERATE_HALF_AND_HALF:
               if ( random_double() < 0.5 )
                    generate_random_grow_tree_wrapper_lvl1_mwh ( 1, depth, fset+tree_map[t].fset,REGROW_SOURCE );  // calls mwh version for lvl 1 information
               else
                    generate_random_full_tree_wrapper_lvl1_mwh ( 1, depth, fset+tree_map[t].fset,REGROW_SOURCE );  // calls wrapper for lvl 1 information
               break;
          }

	  MinDepth_czj = 0;

#ifdef DEBUG_REGROW
          fprintf ( stderr, "the new subtree is: " );
          print_tree ( gensp[1].data, stderr );
#endif

	  /* count the nodes in the new tree. */
          ns = ps - tree_nodes ( replace[0] ) + tree_nodes ( gensp[1].data );
          totalnodes = ns;

	  /* check the regrowd tree against node count and/or size limits. */
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

	  /* check regrowd tree against whole-individual node limits. */
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
#ifdef DEBUG_REGROW
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
#ifdef DEBUG_REGROW
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
               
#ifdef DEBUG_REGROW
               fprintf ( stderr, "    the regrowd individual is: " );
               print_individual ( newpop->ind+newpop->next, stderr );
#endif
          }

          ++newpop->next;
          break;
     }

#ifdef DEBUG_REGROW
     printf ( "MUTATION COMPLETE.\n\n\n" );
#endif
}
