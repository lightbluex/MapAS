#include<time.h>
#include<stdio.h>
#include<stdlib.h>
#include <math.h>
#include <string.h>
#include <limits.h>
#include <assert.h>

#include"AS.h"
#include"IO.h"
#include"utilities.h"
#include"MapCross.h"

clock_t start_time;
clock_t start_time_i;
double elapsed;
double elapsed_t;


double *best_in_try;		/* best length of each try */
long int *best_found_at;	/* iteration when best length is found */
double   *time_best_found;	/* time when best length is found  */
double   *time_total_run;	/* total time of a try */
long int iter_to_best;
long int restart_found_best;/* iteration in which restart-best solution is found */

ant_struct *ant;
ant_struct *best_so_far_ant;



long int jump_out;
long int jump_flag;

double   *prob_of_selection;/*the probability point for the roullete*/

long int n_ants;			/* number of ants */

long int n_try;

long int ras_ranks; 

long int acs_flag;
long int ras_flag;
long int eas_flag;	//
long int mmas_flag;	//

double c_p;

long int elitist_ants;	//

double   trail_max;       /* maximum pheromone trail in MMAS */
double   trail_min;       /* minimum pheromone trail in MMAS */
long int u_gb;            /* every u_gb iterations update with best-so-far ant */

long int iteration;         /* iteration counter */
long int max_tries;         /* maximum number of independent tries */
long int max_iteration;
double   max_time;          /* maximal allowed run time of a try  */
long int optimal;           /* optimal solution or bound to find */

double alpha;         /* importance of trail */
double beta;          /* importance of heuristic evaluate */
double rho;           /* parameter for evaporation */
double Q;			  /* pheromone deposition weight */

double   trail_0;	/*pheromone initialization */

//double **pheromone;
double **heuristic;
double **total;


char name_buf[LINE_BUF_LEN]; /* directory (include file name) of data file */
struct problem instance;

int main(int argc, char *argv[]) 
{
	long int i;

	double cp;


	start_timers();
	init_program(argc, argv);

	init_map(max_map);		/* map version */
	instance.nn_list = compute_nn_lists();

	elapsed = elapsed_time();
    printf("Initialization took %.10f seconds\n",elapsed);


	for ( n_try = 0 ; n_try < max_tries ; n_try++ )
	{
		init_try(n_try);
		
		for ( n_cross = 0 ; n_cross < max_cross ; n_cross++ )
		{

			for ( n_map = 0 ; n_map < max_map ; n_map++ )
			{
				
				start_timers_i();		/* map version */

				compute_total_information();



			cp=ran01(&seed);
			if( cp < c_p ){
				acs_flag=1;
				ras_flag=0;
				printf("acs in \t");
				fprintf(report,"acs in \t");
			}
			else{
				acs_flag=0;
				ras_flag=1;
				printf("ras in \t");
				fprintf(report,"ras in \t");
			}


				while ( !termination_condition() ) 
				{
					
					construct_solutions();

					update_statistics(n_try);

					pheromone_trail_update();  

					iteration++; 
				}

				end_map(n_map);
				
			}
			cross_map();
			
		}
	
		exit_try(n_try);
		
 
	}

	exit_program();

    free( instance.distance );
    free( instance.nn_list );
//    free( pheromone );
    free( heuristic );
    free( total );

    free( time_best_found );
    free( best_found_at );
    free( best_in_try );
    free( time_total_run );
    for ( i = 0 ; i < n_ants ; i++ )
		{
		free( ant[i].tour );
		free( ant[i].visited );
    }
    free( ant );
    free( (*best_so_far_ant).tour );
    free( (*best_so_far_ant).visited );
	free( best_so_far_ant );
    free( prob_of_selection );

	return(0);
}




void init_program( long int argc, char *argv[] ) 
{
	char temp_buffer[LINE_BUF_LEN];

	set_default_parameters();

	best_in_try = calloc(max_tries, sizeof(double));
	best_found_at = calloc(max_tries, sizeof(long int));
	best_found_at_map = calloc(max_map, sizeof(long int));
	time_best_found = calloc(max_tries, sizeof(double));
	time_best_found_map = calloc(max_tries, sizeof(double));
	time_total_run_map = calloc(max_map, sizeof(double));
	time_total_run = calloc(max_tries, sizeof(double));

	instance.nodeptr = read_etsp(name_buf);

	sprintf(temp_buffer,"all.%s",instance.name);
	TRACE ( printf("%s\n",temp_buffer); )
	report = fopen(temp_buffer, "w");

	sprintf(temp_buffer,"try.%s",instance.name);
	TRACE ( printf("%s\n",temp_buffer); )
	report_t = fopen(temp_buffer, "w");

	instance.distance = compute_distances();
	write_params( ); 
	allocate_ants();

  	//pheromone = generate_double_matrix( n, n );
	heuristic = generate_double_matrix( n, n );
	//total = generate_double_matrix( n, n );

	compute_heuristic();

	if ( eas_flag && elitist_ants == 0 )
		elitist_ants = n;


}

void set_default_parameters() 
/*    
      FUNCTION: set default parameter settings
      INPUT:    none
      OUTPUT:   none
      COMMENTS: none
*/
{
    n_ants         = 55;    /* number of ants */
    max_tries      = 100;
	max_iteration  = 500;
    max_time       = 200;
    optimal        = 1;//

	acs_flag	   = TRUE;
	ras_flag       = FALSE;
	eas_flag       = FALSE;
    mmas_flag      = FALSE;

	c_p			   = 1;//probability of acs



	ras_ranks      = 30;


	cross_number   = 3;

	max_cross      = 30;
	max_map        = 20;
 
    alpha          = 1.0;
    beta           = 2.0;
    rho            = 0.5;
    Q			   = 100; 	/* pheromone deposition weight */

	instance.n_near= 30;

//	strncpy( name_buf,"C:\\tspApp\\tsp_data\\oliver30.tsp\\oliver30.tsp", LINE_BUF_LEN);
//	strncpy( name_buf,"C:\\tspApp\\tsp_data\\eil51.tsp\\eil51.tsp", LINE_BUF_LEN);
//	strncpy( name_buf,"C:\\tspApp\\tsp_data\\eil75.tsp\\eil75.tsp", LINE_BUF_LEN);
//	strncpy( name_buf,"C:\\tspApp\\tsp_data\\KroA100.tsp\\KroA100.tsp", LINE_BUF_LEN);
//	strncpy( name_buf,"C:\\tspApp\\tsp_data\\KroA200.tsp\\KroA200.tsp", LINE_BUF_LEN);
	strncpy( name_buf,"C:\\tspApp\\tsp_data\\rat575.tsp\\rat575.tsp", LINE_BUF_LEN);
//	strncpy( name_buf,"C:\\tspApp\\tsp_data\\rat783.tsp\\rat783.tsp", LINE_BUF_LEN);
}


long int ** compute_nn_lists( void )	/*problem maybe exist...*/
/*    
      FUNCTION: computes nearest neighbor lists of depth nn for each city
      INPUT:    none
      OUTPUT:   pointer to the nearest neighbor lists
*/
{
    long int i, node, nn;
    double *distance_vector;
    long int *help_vector;
    long int **m_nnear;
 
    TRACE ( printf("\n computing nearest neighbor lists, "); )

	nn = instance.n_near;
    if((m_nnear = malloc(sizeof(long int) * n * nn
			     + n * sizeof(long int *))) == NULL){
	exit(EXIT_FAILURE);
    }
    distance_vector = calloc(n, sizeof(double));
    help_vector = calloc(n, sizeof(long int));
 
    for ( node = 0 ; node < n ; node++ )
	{/* compute cnd-sets for all node */
		m_nnear[node] = (long int *)(m_nnear + n) + node * nn;

		for ( i = 0 ; i < n ; i++ )
		{  /* Copy distances from nodes to the others */
			distance_vector[i] = instance.distance[node][i];
			help_vector[i] = i;
		}
		distance_vector[node] = LONG_MAX;  /* city is not nearest neighbour */
		sort2(distance_vector, help_vector, 0, n-1);
		for ( i = 0 ; i < nn ; i++ ) 
		{
			m_nnear[node][i] = help_vector[i];
		}
    }
    free(distance_vector);
    free(help_vector);
    TRACE ( printf("\n    .. done\n"); )
    return m_nnear;
}



double ** compute_distances(void)
/*    
      FUNCTION: computes the matrix of all intercity distances
      INPUT:    none
      OUTPUT:   pointer to distance matrix, has to be freed when program stops
*/
{
    long int     i, j;
    double     **matrix;

    if((matrix = malloc(sizeof(double) * n * n +
			sizeof(double *) * n	 )) == NULL){
	fprintf(stderr,"Out of memory_matrix, exit.");
	exit(1);
    }
    for ( i = 0 ; i < n ; i++ ) {
	matrix[i] = (double *)(matrix + n) + i*n;
	for ( j = 0  ; j < n ; j++ ) {
	    matrix[i][j] = distance(i, j);
	}
    }
    return matrix;
	free(matrix);
}

void allocate_ants ( void )
/*    
      FUNCTION:       allocate the memory for the ant colony, the best-so-far and 
                      the iteration best ant
      INPUT:          none
      OUTPUT:         none
      (SIDE)EFFECTS:  allocation of memory for the ant colony and two ants that 
                      store intermediate tours

*/
{
    long int i;
  
    if((ant = malloc(sizeof( ant_struct ) * n_ants +
		     sizeof(ant_struct *) * n_ants	 )) == NULL){
	printf("Out of memory_ant, exit.");
	exit(1);
    }
    for ( i = 0 ; i < n_ants ; i++ ) {
	ant[i].tour        = calloc(n+1, sizeof(long int));
	ant[i].visited     = calloc(n, sizeof(char));
    }

    if((best_so_far_ant = malloc(sizeof( ant_struct ) )) == NULL){
	printf("Out of memory_bestsofarant, exit.");
	exit(1);
    }
    (*best_so_far_ant).tour        = calloc(n+1, sizeof(long int));
    (*best_so_far_ant).visited     = calloc(n, sizeof(char));

    if((prob_of_selection = malloc(sizeof( double ) * instance.n_near )) == NULL){
	printf("Out of memory_prob, exit.");
	exit(1);
    }
}


void init_try( long int ntry ) 
/*    
      FUNCTION: initilialize variables appropriately when starting a trial
      INPUT:    trial number
      OUTPUT:   none
      COMMENTS: none
*/
{ 


	start_timers();
    elapsed = elapsed_time();

	n_map = 0;
    iteration    = 1;         
    (*best_so_far_ant).tour_length = INFTY;
	best_in_try[ntry]= INFTY;



	if ( !(acs_flag || mmas_flag) ) {
		trail_0 = 1. / ( (rho) * nn_tour() );
		/* in the original papers on Ant System, Elitist Ant System, and
			Rank-based Ant System it is not exactly defined what the
			initial value of the pheromones is. Here we set it to some
			small constant, analogously as done in MAX-MIN Ant System.  
		*/
		init_pheromone_trails( trail_0 );
    } 
    
    if ( mmas_flag ) {
		trail_max = 1. / ( (rho) * nn_tour() );
		trail_min = trail_max / ( 2. * n );
		init_pheromone_trails( trail_max );   
    }

    if ( acs_flag ) {
		trail_0 = 1. / ( (double) n * (double) nn_tour( ) ) ;
		init_pheromone_trails( trail_0 );
    }
  
    /* Calculate combined information pheromone times heuristic information */
    //compute_total_information();
    


    fprintf(report,"begin try %li \n",ntry);
}


double nn_tour( void )
/*    
      FUNCTION:       generate some nearest neighbor tour and compute tour length
      INPUT:          none
      OUTPUT:         none
      (SIDE)EFFECTS:  needs ant colony and one statistic ants
*/
{
    long int phase;
	double   help;

    ant_empty_memory( &ant[0] );

    phase = 0; /* counter of the construction steps */
    place_ant( &ant[0], phase);

    while ( phase < n-1 ) {
	phase++;
	choose_closest_next( &ant[0],phase);
    }
    phase = n;
    ant[0].tour[n] = ant[0].tour[0];
/*   copy_from_to( &ant[0], best_so_far_ant ); */
    ant[0].tour_length = compute_tour_length( ant[0].tour );
    help = ant[0].tour_length;
    ant_empty_memory( &ant[0] );
    return help;
}

void init_pheromone_trails( double initial_trail )
/*    
      FUNCTION:      initialize pheromone trails
      INPUT:         initial value of pheromone trails "initial_trail"
      OUTPUT:        none
      (SIDE)EFFECTS: pheromone matrix is reinitialized
*/
{
    long int i, j;

    TRACE ( printf(" init trails with %.15f\n",initial_trail); );

    /* Initialize pheromone trails */
	for ( n_map = 0 ; n_map < max_map ; n_map++){
    for ( i = 0 ; i < n ; i++ ) {
	for ( j =0 ; j <= i ; j++ ) {
	    m_map[n_map].pheromone[i][j] = initial_trail;
	    m_map[n_map].pheromone[j][i] = initial_trail;
	    m_map[n_map].total[i][j] = initial_trail;
	    m_map[n_map].total[j][i] = initial_trail;
	}
    }
	}
}

void compute_total_information( void )
/*    
      FUNCTION: calculates heuristic info times pheromone for each arc
      INPUT:    none  
      OUTPUT:   none
*/
{

    long int     i, j;


    TRACE ( printf("compute total information\n"); );
		jump_flag=0;

    for ( i = 0 ; i < n ; i++ ) {
	for ( j = 0 ; j < i ; j++ ) {
	    m_map[n_map].total[i][j] = pow(m_map[n_map].pheromone[i][j],alpha) * pow(heuristic[i][j],beta);
	    m_map[n_map].total[j][i] = m_map[n_map].total[i][j];
	}
    }
}

void compute_heuristic()
{
	long int i, j;

    for ( i = 0 ; i < n ; i++ )
	{
		for ( j =0 ; j <= i ; j++ )
		{
			heuristic[i][j] = (1.0 / ((double) instance.distance[i][j] + 0.1));
			heuristic[j][i] = heuristic[i][j];
		}
    }
}

void place_ant( ant_struct *a , long int step )
/*    
      FUNCTION:      place an ant on a randomly chosen initial city
      INPUT:         pointer to ant and the number of construction steps 
      OUTPUT:        none
      (SIDE)EFFECT:  ant is put on the chosen city
*/
{
    long int     rnd;

    rnd = (long int) (ran01( &seed ) * (double) n); /* random number between 0 .. n-1 */
    (*a).tour[step] = rnd; 
    (*a).visited[rnd] = TRUE;
}

void choose_closest_next( ant_struct *a, long int phase )
/*    
      FUNCTION:      Choose-
	  for an ant the closest city as the next one 
      INPUT:         pointer to ant and the construction step "phase" 
      OUTPUT:        none 
      (SIDE)EFFECT:  ant moves to the chosen city
*/
{ 
    long int city, current_city, next_city;
	double min_distance;
  
    next_city = n;
    DEBUG( assert ( phase > 0 && phase < n ); );
    current_city = (*a).tour[phase-1];
    min_distance = INFTY;             /* Search shortest edge */    
    for ( city = 0 ; city < n ; city++ ) {
	if ( (*a).visited[city] ) 
	    ; /* city already visited */
	else {
	    if ( instance.distance[current_city][city] < min_distance) {
		next_city = city;
		min_distance = instance.distance[current_city][city];
	    }
	} 
    }
    DEBUG( assert ( 0 <= next_city && next_city < n); );
    (*a).tour[phase] = next_city;
    (*a).visited[next_city] = TRUE;
}

void ant_empty_memory( ant_struct *a ) 
/*    
      FUNCTION:       empty the ants's memory regarding visited cities
      INPUT:          ant identifier
      OUTPUT:         none
      (SIDE)EFFECTS:  vector of visited cities is reinitialized to FALSE
*/
{
    long int   i;

    for( i = 0 ; i < n ; i++ ) {
	(*a).visited[i]=FALSE;
    }
}


long int termination_condition( void )
/*    
      FUNCTION:       checks whether termination condition is met 
      INPUT:          none
      OUTPUT:         0 if condition is not met, number neq 0 otherwise
      (SIDE)EFFECTS:  none
*/
{
	return ( ((iteration >= max_iteration) || (elapsed_time_i() >= max_time)) || 
	  ((*best_so_far_ant).tour_length <= optimal)); 

/******************************stable state of the tour************************
	long int flag;
	
	if(jump_out==(*best_so_far_ant).tour_length)
		jump_flag++;
	else
		jump_out=(*best_so_far_ant).tour_length;

	if(jump_flag>=5)
		flag=1;
	else
		flag=0;


	return(flag);

******************************stable state of the tour************************/

}


void construct_solutions( void )
/*    
      FUNCTION:       manage the solution construction phase
      INPUT:          none
      OUTPUT:         none
      (SIDE)EFFECTS:  when finished, all ants of the colony have constructed a solution  
*/
{
    long int k;        /* counter variable */
    long int step;    /* counter of the number of construction steps */

    /* Mark all cities as unvisited */
    for ( k = 0 ; k < n_ants ; k++)
	{
		ant_empty_memory( &ant[k] );
    }
    
    step = 0; 
    /* Place the ants on same initial city */
    for ( k = 0 ; k < n_ants ; k++ )
		place_ant( &ant[k], step);

    while ( step < n-1 )
	{
		step++;
		for ( k = 0 ; k < n_ants ; k++ )
		{
		   move_to_next( &ant[k], step);

		   if ( acs_flag )
			local_acs_pheromone_update( &ant[k], step );
		}
    }

    step = n;
    for ( k = 0 ; k < n_ants ; k++ )
	{
		ant[k].tour[n] = ant[k].tour[0];
		ant[k].tour_length = compute_tour_length( ant[k].tour );

		if ( acs_flag )
	    local_acs_pheromone_update( &ant[k], step );
	}
}



void move_to_next( ant_struct *a , long int phase )
/*    
      FUNCTION:      Choose for an ant probabilistically a next city among all 
                     unvisited cities in the current city's candidate list. 
		     If this is not possible, choose the closest next
      INPUT:         pointer to ant the construction step "phase" 
      OUTPUT:        none 
      (SIDE)EFFECT:  ant moves to the chosen city
*/
{ 
    long int i, help; 
    long int current_city;
    double   rnd, partial_sum = 0., sum_prob = 0.0;
    /*  double   *prob_of_selection; */ /* stores the selection probabilities 
	of the nearest neighbor cities */
    double   *prob_ptr;

    prob_ptr = prob_of_selection;

    current_city = (*a).tour[phase-1]; /* current_city city of ant k */
	for ( i = 0 ; i < instance.n_near ; i++ )
	{
		if ( (*a).visited[instance.nn_list[current_city][i]] ) 
		    prob_ptr[i] = 0.0;   /* city already visited */
		else
		{
			prob_ptr[i] = m_map[n_map].total[current_city][instance.nn_list[current_city][i]];
		    sum_prob += prob_ptr[i];
		} 
    }

    if (sum_prob <= 0.0)
	{
	/* All cities from the candidate set are tabu */
	choose_best_next( a, phase );
    }     
    else
	{

		/* at least one neighbor is eligible, chose one according to the
		   selection probabilities */
		rnd = ran01( &seed );
		rnd *= sum_prob;
		i = 0;
		partial_sum = prob_ptr[i];
		while ( partial_sum <= rnd )
		{
			i++;
			partial_sum += prob_ptr[i]; 
		}
		help = instance.nn_list[current_city][i];
		(*a).tour[phase] = help; /* instance.nn_list[current_city][i]; */
		(*a).visited[help] = TRUE;
    }
}


void choose_best_next( ant_struct *a, long int phase )
/*    
      FUNCTION:      chooses for an ant as the next city the one with
                     maximal value of heuristic information times pheromone 
      INPUT:         pointer to ant and the construction step
      OUTPUT:        none 
      (SIDE)EFFECT:  ant moves to the chosen city
*/
{ 
    long int city, current_city, next_city;
    double   value_best;

    next_city = n;
    DEBUG( assert ( phase > 0 && phase < n ); );
    current_city = (*a).tour[phase-1];
    value_best = -1.;             /* values in total matrix are always >= 0.0 */    
    for ( city = 0 ; city < n ; city++ ) {
	if ( (*a).visited[city] ) 
	    ; /* city already visited, do nothing */
	else {
	    if ( m_map[n_map].total[current_city][city] > value_best ) {
		next_city = city;
		value_best = m_map[n_map].total[current_city][city];
	    }
	} 
    }
    DEBUG( assert ( 0 <= next_city && next_city < n); );
    DEBUG( assert ( value_best > 0.0 ); )
    DEBUG( assert ( (*a).visited[next_city] == FALSE ); )
    (*a).tour[phase] = next_city;
    (*a).visited[next_city] = TRUE;
}

void update_statistics( long int n_try )
/*    
      FUNCTION:       manage some statistical information about the trial, especially
                      if a new best solution (best-so-far or restart-best) is found and
                      adjust some parameters if a new best solution is found
      INPUT:          none
      OUTPUT:         none
      (SIDE)EFFECTS:  restart-best and best-so-far ant may be updated; trail_min 
                      and trail_max used by MMAS may be updated
*/
{

    long int iteration_best_ant;

    iteration_best_ant = find_best(); /* iteration_best_ant is a global variable */

	

    if ( ant[iteration_best_ant].tour_length < (*best_so_far_ant).tour_length )
	{

		if ( mmas_flag ) {
	    trail_max = 1. / ( (rho) * (*best_so_far_ant).tour_length );
		trail_min = trail_max / ( 2. * n );
		trail_0 = trail_max;
		}
		elapsed = elapsed_time(); /* best sol found after time_used */
		elapsed_map = elapsed_time_i();
		copy_from_to( &ant[iteration_best_ant], best_so_far_ant );
		iter_to_best = iteration;
		restart_found_best = iteration;
	}

}



void copy_from_to(ant_struct *a1, ant_struct *a2) 
{
/*    
      FUNCTION:       copy solution from ant a1 into ant a2
      INPUT:          pointers to the two ants a1 and a2 
      OUTPUT:         none
      (SIDE)EFFECTS:  a2 is copy of a1
*/
    int   i;
  
    (*a2).tour_length = (*a1).tour_length;
    for ( i = 0 ; i < n ; i++ ) {
	(*a2).tour[i] = (*a1).tour[i];
    }
    (*a2).tour[n] = (*a2).tour[0];
}


void pheromone_trail_update( void )  
/*    
      FUNCTION:       manage global pheromone trail update for the ACO algorithms
      INPUT:          none
      OUTPUT:         none
      (SIDE)EFFECTS:  pheromone trails are evaporated and pheromones are deposited 
                      according to the rules defined by the various ACO algorithms.
*/
{


	    evaporation();

		if ( ras_flag )
			ras_update();
    
		else if ( acs_flag )
			acs_global_update();
		
		else if ( eas_flag )
			eas_update();

		else if ( mmas_flag )
			mmas_update();


		else
		{as_update(); 
		}

	    compute_total_information();
}


void evaporation( void )
/*    
      FUNCTION:      implements the pheromone trail evaporation
      INPUT:         none
      OUTPUT:        none
      (SIDE)EFFECTS: pheromones are reduced by factor rho
*/
{ 
    long int    i, j;

    TRACE ( printf("pheromone evaporation\n"); );

    for ( i = 0 ; i < n ; i++ ) {
	for ( j = 0 ; j <= i ; j++ ) {
	    m_map[n_map].pheromone[i][j] = (1 - rho) * m_map[n_map].pheromone[i][j];
	    m_map[n_map].pheromone[j][i] = m_map[n_map].pheromone[i][j];
	}
    }
}

void ras_update( void )
/*    
      FUNCTION:       manage global pheromone deposit for Rank-based Ant System
      INPUT:          none
      OUTPUT:         none
      (SIDE)EFFECTS:  the ras_ranks-1 best ants plus the best-so-far ant deposit pheromone 
                      on matrix "pheromone"
      COMMENTS:       this procedure could be implemented slightly faster, but it is 
                      anyway not critical w.r.t. CPU time given that ras_ranks is 
		      typically very small.
*/
{
    long int i, k, target;
	double b;
    double  *help_b;

    TRACE ( printf("Rank-based Ant System pheromone deposit\n"); );

    help_b = malloc( n_ants  * sizeof(double) );
    for ( k = 0 ; k < n_ants ; k++ )
	help_b[k] = ant[k].tour_length;

    for ( i = 0 ; i < ras_ranks-1 ; i++ ) {
	b = help_b[0]; target = 0;
	for ( k = 0 ; k < n_ants ; k++ ) {
	    if ( help_b[k] < b ) {
		b = help_b[k]; target = k;
	    }
	}
	help_b[target] = LONG_MAX;
	global_update_pheromone_weighted( &ant[target], ras_ranks-i-1 );
    }
    global_update_pheromone_weighted( best_so_far_ant, ras_ranks );
    free ( help_b );
}




void as_update( void )
/*    
      FUNCTION:       manage global pheromone deposit for Ant System
      INPUT:          none
      OUTPUT:         none
      (SIDE)EFFECTS:  all ants deposit pheromones on matrix "pheromone"
*/
{
    long int   k;

    TRACE ( printf("Ant System pheromone deposit\n"); );

    for ( k = 0 ; k < n_ants ; k++ )
	global_update_pheromone( &ant[k] );
}


void local_acs_pheromone_update( ant_struct *a, long int phase )
/*    
      FUNCTION:      removes some pheromone on edge just passed by the ant
      INPUT:         pointer to ant and number of constr. phase
      OUTPUT:        none
      (SIDE)EFFECTS: pheromones of arcs in ant k's tour are increased
      COMMENTS:      I did not do experiments with with different values of the parameter 
                     xi for the local pheromone update; therefore, here xi is fixed to 0.1 
		     as suggested by Gambardella and Dorigo for the TSP. If you wish to run 
		     experiments with that parameter it may be reasonable to use it as a 
		     commandline parameter
*/
{  
    long int  h, j;
    
    DEBUG ( assert ( phase > 0 && phase <= n ); )
	j = (*a).tour[phase];

    h = (*a).tour[phase-1];
    DEBUG ( assert ( 0 <= j && j < n ); )
	DEBUG ( assert ( 0 <= h && h < n ); )
	/* still additional parameter has to be introduced */
	m_map[n_map].pheromone[h][j] = (1. - 0.1) * (m_map[n_map].pheromone[h][j]) + 0.1 * trail_0;
    m_map[n_map].pheromone[j][h] = m_map[n_map].pheromone[h][j];
    m_map[n_map].total[h][j] = pow(m_map[n_map].pheromone[h][j], alpha) * pow(HEURISTIC(h,j),beta);
    m_map[n_map].total[j][h] = m_map[n_map].total[h][j];
}


void acs_global_update( void )
/*    
      FUNCTION:       manage global pheromone deposit for Ant Colony System
      INPUT:          none
      OUTPUT:         none
      (SIDE)EFFECTS:  the best-so-far ant deposits pheromone on matrix "pheromone"
      COMMENTS:       global pheromone deposit in ACS is done per default using 
                      the best-so-far ant; Gambardella & Dorigo examined also iteration-best
		      update (see their IEEE Trans. on Evolutionary Computation article), 
		      but did not use it for the published computational results.
*/
{
    TRACE ( printf("Ant colony System global pheromone deposit\n"); );

    global_acs_pheromone_update( best_so_far_ant );
}


void global_acs_pheromone_update( ant_struct *a )
/*    
      FUNCTION:      reinforces the edges used in ant's solution as in ACS
      INPUT:         pointer to ant that updates the pheromone trail 
      OUTPUT:        none
      (SIDE)EFFECTS: pheromones of arcs in ant k's tour are increased
*/
{  
    long int i, j, h;
    double   d_tau;

    TRACE ( printf("acs specific: global pheromone update\n"); );

    d_tau = 1.0 / (double) (*a).tour_length;

    for ( i = 0 ; i < n ; i++ ) {
	j = (*a).tour[i];
	h = (*a).tour[i+1];

	m_map[n_map].pheromone[j][h] = ((1. - rho) * m_map[n_map].pheromone[j][h]) + rho * d_tau;
	m_map[n_map].pheromone[h][j] = m_map[n_map].pheromone[j][h];

	m_map[n_map].total[h][j] = pow(m_map[n_map].pheromone[h][j], alpha) * pow(HEURISTIC(h,j),beta);
	m_map[n_map].total[j][h] = m_map[n_map].total[h][j];
    }
}


void eas_update( void )
/*    
      FUNCTION:       manage global pheromone deposit for Elitist Ant System
      INPUT:          none
      OUTPUT:         none
      (SIDE)EFFECTS:  all ants plus elitist ant deposit pheromones on matrix "pheromone"
*/
{
    long int   k;

    TRACE ( printf("Elitist Ant System pheromone deposit\n"); );

    for ( k = 0 ; k < n_ants ; k++ )
	global_update_pheromone( &ant[k] );
    global_update_pheromone_weighted( best_so_far_ant, elitist_ants );
}

void mmas_update( void )
/*    
      FUNCTION:       manage global pheromone deposit for MAX-MIN Ant System
      INPUT:          none
      OUTPUT:         none
      (SIDE)EFFECTS:  either the iteration-best or the best-so-far ant deposit pheromone 
                      on matrix "pheromone"
*/
{
    /* we use default upper pheromone trail limit for MMAS and hence we
       do not have to worry regarding keeping the upper limit */

    long int iteration_best_ant;

    TRACE ( printf("MAX-MIN Ant System pheromone deposit\n"); );

    if ( iteration % u_gb ) {
	iteration_best_ant = find_best();
	global_update_pheromone( &ant[iteration_best_ant] );
    }
    else {
        if ( u_gb == 1 && (iteration - restart_found_best > 50))
	    global_update_pheromone( best_so_far_ant );
        else 
	    ;
    }
	u_gb = 25;
}



void global_update_pheromone( ant_struct *a )
/*    
      FUNCTION:      reinforces edges used in ant k's solution
      INPUT:         pointer to ant that updates the pheromone trail 
      OUTPUT:        none
      (SIDE)EFFECTS: pheromones of arcs in ant k's tour are increased
*/
{  
    long int i, j, h;
    double   d_tau;

    TRACE ( printf("global pheromone update\n"); );

    d_tau = Q / (double) (*a).tour_length;
    for ( i = 0 ; i < n ; i++ ) {
	j = (*a).tour[i];
	h = (*a).tour[i+1];
	m_map[n_map].pheromone[j][h] += d_tau;
	m_map[n_map].pheromone[h][j] = m_map[n_map].pheromone[j][h];
    }
}



void global_update_pheromone_weighted( ant_struct *a, long int weight )
/*    
      FUNCTION:      reinforces edges of the ant's tour with weight "weight"
      INPUT:         pointer to ant that updates pheromones and its weight  
      OUTPUT:        none
      (SIDE)EFFECTS: pheromones of arcs in the ant's tour are increased
*/
{  
    long int      i, j, h;
    double        d_tau;

    TRACE ( printf("global pheromone update weighted\n"); );

    d_tau = (double) weight / (double) (*a).tour_length;
    for ( i = 0 ; i < n ; i++ ) {
	j = (*a).tour[i];
	h = (*a).tour[i+1];
	m_map[n_map].pheromone[j][h] += d_tau;
	m_map[n_map].pheromone[h][j] = m_map[n_map].pheromone[j][h];
    }       
}






void exit_try( long int ntry ) 
/*    
      FUNCTION:       save some statistical information on a trial once it finishes
      INPUT:          trial number
      OUTPUT:         none
      COMMENTS:       
*/
{
  checkTour( (*best_so_far_ant).tour );
/*    printTourFile( (*best_so_far_ant).tour ); */

//  best_in_try[ntry] = (*best_so_far_ant).tour_length;
  time_best_found[ntry] = elapsed_t;
//  best_found_at[ntry] = iter_to_best;

  time_total_run[ntry] = elapsed_time();
  printf("In try%ld\t Best: %lf\t Cross:%ld\t Time%lf\t Tot.time%lf\n",ntry, best_in_try[ntry], best_found_at[ntry], time_best_found[ntry],time_total_run[ntry]);
  fprintf(report,"In try%ld\t Best: %lf\t Cross:%ld\t Time%lf\t Tot.time%lf\n",ntry, best_in_try[ntry], best_found_at[ntry], time_best_found[ntry],time_total_run[ntry]);
  fprintf(report_t,"In try%ld\t Best: %lf\t Cross:%ld\t Time%lf\t Tot.time%lf\n",ntry, best_in_try[ntry], best_found_at[ntry], time_best_found[ntry],time_total_run[ntry]);
 
  fprintf(report,"end try %ld\n\n",ntry);
 
  TRACE (output_solution();)
  fflush(report); 
  fflush(report_t);

  
   
 }


void exit_program( void ) 
/*    
      FUNCTION:       save some final statistical information on a trial once it finishes
      INPUT:          none
      OUTPUT:         none
      COMMENTS:       
*/
{
  double best_tour_length, worst_tour_length;
  double   t_avgbest, t_stdbest, t_avgtotal, t_stdtotal;
  double   avg_sol_quality = 0.0, avg_iter_to_bst = 0.0, stddev_best, stddev_iterations;

  best_tour_length = best_of_vector( best_in_try ,max_tries );
  worst_tour_length = worst_of_vector( best_in_try , max_tries );/**/

  avg_iter_to_bst = mean( best_found_at , max_tries );
  stddev_iterations = std_deviation( best_found_at, max_tries, avg_iter_to_bst );

  avg_sol_quality = meanr( best_in_try , max_tries );
  stddev_best = std_deviationr( best_in_try, max_tries, avg_sol_quality);

  t_avgbest = meanr( time_best_found, max_tries );
  t_stdbest = std_deviationr( time_best_found, max_tries, t_avgbest);

  t_avgtotal = meanr( time_total_run, max_tries );
  t_stdtotal = std_deviationr( time_total_run, max_tries, t_avgtotal);
 
  printf(" t_avgbest = %f\n", t_avgbest );
  printf(" t_avgtotal = %f\n", t_avgtotal );

  fprintf(report,"\nBest try: %lf\t\t Worst try: %lf\n", best_tour_length, worst_tour_length);
  fprintf(report,"\nAverage-Best: %.2f\t Stddev-Best: %.2f", avg_sol_quality, stddev_best);
  fprintf(report,"\nAverage-Iterations: %.2f \t Stddev Iterations: %.2f", avg_iter_to_bst, stddev_iterations);
  fprintf(report,"\nAvg.time-best: %.2f stddev.time-best: %.2f\n", t_avgbest, t_stdbest);  
  fprintf(report,"\nAvg.time-total: %.2f stddev.time-total: %.2f\n", t_avgtotal, t_stdtotal); 

  fprintf(report_t,"\nBest try: %lf\t\t Worst try: %lf\n", best_tour_length, worst_tour_length);
  fprintf(report_t,"\nAverage-Best: %.2f\t Stddev-Best: %.2f", avg_sol_quality, stddev_best);
  fprintf(report_t,"\nAverage-Iterations: %.2f \t Stddev Iterations: %.2f", avg_iter_to_bst, stddev_iterations);
  fprintf(report_t,"\nAvg.time-best: %.2f stddev.time-best: %.2f\n", t_avgbest, t_stdbest);  
  fprintf(report_t,"\nAvg.time-total: %.2f stddev.time-total: %.2f\n", t_avgtotal, t_stdtotal); 

  fflush(report);
  fflush(report_t);
}