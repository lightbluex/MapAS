#include<stdio.h>
#include<stdlib.h>
#include <math.h>
#include<time.h>

#include"MapCross.h"
#include"utilities.h"
#include"AS.h"
#include"IO.h"

map_struct *m_map;

long int n_cross;
long int n_map;

long int *m;
double ** matrix;

long int max_cross;
long int max_map;
long int cross_number;

long int route_number;

double elapsed_map;
double *time_best_found_map;

long int *best_found_at_map;
double *time_total_run_map;


void init_map( long int nmap )
/*    
      FUNCTION: initial the maps' statement
      INPUT:    unknown yet
      OUTPUT:   none
*/
{
	long int i,j;

	/*	matrix=generate_double_matrix(n,n);*/
	if((matrix = malloc(sizeof(double) * n * n +
		      sizeof(double *) * n	 )) == NULL){
    printf("Out of memory_matrix_m_cross, exit.");
	}
	for ( i = 0 ; i < n ; i++ ) {
    matrix[i] = (double *)(matrix + n) + i*n;
	}

	if((m_map = malloc(sizeof( map_struct ) * nmap +
		sizeof(map_struct *) * nmap	 )) == NULL)
	{
		printf("Out of memory_map, exit.");
		exit(1);
    }

    for ( i = 0 ; i < nmap ; i++ ) 
	{
		m_map[i].pheromone = generate_double_matrix( n, n );
		m_map[i].total = generate_double_matrix( n, n );

    }
	
	route_number=0;
	for ( i = 0 ; i < n ; i++ ) {
	for ( j = 0 ; j < i ; j++ ) {
	   route_number++;
	}
    }
}





void cross_map()
/*    
      FUNCTION: crossover the pheromone on maps
      INPUT:    none
      OUTPUT:   none
*/
{	
	double sum_fitness=0;

	double ***son_pheromone;

	long int i,j;
	long int y,z;
	
	
//	long int m1;
//	long int m2;

	m = calloc(cross_number, sizeof(long int));

	if((son_pheromone = malloc(sizeof(double)* max_map * n * n )) == NULL){
				  printf("Out of memory_sonphe, exit.");}	
	 
	for(i = 0 ; i < max_map ; i++)
	{
		m_map[i].ave_pheromone=0;
			
		for ( y = 0 ; y < n ; y++ ) {
		for ( z = 0 ; z < y ; z++ ) {
			m_map[i].ave_pheromone+=m_map[i].pheromone[y][z];
		}
		}
		m_map[i].ave_pheromone=(m_map[i].ave_pheromone)/route_number;
		m_map[i].weight=0.5/(m_map[i].ave_pheromone);

		m_map[i].fitness = (1./m_map[i].mapbest)*(1./m_map[i].mapbest)* (1./m_map[i].mapbest)*(1./m_map[i].mapbest)* (1./m_map[i].mapbest)*(1./m_map[i].mapbest)* (1./m_map[i].mapbest)*(1./m_map[i].mapbest);	//Try to enlarge the pheromone xia
		sum_fitness+=m_map[i].fitness;
/*		son_pheromone[i]=generate_double_matrix(n,n);*/
	}
	
	for(i = 0 ; i < max_map ; i++)
	{
		m_map[i].prob=(m_map[i].fitness)/(sum_fitness);
		XIADB(printf("map %ld\tbest %lf\tfit %lf\tprob %lf\t\n",i,m_map[i].mapbest,m_map[i].fitness,m_map[i].prob);)
		XIADB(fprintf(report,"map %ld\tbest %lf\tfit %lf\tprob %lf\t\n",i,m_map[i].mapbest,m_map[i].fitness,m_map[i].prob);)
	}


	for(i = 0 ; i < max_map ; i++)
	{
		for(j = 0 ; j < cross_number ; j++)
		{
			m[j]=choose_map();
//			XIADB(printf("map%ld*%ld\n",j,m[j]);)
		}
//		m1=choose_map();
//		m2=choose_map();
//		XIADB(printf("parents %ldpair\tp1-%ld\t[0-1]phe%lf\tweight%lf\tp2-%ld\t[0-1]phe%lf\tweight%lf",(i+1),m1,m_map[m1].ave_pheromone,m_map[m1].weight,m2,m_map[m2].ave_pheromone,m_map[m2].weight);)
//		XIADB(fprintf(report,"parents %ldpair\tp1-%ld\t[0-1]phe%lf\tweight%lf\tp2-%ld\t[0-1]phe%lf\tweight%lf",(i+1),m1,m_map[m1].pheromone[0][1],m_map[m1].weight,m2,m_map[m2].pheromone[0][1],m_map[m2].weight);)

		son_pheromone[i]=m_cross_over(cross_number);
/*		son_pheromone[i]=cross_over(m1,m2);         *///important!!

//		XIADB(printf("%lf%\n",son_pheromone[i][0][1]);)
	}
	
	for(i = 0 ; i < max_map ; i++)
	{
		copy_pheromone(i,son_pheromone[i]);
//		XIADB(printf("after\t%lf\n",m_map[i].pheromone[0][1]);)
	}

	free(son_pheromone);
	printf("Cross over %ld times...\n\n",(n_cross+1));
	fprintf(report,"Cross over %ld times...\n\n",(n_cross+1));
}



long int choose_map()
/*    
      FUNCTION: choose a map by the probility from all the maps
      INPUT:    none
      OUTPUT:   the number of this choosed map
*/

{	
	long int i;
	long int rnd;

	double finger;
	double roulette=0;

	finger=ran01(&seed);

	i=0;

	while( roulette < finger && i < max_map )
	{
		rnd = (long int) (ran01( &seed ) * (double) max_map);
		roulette+=m_map[rnd].prob;
		i++;
	}

	return rnd;
}


double ** m_cross_over(long int x)
{
	long int i,j,k;

	
	long int change_number;
	
	change_number = (long int) (ran01( &seed ) * (double) route_number);

//	XIADB(printf("route number:%ld\tchange number:%ld\n",route_number,change_number);)
//	XIADB(fprintf(report,"route number:%ld\tchange number:%ld\n",route_number,change_number);)

	
	for ( i = 0 ; i < n ; i++ ) {
	for ( j = 0 ; j < i ; j++ ) {
		for ( k = 0 ; k < x ; k++ ) {
			matrix[i][j]+=((m_map[m[k]].weight)*(m_map[m[k]].pheromone[i][j]));
		}
	matrix[i][j]=matrix[i][j]/(double)x;
	matrix[j][i]=matrix[j][i];
	}
    }

	return matrix;
}



double ** cross_over( long int y, long int z)
{
	long int i,j;
	double ** matrix;
	
	long int change_number;

/*	matrix=generate_double_matrix(n,n);*/
	if((matrix = malloc(sizeof(double) * n * n +
		      sizeof(double *) * n	 )) == NULL){
    printf("Out of memory_matrix_cross_over, exit.");
	}
	for ( i = 0 ; i < n ; i++ ) {
    matrix[i] = (double *)(matrix + n) + i*n;
	}
	
	change_number = (long int) (ran01( &seed ) * (double) route_number);

	XIADB(printf("route number:%ld\tchange number:%ld\n",route_number,change_number);)
	XIADB(fprintf(report,"route number:%ld\tchange number:%ld\n",route_number,change_number);)

	for ( i = 0 ; i < n ; i++ ) {
	for ( j = 0 ; j < i ; j++ ) {
		matrix[i][j]=(m_map[y].weight*m_map[y].pheromone[i][j]+m_map[z].weight*m_map[z].pheromone[i][j])/2.;
		matrix[j][i]=matrix[i][j];
	}
    }

	return matrix;
	free(matrix);
}



void copy_pheromone(long int z, double ** tem)
{
	long int i,j;
	for ( i = 0 ; i < n ; i++ ) {
	for ( j = 0 ; j < i ; j++ ) {
		m_map[z].pheromone[j][i]=m_map[z].pheromone[i][j]=tem[i][j];
		
	}
    }
}




void end_map( long int nmap )
/*    
      FUNCTION:       save some statistical information on a trial once it finishes
      INPUT:          trial number
      OUTPUT:         none
      COMMENTS:       
*/
{
  checkTour( (*best_so_far_ant).tour );
/*    printTourFile( (*best_so_far_ant).tour ); */

  m_map[nmap].mapbest = (*best_so_far_ant).tour_length;
  best_found_at_map[nmap] = iter_to_best;

  time_best_found_map[nmap]=elapsed_map;
  time_total_run_map[nmap] = elapsed_time_i();

  if(m_map[nmap].mapbest < best_in_try[n_try])
  {
		best_in_try[n_try] = m_map[nmap].mapbest;
		best_found_at[n_try] = n_cross;
		elapsed_t = elapsed; /* best sol found after time_used */
  }

  printf("In map %ld\t Best: %lf\t Iterations: %ld\t BestinTime:  %lf\t MapTotaltime %lf\n",nmap, m_map[nmap].mapbest, best_found_at_map[nmap],time_best_found_map[nmap], time_total_run_map[nmap]);
  fprintf(report,"In map %ld\t Best: %lf\t Iterations: %ld\t BestinTime:  %lf\t MapTotaltime %lf\n",nmap, m_map[nmap].mapbest, best_found_at_map[nmap],time_best_found_map[nmap], time_total_run_map[nmap]);
 
  fprintf(report,"end map %ld\n\n",nmap);
 
  TRACE (output_solution();)
  fflush(report); 


  iteration    = 1;         
  (*best_so_far_ant).tour_length = INFTY;
 }
