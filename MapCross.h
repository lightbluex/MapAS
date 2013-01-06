
typedef struct 
{
	double fitness;
	double prob;
	double ave_best;
	double mapbest;
	double **pheromone;
	double **total;
	double ave_pheromone;
	double weight;
} map_struct;


extern map_struct * m_map;

extern long int n_cross;
extern long int n_map;
extern long int max_cross;
extern long int max_map;
extern long int route_number;

extern long int cross_number;

extern long int *best_found_at_map;
extern double *time_total_run_map;

extern double elapsed_map;
extern double *time_best_found_map;

void init_map( long int nmap );
void cross_map();
void end_map( long int nmap );
long int choose_map();

double ** cross_over( long int n, long int m);
double ** m_cross_over(long int x);
void copy_pheromone();