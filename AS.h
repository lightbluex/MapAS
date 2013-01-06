#define LINE_BUF_LEN     100
#define INFTY            LONG_MAX
#define TRUE  1
#define FALSE 0

#define HEURISTIC(m,n)     (1.0 / ((double) instance.distance[m][n] + 0.1))
/* add a small constant to avoid division by zero if a distance is 
zero */

struct point {
  double x;
  double y;
};


struct problem
{
  char          name[LINE_BUF_LEN];      	 /* instance name */
  double        optimum;					/* optimal tour length if known, otherwise a bound */
  long int      n;                      /* number of cities */
  long int      n_near;                 /* number of nearest neighbors */
  struct point  *nodeptr;               /* array of structs containing coordinates of nodes */
  double        **distance;			    /* distance matrix: distance[i][j] gives distance between city i und j */
  long int      **nn_list;              /* nearest neighbor list; contains for each node i a
                                           sorted list of n_near nearest neighbors */
};




typedef struct 
{
  long int  *tour;
  char      *visited;
  double  tour_length;
} ant_struct;

//extern clock_t start_time;
extern double elapsed;
extern double elapsed_t;

extern double   *best_in_try;		/* best length of each try */
extern long int *best_found_at;		/* iteration when best length is found */
extern double   *time_best_found;	/* time when best length is found  */
extern double   *time_total_run;	/* total time of a try */
extern long int iter_to_best;
extern long int restart_found_best;/* iteration in which restart-best solution is found */

extern ant_struct *ant;
extern ant_struct *best_so_far_ant;

extern long int jump_out;
extern long int jump_flag;

extern double   *prob_of_selection;

extern long int n_ants;			/* number of ants */

extern long int n_try;

extern long int ras_ranks; 

extern long int acs_flag;
extern long int ras_flag;
extern long int eas_flag;	//
extern long int mmas_flag;	//

extern double c_p;

extern double   trail_max;       /* maximum pheromone trail in MMAS */
extern double   trail_min;       /* minimum pheromone trail in MMAS */
extern long int u_gb;            /* every u_gb iterations update with best-so-far ant */

extern long int elitist_ants;	//

extern long int iteration;         /* iteration counter */
extern long int max_iteration;
extern long int max_tries;         /* maximum number of independent tries */
extern double   max_time;          /* maximal allowed run time of a try  */
extern long int optimal;           /* optimal solution or bound to find */

extern double alpha;         /* importance of trail */
extern double beta;          /* importance of heuristic evaluate */
extern double rho;           /* parameter for evaporation */
extern double Q;			  /* pheromone deposition weight */

extern char name_buf[LINE_BUF_LEN]; /* directory (include file name) of data file */

extern struct problem instance;
extern long int n;          /* number of cities in the instance to be solved */


extern	double   trail_0;	/*pheromone initialization */


//extern double **pheromone;
extern double **heuristic;
extern double **total;
extern double  (*distance)(long int, long int);  /* pointer to function returning distance */


void init_program( long int argc, char *argv[] );
void init_try( long int ntry );
void init_pheromone_trails( double initial_trail );

void set_default_parameters();
double ** compute_distances(void);
long int ** compute_nn_lists( void );
void allocate_ants ( void );
double nn_tour( void );
void compute_total_information( void );
void compute_heuristic();
void place_ant( ant_struct *a , long int step );
void choose_closest_next( ant_struct *a, long int phase );
void ant_empty_memory( ant_struct *a );
long int termination_condition( void );
void construct_solutions( void );
void move_to_next( ant_struct *a , long int phase );
void choose_best_next( ant_struct *a, long int phase );
void update_statistics( long int );
void copy_from_to(ant_struct *a1, ant_struct *a2);
void pheromone_trail_update( void );
void evaporation( void );
void as_update( void );
void ras_update( void );
void local_acs_pheromone_update( ant_struct *a, long int phase );
void acs_global_update( void );
void global_acs_pheromone_update( ant_struct *a );
void eas_update( void );
void mmas_update( void );
void global_update_pheromone( ant_struct *a );
void global_update_pheromone_weighted( ant_struct *a, long int weight );
void exit_try( long int ntry );
void exit_program( void );








