#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define MASK 123459876

extern clock_t start_time;
extern clock_t start_time_i;
extern long int seed;

void sort2(double v[], long int v2[], long int left, long int right);
void swap2(double v[], long int v2[], long int i, long int j);
double ** generate_double_matrix( long int n, long int m);
void start_timers();
double elapsed_time();
void start_timers_i();
double elapsed_time_i();



double ran01( long *idum );
double compute_tour_length( long int *t );
long int find_best( void );
void output_solution( void );
void checkTour( long int *t);
void printTour( long int *t );

double mean( long int *values, long int max );
double meanr( double *values, long int max );
double std_deviation( long int *values, long int max, double mean );
double std_deviationr( double *values, long int max, double mean );
double best_of_vector( double *values, long int l );
double worst_of_vector( double *values, long int l );


typedef enum type_timer {REAL, VIRTUAL} TIMER_TYPE;













