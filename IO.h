#define TRACE(x)
#define DEBUG(x)
#define XIADB(x) x

struct point * read_etsp(const char *tsp_file_name);
extern FILE *report,*report_t;

double round_distance (long int i, long int j);
void write_params( void ); 