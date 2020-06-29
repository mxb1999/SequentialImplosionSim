int beam;
int raynum;
int thisx;
int thisz;
//Launch Ray Values
int thisx_0;

int thisx_00;
int thisz_0;
int thisz_00;
int finalt;
double ztarg;
double slope;
double xtarg;
double k;
double knorm;
double cs;
int ray1num;
//Pointers for necessary arrays
double* uray; //nt
double* rayx;//nt
double* rayz;//nt
double* amp_norm;//nt
double** intersections; //nx nz
int** marked; //nx nz numstored nbeams
double** dedendx; //nx nz
double** dedendz; //nx nz
double* x; //nx
double* z; //nz
double** eden; //nx nz

double*** edep; //nx+2 nz+2 nbeams
int*** present; //nx nz nbeams
double** machnum; //nx nz
int**** boxes; //nbeams nrays nx*3 2
bool**** boxTrack;
double*** W1_storage; //nx nz numstored
double*** W2_storage; //nx nz numstored
double** u_flow; //nx nz
double*** dkx; //nbeams nrays 2
double*** dkz; //nbeams nrays 2
double*** dkmag; //nbeams nrays 2
double** W1;//nx nz
double** W2;//nx nz
double** W1_init;//nx nz
//Launch_Ray_XZ specific arrays (all have a length of nt)
double* myx; //nt
double* mytime;//nt
double* myz;//nt
double* mykx;//nt
double* mykz;//nt
double* myvx;//nt
double* myvz;//nt
double* amplitude_norm;//nt
double* markingx;//nt
double* markingz;//nt
double* nuei;//nt
//CBET specific arrays
double** W2_init;//nx nz
double** W1_new;//nx nz
double** W2_new;//nx n
double** i_b1;//nx nz
double** i_b2;//nx nz
double** i_b1_new;//nx nz
double** i_b2_new;//nx nz
double** wpe; //nx nz
double*** crossesz; //nbeams nrays ncrossings
double*** crossesx; //nbeams nrays ncrossings
int*** ints; //nbeams nrays ncrossings
//arrays used only for plotting
double** i_bplot;//nx nz
double** i_b_newplot;//nx nz
double** edenplot; //the array is eden/ncrit,  nx nz
double** edepplot; //nx nz
