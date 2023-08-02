#include <stdlib.h> 
#include <string.h> 
#include <math.h>
int niter;
double dt;
int block0np0;
int block0np1;
int block0np2;
double Delta0block0;
double Delta1block0;
double Delta2block0;
double Twall;
double Minf;
double SuthT;
double gama;
double RefT;
double Pr;
double Re;
double gamma_m1;
double inv_0;
double inv_1;
double inv_2;
double inv_3;
double inv_4;
double inv_5;
double beta_0;
double xts;
double k_0;
double omega_0;
double tripA;
double a;
double b;
double H;
double L;
double rcinv6;
double rcinv7;
double rc8;
double rc9;
double rcinv10;
double rcinv11;
double rcinv12;
double rc13;
double rc14;
double rc15;
double rc16;
double rcinv17;
double rcinv18;
double rc19;
double rc20;
double rc21;
double rc22;
double rc23;
double rc24;
double rc25;
double rc26;
double rc27;
double rc28;
double rc29;
double rc30;
double rc31;
double rc32;
double rc33;
double rc34;
double rc35;
double rcinv36;
double rcinv37;
double rc38;
double rcinv39;
#define OPS_3D
#include "ops_seq.h"
#include "opensbliblock00_kernels.h"
int main(int argc, char **argv) 
{
niter = 200000;
double rkB[] = {1.0/3.0, 15.0/16.0, 8.0/15.0};
double rkA[] = {0, -5.0/9.0, -153.0/128.0};
dt = 0.003;
block0np0 = 400.0;
block0np1 = 200.0;
block0np2 = 100;
Delta0block0 = 400.0/(block0np0-1);
Delta1block0 = 115.0/(block0np1-1);
Delta2block0 = 20.0/(block0np2-1);
Twall = 1.68;
Minf = 2.0;
SuthT = 110.4;
gama = 1.4;
RefT = 288.0;
Pr = 0.72;
Re = 950.0;
gamma_m1 = gama - 1;
inv_0 = 1.0/Delta1block0;
inv_1 = 1.0/Delta2block0;
inv_2 = 1.0/Delta0block0;
inv_3 = pow(Delta2block0, -2);
inv_4 = pow(Delta1block0, -2);
inv_5 = pow(Delta0block0, -2);
beta_0 = 0.41;
xts = 50.0;
k_0 = 0.02;
omega_0 = 0.1;
tripA = 0.2;
a = 20.0;
b = 5.0;
H = 115.0;
L = 400.0;
rcinv6 = 1.0/H;
rcinv7 = pow(a, -2);
rc8 = 1.0/20.0;
rc9 = 1.0/2.0;
rcinv10 = 1.0/(gama - 1);
rcinv11 = pow(Minf, -2);
rcinv12 = 1.0/gama;
rc13 = 1.0/4.0;
rc14 = 4.0/3.0;
rc15 = 25.0/12.0;
rc16 = 1.0/12.0;
rcinv17 = pow(Minf, -2.0);
rcinv18 = 1.0/RefT;
rc19 = -1.0/2.0;
rc20 = 5.0/3.0;
rc21 = 31.0/6.0;
rc22 = 11.0/6.0;
rc23 = 25.0/6.0;
rc24 = 19.0/6.0;
rc25 = 1.0/3.0;
rc26 = 2.0/3.0;
rc27 = 13.0/6.0;
rc28 = 5.0/6.0;
rc29 = 3.0/10.0;
rc30 = 3.0/5.0;
rc31 = 1.0/10.0;
rc32 = 1.0/6.0;
rc33 = 5.0/12.0;
rc34 = 7.0/12.0;
rc35 = 11.0/12.0;
rcinv36 = 1.0/gamma_m1;
rcinv37 = 1.0/Re;
rc38 = 8.0/3.0;
rcinv39 = 1.0/Pr;
// Initializing OPS 
ops_init(argc,argv,1);
ops_decl_const("niter" , 1, "int", &niter);
ops_decl_const("dt" , 1, "double", &dt);
ops_decl_const("block0np0" , 1, "int", &block0np0);
ops_decl_const("block0np1" , 1, "int", &block0np1);
ops_decl_const("block0np2" , 1, "int", &block0np2);
ops_decl_const("Delta0block0" , 1, "double", &Delta0block0);
ops_decl_const("Delta1block0" , 1, "double", &Delta1block0);
ops_decl_const("Delta2block0" , 1, "double", &Delta2block0);
ops_decl_const("Twall" , 1, "double", &Twall);
ops_decl_const("Minf" , 1, "double", &Minf);
ops_decl_const("SuthT" , 1, "double", &SuthT);
ops_decl_const("gama" , 1, "double", &gama);
ops_decl_const("RefT" , 1, "double", &RefT);
ops_decl_const("Pr" , 1, "double", &Pr);
ops_decl_const("Re" , 1, "double", &Re);
ops_decl_const("gamma_m1" , 1, "double", &gamma_m1);
ops_decl_const("inv_0" , 1, "double", &inv_0);
ops_decl_const("inv_1" , 1, "double", &inv_1);
ops_decl_const("inv_2" , 1, "double", &inv_2);
ops_decl_const("inv_3" , 1, "double", &inv_3);
ops_decl_const("inv_4" , 1, "double", &inv_4);
ops_decl_const("inv_5" , 1, "double", &inv_5);
ops_decl_const("beta_0" , 1, "double", &beta_0);
ops_decl_const("xts" , 1, "double", &xts);
ops_decl_const("k_0" , 1, "double", &k_0);
ops_decl_const("omega_0" , 1, "double", &omega_0);
ops_decl_const("tripA" , 1, "double", &tripA);
ops_decl_const("a" , 1, "double", &a);
ops_decl_const("b" , 1, "double", &b);
ops_decl_const("H" , 1, "double", &H);
ops_decl_const("L" , 1, "double", &L);
ops_decl_const("rcinv6" , 1, "double", &rcinv6);
ops_decl_const("rcinv7" , 1, "double", &rcinv7);
ops_decl_const("rc8" , 1, "double", &rc8);
ops_decl_const("rc9" , 1, "double", &rc9);
ops_decl_const("rcinv10" , 1, "double", &rcinv10);
ops_decl_const("rcinv11" , 1, "double", &rcinv11);
ops_decl_const("rcinv12" , 1, "double", &rcinv12);
ops_decl_const("rc13" , 1, "double", &rc13);
ops_decl_const("rc14" , 1, "double", &rc14);
ops_decl_const("rc15" , 1, "double", &rc15);
ops_decl_const("rc16" , 1, "double", &rc16);
ops_decl_const("rcinv17" , 1, "double", &rcinv17);
ops_decl_const("rcinv18" , 1, "double", &rcinv18);
ops_decl_const("rc19" , 1, "double", &rc19);
ops_decl_const("rc20" , 1, "double", &rc20);
ops_decl_const("rc21" , 1, "double", &rc21);
ops_decl_const("rc22" , 1, "double", &rc22);
ops_decl_const("rc23" , 1, "double", &rc23);
ops_decl_const("rc24" , 1, "double", &rc24);
ops_decl_const("rc25" , 1, "double", &rc25);
ops_decl_const("rc26" , 1, "double", &rc26);
ops_decl_const("rc27" , 1, "double", &rc27);
ops_decl_const("rc28" , 1, "double", &rc28);
ops_decl_const("rc29" , 1, "double", &rc29);
ops_decl_const("rc30" , 1, "double", &rc30);
ops_decl_const("rc31" , 1, "double", &rc31);
ops_decl_const("rc32" , 1, "double", &rc32);
ops_decl_const("rc33" , 1, "double", &rc33);
ops_decl_const("rc34" , 1, "double", &rc34);
ops_decl_const("rc35" , 1, "double", &rc35);
ops_decl_const("rcinv36" , 1, "double", &rcinv36);
ops_decl_const("rcinv37" , 1, "double", &rcinv37);
ops_decl_const("rc38" , 1, "double", &rc38);
ops_decl_const("rcinv39" , 1, "double", &rcinv39);
// Define and Declare OPS Block
ops_block opensbliblock00 = ops_decl_block(3, "opensbliblock00");
#include "defdec_data_set.h"
// Define and declare stencils
int stencil_0_22temp[] = {-4, 0, 0, -3, 0, 0, -2, 0, 0, -1, 0, 0, 1, 0, 0, 2, 0, 0, 3, 0, 0, 4, 0, 0};
ops_stencil stencil_0_22 = ops_decl_stencil(3,8,stencil_0_22temp,"stencil_0_22temp");
int stencil_0_01temp[] = {0, -2, 0, 0, -1, 0, 0, 0, -2, 0, 0, -1, 0, 0, 0, 0, 0, 1, 0, 0, 2, 0, 1, 0, 0, 2, 0, 0, 3, 0, 0, 4, 0};
ops_stencil stencil_0_01 = ops_decl_stencil(3,11,stencil_0_01temp,"stencil_0_01temp");
int stencil_0_23temp[] = {0, -3, 0, 0, -2, 0, 0, -1, 0, 0, 1, 0, 0, 2, 0, 0, 3, 0};
ops_stencil stencil_0_23 = ops_decl_stencil(3,6,stencil_0_23temp,"stencil_0_23temp");
int stencil_0_17temp[] = {0, -3, 0, 0, -2, 0, 0, -1, 0, 0, 0, 0, 0, 1, 0};
ops_stencil stencil_0_17 = ops_decl_stencil(3,5,stencil_0_17temp,"stencil_0_17temp");
int stencil_0_08temp[] = {0, -2, 0, 0, -1, 0, 0, 0, 0, 0, 1, 0, 0, 2, 0, 0, 3, 0};
ops_stencil stencil_0_08 = ops_decl_stencil(3,6,stencil_0_08temp,"stencil_0_08temp");
int stencil_0_05temp[] = {0, 0, 0, 1, 0, 0};
ops_stencil stencil_0_05 = ops_decl_stencil(3,2,stencil_0_05temp,"stencil_0_05temp");
int stencil_0_06temp[] = {-2, 0, 0, -1, 0, 0, 0, 0, 0, 1, 0, 0, 2, 0, 0, 3, 0, 0};
ops_stencil stencil_0_06 = ops_decl_stencil(3,6,stencil_0_06temp,"stencil_0_06temp");
int stencil_0_12temp[] = {0, -1, 0, 0, 0, 0};
ops_stencil stencil_0_12 = ops_decl_stencil(3,2,stencil_0_12temp,"stencil_0_12temp");
int stencil_0_21temp[] = {-3, 0, 0, -2, 0, 0, -1, 0, 0, 1, 0, 0, 2, 0, 0, 3, 0, 0};
ops_stencil stencil_0_21 = ops_decl_stencil(3,6,stencil_0_21temp,"stencil_0_21temp");
int stencil_0_04temp[] = {0, -2, 0, 0, -1, 0, 0, 0, 0, 0, 1, 0, 0, 2, 0, 0, 3, 0, 0, 4, 0};
ops_stencil stencil_0_04 = ops_decl_stencil(3,7,stencil_0_04temp,"stencil_0_04temp");
int stencil_0_24temp[] = {0, -4, 0, 0, -3, 0, 0, -2, 0, 0, -1, 0, 0, 1, 0, 0, 2, 0, 0, 3, 0, 0, 4, 0};
ops_stencil stencil_0_24 = ops_decl_stencil(3,8,stencil_0_24temp,"stencil_0_24temp");
int stencil_0_20temp[] = {-4, 0, 0, -3, 0, 0, -2, 0, 0, -1, 0, 0, 0, -2, 0, 0, -1, 0, 0, 0, 0, 0, 1, 0, 0, 2, 0, 0, 3, 0, 0, 4, 0, 1, 0, 0, 2, 0, 0, 3, 0, 0, 4, 0, 0};
ops_stencil stencil_0_20 = ops_decl_stencil(3,15,stencil_0_20temp,"stencil_0_20temp");
int stencil_0_07temp[] = {0, 0, 0, 0, 1, 0};
ops_stencil stencil_0_07 = ops_decl_stencil(3,2,stencil_0_07temp,"stencil_0_07temp");
int stencil_0_19temp[] = {0, -4, 0, 0, -3, 0, 0, -2, 0, 0, -1, 0, 0, 0, 0, 0, 1, 0, 0, 2, 0, 0, 3, 0, 0, 4, 0};
ops_stencil stencil_0_19 = ops_decl_stencil(3,9,stencil_0_19temp,"stencil_0_19temp");
int stencil_0_13temp[] = {-1, 0, 0, 0, 0, 0};
ops_stencil stencil_0_13 = ops_decl_stencil(3,2,stencil_0_13temp,"stencil_0_13temp");
int stencil_0_00temp[] = {0, 0, 0};
ops_stencil stencil_0_00 = ops_decl_stencil(3,1,stencil_0_00temp,"stencil_0_00temp");
int stencil_0_16temp[] = {-1, 0, 0, 0, 0, 0, 1, 0, 0, 2, 0, 0, 3, 0, 0, 4, 0, 0};
ops_stencil stencil_0_16 = ops_decl_stencil(3,6,stencil_0_16temp,"stencil_0_16temp");
int stencil_0_02temp[] = {-4, 0, 0, -3, 0, 0, -2, 0, 0, -1, 0, 0, 0, -2, 0, 0, -1, 0, 0, 0, -2, 0, 0, -1, 0, 0, 0, 0, 0, 1, 0, 0, 2, 0, 1, 0, 0, 2, 0, 0, 3, 0, 0, 4, 0, 1, 0, 0, 2, 0, 0, 3, 0, 0, 4, 0, 0};
ops_stencil stencil_0_02 = ops_decl_stencil(3,19,stencil_0_02temp,"stencil_0_02temp");
int stencil_0_03temp[] = {0, 0, -2, 0, 0, -1, 0, 0, 0, 0, 0, 1, 0, 0, 2};
ops_stencil stencil_0_03 = ops_decl_stencil(3,5,stencil_0_03temp,"stencil_0_03temp");
int stencil_0_18temp[] = {0, -3, 0, 0, -2, 0, 0, -1, 0, 0, 0, 0, 0, 1, 0, 0, 2, 0, 0, 3, 0};
ops_stencil stencil_0_18 = ops_decl_stencil(3,7,stencil_0_18temp,"stencil_0_18temp");
int stencil_0_14temp[] = {-4, 0, 0, -3, 0, 0, -2, 0, 0, -1, 0, 0, 0, 0, 0, 1, 0, 0, 2, 0, 0, 3, 0, 0, 4, 0, 0};
ops_stencil stencil_0_14 = ops_decl_stencil(3,9,stencil_0_14temp,"stencil_0_14temp");
int stencil_0_15temp[] = {0, 0, -2, 0, 0, -1, 0, 0, 1, 0, 0, 2};
ops_stencil stencil_0_15 = ops_decl_stencil(3,4,stencil_0_15temp,"stencil_0_15temp");
int stencil_0_11temp[] = {0, 0, -1, 0, 0, 0};
ops_stencil stencil_0_11 = ops_decl_stencil(3,2,stencil_0_11temp,"stencil_0_11temp");
int stencil_0_09temp[] = {0, 0, 0, 0, 0, 1};
ops_stencil stencil_0_09 = ops_decl_stencil(3,2,stencil_0_09temp,"stencil_0_09temp");
int stencil_0_10temp[] = {0, 0, -2, 0, 0, -1, 0, 0, 0, 0, 0, 1, 0, 0, 2, 0, 0, 3};
ops_stencil stencil_0_10 = ops_decl_stencil(3,6,stencil_0_10temp,"stencil_0_10temp");
#include "bc_exchanges.h"
// Init OPS partition
ops_partition("");

int iteration_range_67_block0[] = {-3, block0np0 + 4, -3, block0np1 + 4, -3, block0np2 + 4};
ops_par_loop(opensbliblock00Kernel067, "Grid_based_initialisation0", opensbliblock00, 3, iteration_range_67_block0,
ops_arg_dat(rhoE_B0, 1, stencil_0_00, "double", OPS_WRITE),
ops_arg_dat(x2_B0, 1, stencil_0_00, "double", OPS_WRITE),
ops_arg_dat(rho_B0, 1, stencil_0_00, "double", OPS_WRITE),
ops_arg_dat(rhou1_B0, 1, stencil_0_00, "double", OPS_WRITE),
ops_arg_dat(rhou0_B0, 1, stencil_0_00, "double", OPS_WRITE),
ops_arg_dat(x0_B0, 1, stencil_0_00, "double", OPS_WRITE),
ops_arg_dat(rhou2_B0, 1, stencil_0_00, "double", OPS_WRITE),
ops_arg_dat(x1_B0, 1, stencil_0_00, "double", OPS_RW),
ops_arg_idx());


int iteration_range_70_block0[] = {0, block0np0, 0, block0np1, 0, block0np2};
ops_par_loop(opensbliblock00Kernel070, "MetricsEquation evaluation", opensbliblock00, 3, iteration_range_70_block0,
ops_arg_dat(x1_B0, 1, stencil_0_20, "double", OPS_READ),
ops_arg_dat(x0_B0, 1, stencil_0_20, "double", OPS_READ),
ops_arg_dat(D00_B0, 1, stencil_0_00, "double", OPS_WRITE),
ops_arg_dat(D01_B0, 1, stencil_0_00, "double", OPS_WRITE),
ops_arg_dat(D11_B0, 1, stencil_0_00, "double", OPS_WRITE),
ops_arg_dat(D10_B0, 1, stencil_0_00, "double", OPS_WRITE),
ops_arg_dat(detJ_B0, 1, stencil_0_00, "double", OPS_WRITE),
ops_arg_idx());


int iteration_range_71_block0[] = {0, 1, -3, block0np1 + 4, -3, block0np2 + 4};
ops_par_loop(opensbliblock00Kernel071, "Metric boundary dir0 side0", opensbliblock00, 3, iteration_range_71_block0,
ops_arg_dat(D00_B0, 1, stencil_0_21, "double", OPS_RW),
ops_arg_dat(D01_B0, 1, stencil_0_21, "double", OPS_RW),
ops_arg_dat(D11_B0, 1, stencil_0_21, "double", OPS_RW),
ops_arg_dat(D10_B0, 1, stencil_0_21, "double", OPS_RW),
ops_arg_dat(detJ_B0, 1, stencil_0_21, "double", OPS_RW));


int iteration_range_72_block0[] = {block0np0 - 1, block0np0, -3, block0np1 + 4, -3, block0np2 + 4};
ops_par_loop(opensbliblock00Kernel072, "Metric boundary dir0 side1", opensbliblock00, 3, iteration_range_72_block0,
ops_arg_dat(D00_B0, 1, stencil_0_22, "double", OPS_RW),
ops_arg_dat(D01_B0, 1, stencil_0_22, "double", OPS_RW),
ops_arg_dat(D11_B0, 1, stencil_0_22, "double", OPS_RW),
ops_arg_dat(D10_B0, 1, stencil_0_22, "double", OPS_RW),
ops_arg_dat(detJ_B0, 1, stencil_0_22, "double", OPS_RW));


int iteration_range_73_block0[] = {-3, block0np0 + 4, 0, 1, -3, block0np2 + 4};
ops_par_loop(opensbliblock00Kernel073, "Metric boundary dir1 side0", opensbliblock00, 3, iteration_range_73_block0,
ops_arg_dat(D00_B0, 1, stencil_0_23, "double", OPS_RW),
ops_arg_dat(D01_B0, 1, stencil_0_23, "double", OPS_RW),
ops_arg_dat(D11_B0, 1, stencil_0_23, "double", OPS_RW),
ops_arg_dat(D10_B0, 1, stencil_0_23, "double", OPS_RW),
ops_arg_dat(detJ_B0, 1, stencil_0_23, "double", OPS_RW));


int iteration_range_74_block0[] = {-3, block0np0 + 4, block0np1 - 1, block0np1, -3, block0np2 + 4};
ops_par_loop(opensbliblock00Kernel074, "Metric boundary dir1 side1", opensbliblock00, 3, iteration_range_74_block0,
ops_arg_dat(D00_B0, 1, stencil_0_24, "double", OPS_RW),
ops_arg_dat(D01_B0, 1, stencil_0_24, "double", OPS_RW),
ops_arg_dat(D11_B0, 1, stencil_0_24, "double", OPS_RW),
ops_arg_dat(D10_B0, 1, stencil_0_24, "double", OPS_RW),
ops_arg_dat(detJ_B0, 1, stencil_0_24, "double", OPS_RW));


ops_halo_transfer(exchange75_block0);
ops_halo_transfer(exchange76_block0);
int iteration_range_81_block0[] = {0, block0np0, 0, block0np1, 0, block0np2};
ops_par_loop(opensbliblock00Kernel081, "MetricsEquation evaluation", opensbliblock00, 3, iteration_range_81_block0,
ops_arg_dat(D00_B0, 1, stencil_0_20, "double", OPS_READ),
ops_arg_dat(D01_B0, 1, stencil_0_20, "double", OPS_READ),
ops_arg_dat(D11_B0, 1, stencil_0_20, "double", OPS_READ),
ops_arg_dat(D10_B0, 1, stencil_0_20, "double", OPS_READ),
ops_arg_dat(SD101_B0, 1, stencil_0_00, "double", OPS_WRITE),
ops_arg_dat(SD100_B0, 1, stencil_0_00, "double", OPS_WRITE),
ops_arg_dat(SD000_B0, 1, stencil_0_00, "double", OPS_WRITE),
ops_arg_dat(SD001_B0, 1, stencil_0_00, "double", OPS_WRITE),
ops_arg_dat(SD010_B0, 1, stencil_0_00, "double", OPS_WRITE),
ops_arg_dat(SD110_B0, 1, stencil_0_00, "double", OPS_WRITE),
ops_arg_dat(SD111_B0, 1, stencil_0_00, "double", OPS_WRITE),
ops_arg_dat(SD011_B0, 1, stencil_0_00, "double", OPS_WRITE),
ops_arg_idx());


double cpu_start0, elapsed_start0;
ops_timers(&cpu_start0, &elapsed_start0);
for(int iter=0; iter<=niter - 1; iter++)
{
if(fmod(iter+1, 50) == 0){
        ops_printf("Iteration is %d\n", iter+1); 
        ops_NaNcheck(rho_B0);
}
int iteration_range_61_block0[] = {0, 1, -3, block0np1 + 4, -3, block0np2 + 4};
ops_par_loop(opensbliblock00Kernel061, "InletTransfer boundary dir0 side0", opensbliblock00, 3, iteration_range_61_block0,
ops_arg_dat(rhou1_B0, 1, stencil_0_13, "double", OPS_RW),
ops_arg_dat(rhoE_B0, 1, stencil_0_13, "double", OPS_RW),
ops_arg_dat(rhou0_B0, 1, stencil_0_13, "double", OPS_RW),
ops_arg_dat(rho_B0, 1, stencil_0_13, "double", OPS_RW),
ops_arg_dat(rhou2_B0, 1, stencil_0_13, "double", OPS_RW));


int iteration_range_62_block0[] = {block0np0 - 1, block0np0, -3, block0np1 + 4, -3, block0np2 + 4};
ops_par_loop(opensbliblock00Kernel062, "Extrapolation boundary dir0 side1", opensbliblock00, 3, iteration_range_62_block0,
ops_arg_dat(rhoE_B0, 1, stencil_0_16, "double", OPS_RW),
ops_arg_dat(rhou1_B0, 1, stencil_0_16, "double", OPS_RW),
ops_arg_dat(rhou0_B0, 1, stencil_0_16, "double", OPS_RW),
ops_arg_dat(rho_B0, 1, stencil_0_16, "double", OPS_RW),
ops_arg_dat(rhou2_B0, 1, stencil_0_16, "double", OPS_RW));


int iteration_range_63_block0[] = {-3, block0np0 + 4, 0, 1, -3, block0np2 + 4};
ops_par_loop(opensbliblock00Kernel063, "ForcingStripWall boundary dir1 side0", opensbliblock00, 3, iteration_range_63_block0,
ops_arg_dat(x2_B0, 1, stencil_0_00, "double", OPS_READ),
ops_arg_dat(x0_B0, 1, stencil_0_00, "double", OPS_READ),
ops_arg_dat(rhoE_B0, 1, stencil_0_17, "double", OPS_RW),
ops_arg_dat(rhou1_B0, 1, stencil_0_18, "double", OPS_RW),
ops_arg_dat(rhou0_B0, 1, stencil_0_18, "double", OPS_RW),
ops_arg_dat(rho_B0, 1, stencil_0_18, "double", OPS_RW),
ops_arg_dat(rhou2_B0, 1, stencil_0_18, "double", OPS_RW),
ops_arg_gbl(&iter, 1, "int", OPS_READ));


int iteration_range_64_block0[] = {-3, block0np0 + 4, block0np1 - 1, block0np1, -3, block0np2 + 4};
ops_par_loop(opensbliblock00Kernel064, "ZeroGradientOutlet boundary dir1 side1", opensbliblock00, 3, iteration_range_64_block0,
ops_arg_dat(rhoE_B0, 1, stencil_0_19, "double", OPS_RW),
ops_arg_dat(rhou1_B0, 1, stencil_0_19, "double", OPS_RW),
ops_arg_dat(rhou0_B0, 1, stencil_0_19, "double", OPS_RW),
ops_arg_dat(rho_B0, 1, stencil_0_19, "double", OPS_RW),
ops_arg_dat(rhou2_B0, 1, stencil_0_19, "double", OPS_RW));


ops_halo_transfer(exchange65_block0);
ops_halo_transfer(exchange66_block0);
for(int stage=0; stage<=2; stage++)
{
int iteration_range_3_block0[] = {-3, block0np0 + 4, -3, block0np1 + 4, -3, block0np2 + 4};
ops_par_loop(opensbliblock00Kernel003, "CRu2", opensbliblock00, 3, iteration_range_3_block0,
ops_arg_dat(rho_B0, 1, stencil_0_00, "double", OPS_READ),
ops_arg_dat(rhou2_B0, 1, stencil_0_00, "double", OPS_READ),
ops_arg_dat(u2_B0, 1, stencil_0_00, "double", OPS_WRITE));


int iteration_range_10_block0[] = {-3, block0np0 + 4, -3, block0np1 + 4, -3, block0np2 + 4};
ops_par_loop(opensbliblock00Kernel010, "CRu0", opensbliblock00, 3, iteration_range_10_block0,
ops_arg_dat(rhou0_B0, 1, stencil_0_00, "double", OPS_READ),
ops_arg_dat(rho_B0, 1, stencil_0_00, "double", OPS_READ),
ops_arg_dat(u0_B0, 1, stencil_0_00, "double", OPS_WRITE));


int iteration_range_5_block0[] = {-3, block0np0 + 4, -3, block0np1 + 4, -3, block0np2 + 4};
ops_par_loop(opensbliblock00Kernel005, "CRu1", opensbliblock00, 3, iteration_range_5_block0,
ops_arg_dat(rhou1_B0, 1, stencil_0_00, "double", OPS_READ),
ops_arg_dat(rho_B0, 1, stencil_0_00, "double", OPS_READ),
ops_arg_dat(u1_B0, 1, stencil_0_00, "double", OPS_WRITE));


int iteration_range_12_block0[] = {0, block0np0, -3, block0np1 + 4, 0, block0np2};
ops_par_loop(opensbliblock00Kernel012, "CRU1", opensbliblock00, 3, iteration_range_12_block0,
ops_arg_dat(u0_B0, 1, stencil_0_00, "double", OPS_READ),
ops_arg_dat(D11_B0, 1, stencil_0_00, "double", OPS_READ),
ops_arg_dat(u1_B0, 1, stencil_0_00, "double", OPS_READ),
ops_arg_dat(D10_B0, 1, stencil_0_00, "double", OPS_READ),
ops_arg_dat(U1_B0, 1, stencil_0_00, "double", OPS_WRITE));


int iteration_range_17_block0[] = {-3, block0np0 + 4, -3, block0np1 + 4, -3, block0np2 + 4};
ops_par_loop(opensbliblock00Kernel017, "CRp", opensbliblock00, 3, iteration_range_17_block0,
ops_arg_dat(rhoE_B0, 1, stencil_0_00, "double", OPS_READ),
ops_arg_dat(rho_B0, 1, stencil_0_00, "double", OPS_READ),
ops_arg_dat(u0_B0, 1, stencil_0_00, "double", OPS_READ),
ops_arg_dat(u1_B0, 1, stencil_0_00, "double", OPS_READ),
ops_arg_dat(u2_B0, 1, stencil_0_00, "double", OPS_READ),
ops_arg_dat(p_B0, 1, stencil_0_00, "double", OPS_WRITE));


int iteration_range_14_block0[] = {-3, block0np0 + 4, 0, block0np1, 0, block0np2};
ops_par_loop(opensbliblock00Kernel014, "CRU0", opensbliblock00, 3, iteration_range_14_block0,
ops_arg_dat(D00_B0, 1, stencil_0_00, "double", OPS_READ),
ops_arg_dat(D01_B0, 1, stencil_0_00, "double", OPS_READ),
ops_arg_dat(u0_B0, 1, stencil_0_00, "double", OPS_READ),
ops_arg_dat(u1_B0, 1, stencil_0_00, "double", OPS_READ),
ops_arg_dat(U0_B0, 1, stencil_0_00, "double", OPS_WRITE));


int iteration_range_4_block0[] = {0, block0np0, 0, block0np1, -3, block0np2 + 4};
ops_par_loop(opensbliblock00Kernel004, "CRU2", opensbliblock00, 3, iteration_range_4_block0,
ops_arg_dat(u2_B0, 1, stencil_0_00, "double", OPS_READ),
ops_arg_dat(U2_B0, 1, stencil_0_00, "double", OPS_WRITE));


int iteration_range_31_block0[] = {-2, block0np0 + 2, -2, block0np1 + 2, -2, block0np2 + 2};
ops_par_loop(opensbliblock00Kernel031, "CRT_B0", opensbliblock00, 3, iteration_range_31_block0,
ops_arg_dat(rho_B0, 1, stencil_0_00, "double", OPS_READ),
ops_arg_dat(p_B0, 1, stencil_0_00, "double", OPS_READ),
ops_arg_dat(T_B0, 1, stencil_0_00, "double", OPS_WRITE));


int iteration_range_20_block0[] = {-3, block0np0 + 4, -3, block0np1 + 4, -3, block0np2 + 4};
ops_par_loop(opensbliblock00Kernel020, "CRa", opensbliblock00, 3, iteration_range_20_block0,
ops_arg_dat(rho_B0, 1, stencil_0_00, "double", OPS_READ),
ops_arg_dat(p_B0, 1, stencil_0_00, "double", OPS_READ),
ops_arg_dat(a_B0, 1, stencil_0_00, "double", OPS_WRITE));


int iteration_range_30_block0[] = {-2, block0np0 + 2, -2, block0np1 + 2, -2, block0np2 + 2};
ops_par_loop(opensbliblock00Kernel030, "CRmu_B0", opensbliblock00, 3, iteration_range_30_block0,
ops_arg_dat(T_B0, 1, stencil_0_00, "double", OPS_READ),
ops_arg_dat(mu_B0, 1, stencil_0_00, "double", OPS_WRITE));


int iteration_range_0_block0[] = {-1, block0np0 + 1, 0, block0np1, 0, block0np2};
ops_par_loop(opensbliblock00Kernel000, "LLFWeno_reconstruction_0_direction", opensbliblock00, 3, iteration_range_0_block0,
ops_arg_dat(U0_B0, 1, stencil_0_06, "double", OPS_READ),
ops_arg_dat(rhou0_B0, 1, stencil_0_06, "double", OPS_READ),
ops_arg_dat(D01_B0, 1, stencil_0_06, "double", OPS_READ),
ops_arg_dat(rhou2_B0, 1, stencil_0_06, "double", OPS_READ),
ops_arg_dat(rhoE_B0, 1, stencil_0_06, "double", OPS_READ),
ops_arg_dat(D00_B0, 1, stencil_0_06, "double", OPS_READ),
ops_arg_dat(rho_B0, 1, stencil_0_06, "double", OPS_READ),
ops_arg_dat(p_B0, 1, stencil_0_06, "double", OPS_READ),
ops_arg_dat(detJ_B0, 1, stencil_0_06, "double", OPS_READ),
ops_arg_dat(rhou1_B0, 1, stencil_0_06, "double", OPS_READ),
ops_arg_dat(a_B0, 1, stencil_0_06, "double", OPS_READ),
ops_arg_dat(u2_B0, 1, stencil_0_05, "double", OPS_READ),
ops_arg_dat(u0_B0, 1, stencil_0_06, "double", OPS_READ),
ops_arg_dat(u1_B0, 1, stencil_0_06, "double", OPS_READ),
ops_arg_dat(wk3_B0, 1, stencil_0_00, "double", OPS_WRITE),
ops_arg_dat(wk0_B0, 1, stencil_0_00, "double", OPS_WRITE),
ops_arg_dat(wk1_B0, 1, stencil_0_00, "double", OPS_WRITE),
ops_arg_dat(wk4_B0, 1, stencil_0_00, "double", OPS_WRITE),
ops_arg_dat(wk2_B0, 1, stencil_0_00, "double", OPS_WRITE));


int iteration_range_1_block0[] = {0, block0np0, -1, block0np1 + 1, 0, block0np2};
ops_par_loop(opensbliblock00Kernel001, "LLFWeno_reconstruction_1_direction", opensbliblock00, 3, iteration_range_1_block0,
ops_arg_dat(rhou2_B0, 1, stencil_0_08, "double", OPS_READ),
ops_arg_dat(u1_B0, 1, stencil_0_08, "double", OPS_READ),
ops_arg_dat(u2_B0, 1, stencil_0_07, "double", OPS_READ),
ops_arg_dat(rhoE_B0, 1, stencil_0_08, "double", OPS_READ),
ops_arg_dat(rho_B0, 1, stencil_0_08, "double", OPS_READ),
ops_arg_dat(D11_B0, 1, stencil_0_08, "double", OPS_READ),
ops_arg_dat(D10_B0, 1, stencil_0_08, "double", OPS_READ),
ops_arg_dat(detJ_B0, 1, stencil_0_08, "double", OPS_READ),
ops_arg_dat(rhou1_B0, 1, stencil_0_08, "double", OPS_READ),
ops_arg_dat(a_B0, 1, stencil_0_08, "double", OPS_READ),
ops_arg_dat(rhou0_B0, 1, stencil_0_08, "double", OPS_READ),
ops_arg_dat(u0_B0, 1, stencil_0_08, "double", OPS_READ),
ops_arg_dat(U1_B0, 1, stencil_0_08, "double", OPS_READ),
ops_arg_dat(p_B0, 1, stencil_0_08, "double", OPS_READ),
ops_arg_dat(wk8_B0, 1, stencil_0_00, "double", OPS_WRITE),
ops_arg_dat(wk9_B0, 1, stencil_0_00, "double", OPS_WRITE),
ops_arg_dat(wk6_B0, 1, stencil_0_00, "double", OPS_WRITE),
ops_arg_dat(wk7_B0, 1, stencil_0_00, "double", OPS_WRITE),
ops_arg_dat(wk5_B0, 1, stencil_0_00, "double", OPS_WRITE));


int iteration_range_2_block0[] = {0, block0np0, 0, block0np1, -1, block0np2 + 1};
ops_par_loop(opensbliblock00Kernel002, "LLFWeno_reconstruction_2_direction", opensbliblock00, 3, iteration_range_2_block0,
ops_arg_dat(rhou2_B0, 1, stencil_0_10, "double", OPS_READ),
ops_arg_dat(u2_B0, 1, stencil_0_10, "double", OPS_READ),
ops_arg_dat(rhoE_B0, 1, stencil_0_10, "double", OPS_READ),
ops_arg_dat(rho_B0, 1, stencil_0_10, "double", OPS_READ),
ops_arg_dat(p_B0, 1, stencil_0_10, "double", OPS_READ),
ops_arg_dat(detJ_B0, 1, stencil_0_10, "double", OPS_READ),
ops_arg_dat(rhou1_B0, 1, stencil_0_10, "double", OPS_READ),
ops_arg_dat(a_B0, 1, stencil_0_10, "double", OPS_READ),
ops_arg_dat(rhou0_B0, 1, stencil_0_10, "double", OPS_READ),
ops_arg_dat(u0_B0, 1, stencil_0_09, "double", OPS_READ),
ops_arg_dat(u1_B0, 1, stencil_0_09, "double", OPS_READ),
ops_arg_dat(U2_B0, 1, stencil_0_10, "double", OPS_READ),
ops_arg_dat(wk13_B0, 1, stencil_0_00, "double", OPS_WRITE),
ops_arg_dat(wk10_B0, 1, stencil_0_00, "double", OPS_WRITE),
ops_arg_dat(wk12_B0, 1, stencil_0_00, "double", OPS_WRITE),
ops_arg_dat(wk14_B0, 1, stencil_0_00, "double", OPS_WRITE),
ops_arg_dat(wk11_B0, 1, stencil_0_00, "double", OPS_WRITE));


int iteration_range_21_block0[] = {0, block0np0, 0, block0np1, 0, block0np2};
ops_par_loop(opensbliblock00Kernel021, "LLFWeno Residual", opensbliblock00, 3, iteration_range_21_block0,
ops_arg_dat(wk14_B0, 1, stencil_0_11, "double", OPS_READ),
ops_arg_dat(wk11_B0, 1, stencil_0_11, "double", OPS_READ),
ops_arg_dat(wk8_B0, 1, stencil_0_12, "double", OPS_READ),
ops_arg_dat(wk9_B0, 1, stencil_0_12, "double", OPS_READ),
ops_arg_dat(wk6_B0, 1, stencil_0_12, "double", OPS_READ),
ops_arg_dat(wk4_B0, 1, stencil_0_13, "double", OPS_READ),
ops_arg_dat(wk5_B0, 1, stencil_0_12, "double", OPS_READ),
ops_arg_dat(wk3_B0, 1, stencil_0_13, "double", OPS_READ),
ops_arg_dat(wk0_B0, 1, stencil_0_13, "double", OPS_READ),
ops_arg_dat(wk1_B0, 1, stencil_0_13, "double", OPS_READ),
ops_arg_dat(detJ_B0, 1, stencil_0_00, "double", OPS_READ),
ops_arg_dat(wk13_B0, 1, stencil_0_11, "double", OPS_READ),
ops_arg_dat(wk10_B0, 1, stencil_0_11, "double", OPS_READ),
ops_arg_dat(wk12_B0, 1, stencil_0_11, "double", OPS_READ),
ops_arg_dat(wk7_B0, 1, stencil_0_12, "double", OPS_READ),
ops_arg_dat(wk2_B0, 1, stencil_0_13, "double", OPS_READ),
ops_arg_dat(Residual4_B0, 1, stencil_0_00, "double", OPS_WRITE),
ops_arg_dat(Residual3_B0, 1, stencil_0_00, "double", OPS_WRITE),
ops_arg_dat(Residual0_B0, 1, stencil_0_00, "double", OPS_WRITE),
ops_arg_dat(Residual2_B0, 1, stencil_0_00, "double", OPS_WRITE),
ops_arg_dat(Residual1_B0, 1, stencil_0_00, "double", OPS_WRITE));


int iteration_range_40_block0[] = {0, block0np0, 0, block0np1, -2, block0np2 + 2};
ops_par_loop(opensbliblock00Kernel040, "Derivative evaluation CD u2_B0 xi1 ", opensbliblock00, 3, iteration_range_40_block0,
ops_arg_dat(u2_B0, 1, stencil_0_04, "double", OPS_READ),
ops_arg_dat(wk6_B0, 1, stencil_0_00, "double", OPS_WRITE),
ops_arg_idx());


int iteration_range_35_block0[] = {0, block0np0, -2, block0np1 + 2, -2, block0np2 + 2};
ops_par_loop(opensbliblock00Kernel035, "Derivative evaluation CD u0_B0 xi0 ", opensbliblock00, 3, iteration_range_35_block0,
ops_arg_dat(u0_B0, 1, stencil_0_14, "double", OPS_READ),
ops_arg_dat(wk2_B0, 1, stencil_0_00, "double", OPS_WRITE),
ops_arg_idx());


int iteration_range_43_block0[] = {0, block0np0, -2, block0np1 + 2, -2, block0np2 + 2};
ops_par_loop(opensbliblock00Kernel043, "Derivative evaluation CD u2_B0 xi0 ", opensbliblock00, 3, iteration_range_43_block0,
ops_arg_dat(u2_B0, 1, stencil_0_14, "double", OPS_READ),
ops_arg_dat(wk8_B0, 1, stencil_0_00, "double", OPS_WRITE),
ops_arg_idx());


int iteration_range_42_block0[] = {0, block0np0, 0, block0np1, 0, block0np2};
ops_par_loop(opensbliblock00Kernel042, "Derivative evaluation CD u2_B0 xi2 ", opensbliblock00, 3, iteration_range_42_block0,
ops_arg_dat(u2_B0, 1, stencil_0_15, "double", OPS_READ),
ops_arg_dat(wk7_B0, 1, stencil_0_00, "double", OPS_WRITE));


int iteration_range_34_block0[] = {0, block0np0, 0, block0np1, 0, block0np2};
ops_par_loop(opensbliblock00Kernel034, "Derivative evaluation CD u0_B0 xi2 ", opensbliblock00, 3, iteration_range_34_block0,
ops_arg_dat(u0_B0, 1, stencil_0_15, "double", OPS_READ),
ops_arg_dat(wk1_B0, 1, stencil_0_00, "double", OPS_WRITE));


int iteration_range_44_block0[] = {0, block0np0, 0, block0np1, 0, block0np2};
ops_par_loop(opensbliblock00Kernel044, "Derivative evaluation CD T_B0 xi1 ", opensbliblock00, 3, iteration_range_44_block0,
ops_arg_dat(T_B0, 1, stencil_0_04, "double", OPS_READ),
ops_arg_dat(wk9_B0, 1, stencil_0_00, "double", OPS_WRITE),
ops_arg_idx());


int iteration_range_39_block0[] = {0, block0np0, -2, block0np1 + 2, -2, block0np2 + 2};
ops_par_loop(opensbliblock00Kernel039, "Derivative evaluation CD u1_B0 xi0 ", opensbliblock00, 3, iteration_range_39_block0,
ops_arg_dat(u1_B0, 1, stencil_0_14, "double", OPS_READ),
ops_arg_dat(wk5_B0, 1, stencil_0_00, "double", OPS_WRITE),
ops_arg_idx());


int iteration_range_38_block0[] = {0, block0np0, 0, block0np1, 0, block0np2};
ops_par_loop(opensbliblock00Kernel038, "Derivative evaluation CD u1_B0 xi2 ", opensbliblock00, 3, iteration_range_38_block0,
ops_arg_dat(u1_B0, 1, stencil_0_15, "double", OPS_READ),
ops_arg_dat(wk4_B0, 1, stencil_0_00, "double", OPS_WRITE));


int iteration_range_46_block0[] = {0, block0np0, 0, block0np1, 0, block0np2};
ops_par_loop(opensbliblock00Kernel046, "Derivative evaluation CD T_B0 xi2 ", opensbliblock00, 3, iteration_range_46_block0,
ops_arg_dat(T_B0, 1, stencil_0_15, "double", OPS_READ),
ops_arg_dat(wk10_B0, 1, stencil_0_00, "double", OPS_WRITE));


int iteration_range_47_block0[] = {0, block0np0, -2, block0np1 + 2, 0, block0np2};
ops_par_loop(opensbliblock00Kernel047, "Derivative evaluation CD T_B0 xi0 ", opensbliblock00, 3, iteration_range_47_block0,
ops_arg_dat(T_B0, 1, stencil_0_14, "double", OPS_READ),
ops_arg_dat(wk11_B0, 1, stencil_0_00, "double", OPS_WRITE),
ops_arg_idx());


int iteration_range_36_block0[] = {0, block0np0, 0, block0np1, -2, block0np2 + 2};
ops_par_loop(opensbliblock00Kernel036, "Derivative evaluation CD u1_B0 xi1 ", opensbliblock00, 3, iteration_range_36_block0,
ops_arg_dat(u1_B0, 1, stencil_0_04, "double", OPS_READ),
ops_arg_dat(wk3_B0, 1, stencil_0_00, "double", OPS_WRITE),
ops_arg_idx());


int iteration_range_32_block0[] = {0, block0np0, 0, block0np1, -2, block0np2 + 2};
ops_par_loop(opensbliblock00Kernel032, "Derivative evaluation CD u0_B0 xi1 ", opensbliblock00, 3, iteration_range_32_block0,
ops_arg_dat(u0_B0, 1, stencil_0_04, "double", OPS_READ),
ops_arg_dat(wk0_B0, 1, stencil_0_00, "double", OPS_WRITE),
ops_arg_idx());


int iteration_range_56_block0[] = {0, block0np0, 0, block0np1, 0, block0np2};
ops_par_loop(opensbliblock00Kernel056, "Viscous terms", opensbliblock00, 3, iteration_range_56_block0,
ops_arg_dat(SD100_B0, 1, stencil_0_00, "double", OPS_READ),
ops_arg_dat(SD000_B0, 1, stencil_0_00, "double", OPS_READ),
ops_arg_dat(wk8_B0, 1, stencil_0_01, "double", OPS_READ),
ops_arg_dat(wk5_B0, 1, stencil_0_01, "double", OPS_READ),
ops_arg_dat(u2_B0, 1, stencil_0_02, "double", OPS_READ),
ops_arg_dat(SD110_B0, 1, stencil_0_00, "double", OPS_READ),
ops_arg_dat(u0_B0, 1, stencil_0_02, "double", OPS_READ),
ops_arg_dat(T_B0, 1, stencil_0_02, "double", OPS_READ),
ops_arg_dat(D01_B0, 1, stencil_0_00, "double", OPS_READ),
ops_arg_dat(mu_B0, 1, stencil_0_02, "double", OPS_READ),
ops_arg_dat(SD010_B0, 1, stencil_0_00, "double", OPS_READ),
ops_arg_dat(wk1_B0, 1, stencil_0_00, "double", OPS_READ),
ops_arg_dat(SD001_B0, 1, stencil_0_00, "double", OPS_READ),
ops_arg_dat(D11_B0, 1, stencil_0_00, "double", OPS_READ),
ops_arg_dat(D10_B0, 1, stencil_0_00, "double", OPS_READ),
ops_arg_dat(SD111_B0, 1, stencil_0_00, "double", OPS_READ),
ops_arg_dat(SD011_B0, 1, stencil_0_00, "double", OPS_READ),
ops_arg_dat(u1_B0, 1, stencil_0_02, "double", OPS_READ),
ops_arg_dat(SD101_B0, 1, stencil_0_00, "double", OPS_READ),
ops_arg_dat(wk3_B0, 1, stencil_0_03, "double", OPS_READ),
ops_arg_dat(wk0_B0, 1, stencil_0_03, "double", OPS_READ),
ops_arg_dat(D00_B0, 1, stencil_0_00, "double", OPS_READ),
ops_arg_dat(wk10_B0, 1, stencil_0_00, "double", OPS_READ),
ops_arg_dat(wk11_B0, 1, stencil_0_04, "double", OPS_READ),
ops_arg_dat(wk6_B0, 1, stencil_0_03, "double", OPS_READ),
ops_arg_dat(wk4_B0, 1, stencil_0_00, "double", OPS_READ),
ops_arg_dat(wk7_B0, 1, stencil_0_00, "double", OPS_READ),
ops_arg_dat(wk9_B0, 1, stencil_0_00, "double", OPS_READ),
ops_arg_dat(wk2_B0, 1, stencil_0_01, "double", OPS_READ),
ops_arg_dat(Residual4_B0, 1, stencil_0_00, "double", OPS_RW),
ops_arg_dat(Residual3_B0, 1, stencil_0_00, "double", OPS_RW),
ops_arg_dat(Residual2_B0, 1, stencil_0_00, "double", OPS_RW),
ops_arg_dat(Residual1_B0, 1, stencil_0_00, "double", OPS_RW),
ops_arg_idx());


int iteration_range_82_block0[] = {0, block0np0, 0, block0np1, 0, block0np2};
ops_par_loop(opensbliblock00Kernel082, "Temporal solution advancement", opensbliblock00, 3, iteration_range_82_block0,
ops_arg_dat(Residual4_B0, 1, stencil_0_00, "double", OPS_READ),
ops_arg_dat(Residual3_B0, 1, stencil_0_00, "double", OPS_READ),
ops_arg_dat(Residual0_B0, 1, stencil_0_00, "double", OPS_READ),
ops_arg_dat(Residual2_B0, 1, stencil_0_00, "double", OPS_READ),
ops_arg_dat(Residual1_B0, 1, stencil_0_00, "double", OPS_READ),
ops_arg_dat(tempRK_rhou0_B0, 1, stencil_0_00, "double", OPS_RW),
ops_arg_dat(tempRK_rhou2_B0, 1, stencil_0_00, "double", OPS_RW),
ops_arg_dat(tempRK_rhoE_B0, 1, stencil_0_00, "double", OPS_RW),
ops_arg_dat(rhoE_B0, 1, stencil_0_00, "double", OPS_RW),
ops_arg_dat(rho_B0, 1, stencil_0_00, "double", OPS_RW),
ops_arg_dat(tempRK_rhou1_B0, 1, stencil_0_00, "double", OPS_RW),
ops_arg_dat(rhou1_B0, 1, stencil_0_00, "double", OPS_RW),
ops_arg_dat(tempRK_rho_B0, 1, stencil_0_00, "double", OPS_RW),
ops_arg_dat(rhou0_B0, 1, stencil_0_00, "double", OPS_RW),
ops_arg_dat(rhou2_B0, 1, stencil_0_00, "double", OPS_RW),
ops_arg_gbl(&rkA[stage], 1, "double", OPS_READ),
ops_arg_gbl(&rkB[stage], 1, "double", OPS_READ));


int iteration_range_61_block0[] = {0, 1, -3, block0np1 + 4, -3, block0np2 + 4};
ops_par_loop(opensbliblock00Kernel061, "InletTransfer boundary dir0 side0", opensbliblock00, 3, iteration_range_61_block0,
ops_arg_dat(rhou1_B0, 1, stencil_0_13, "double", OPS_RW),
ops_arg_dat(rhoE_B0, 1, stencil_0_13, "double", OPS_RW),
ops_arg_dat(rhou0_B0, 1, stencil_0_13, "double", OPS_RW),
ops_arg_dat(rho_B0, 1, stencil_0_13, "double", OPS_RW),
ops_arg_dat(rhou2_B0, 1, stencil_0_13, "double", OPS_RW));


int iteration_range_62_block0[] = {block0np0 - 1, block0np0, -3, block0np1 + 4, -3, block0np2 + 4};
ops_par_loop(opensbliblock00Kernel062, "Extrapolation boundary dir0 side1", opensbliblock00, 3, iteration_range_62_block0,
ops_arg_dat(rhoE_B0, 1, stencil_0_16, "double", OPS_RW),
ops_arg_dat(rhou1_B0, 1, stencil_0_16, "double", OPS_RW),
ops_arg_dat(rhou0_B0, 1, stencil_0_16, "double", OPS_RW),
ops_arg_dat(rho_B0, 1, stencil_0_16, "double", OPS_RW),
ops_arg_dat(rhou2_B0, 1, stencil_0_16, "double", OPS_RW));


int iteration_range_63_block0[] = {-3, block0np0 + 4, 0, 1, -3, block0np2 + 4};
ops_par_loop(opensbliblock00Kernel063, "ForcingStripWall boundary dir1 side0", opensbliblock00, 3, iteration_range_63_block0,
ops_arg_dat(x2_B0, 1, stencil_0_00, "double", OPS_READ),
ops_arg_dat(x0_B0, 1, stencil_0_00, "double", OPS_READ),
ops_arg_dat(rhoE_B0, 1, stencil_0_17, "double", OPS_RW),
ops_arg_dat(rhou1_B0, 1, stencil_0_18, "double", OPS_RW),
ops_arg_dat(rhou0_B0, 1, stencil_0_18, "double", OPS_RW),
ops_arg_dat(rho_B0, 1, stencil_0_18, "double", OPS_RW),
ops_arg_dat(rhou2_B0, 1, stencil_0_18, "double", OPS_RW),
ops_arg_gbl(&iter, 1, "int", OPS_READ));


int iteration_range_64_block0[] = {-3, block0np0 + 4, block0np1 - 1, block0np1, -3, block0np2 + 4};
ops_par_loop(opensbliblock00Kernel064, "ZeroGradientOutlet boundary dir1 side1", opensbliblock00, 3, iteration_range_64_block0,
ops_arg_dat(rhoE_B0, 1, stencil_0_19, "double", OPS_RW),
ops_arg_dat(rhou1_B0, 1, stencil_0_19, "double", OPS_RW),
ops_arg_dat(rhou0_B0, 1, stencil_0_19, "double", OPS_RW),
ops_arg_dat(rho_B0, 1, stencil_0_19, "double", OPS_RW),
ops_arg_dat(rhou2_B0, 1, stencil_0_19, "double", OPS_RW));


ops_halo_transfer(exchange65_block0);
ops_halo_transfer(exchange66_block0);
}
if (fmod(iter + 1,2500) == 0){
char name0[80];
sprintf(name0, "opensbli_output_%06d.h5", iter + 1);
ops_fetch_block_hdf5_file(opensbliblock00, name0);
ops_fetch_dat_hdf5_file(rho_B0, name0);
ops_fetch_dat_hdf5_file(rhou0_B0, name0);
ops_fetch_dat_hdf5_file(rhou1_B0, name0);
ops_fetch_dat_hdf5_file(rhou2_B0, name0);
ops_fetch_dat_hdf5_file(rhoE_B0, name0);
ops_fetch_dat_hdf5_file(x0_B0, name0);
ops_fetch_dat_hdf5_file(x1_B0, name0);
ops_fetch_dat_hdf5_file(x2_B0, name0);
ops_fetch_dat_hdf5_file(D11_B0, name0);
ops_fetch_dat_hdf5_file(p_B0, name0);
}
}
double cpu_end0, elapsed_end0;
ops_timers(&cpu_end0, &elapsed_end0);
ops_printf("\nTimings are:\n");
ops_printf("-----------------------------------------\n");
ops_printf("Total Wall time %lf\n",elapsed_end0-elapsed_start0);
char name0[80];
sprintf(name0, "opensbli_output.h5");
ops_fetch_block_hdf5_file(opensbliblock00, name0);
ops_fetch_dat_hdf5_file(rho_B0, name0);
ops_fetch_dat_hdf5_file(rhou0_B0, name0);
ops_fetch_dat_hdf5_file(rhou1_B0, name0);
ops_fetch_dat_hdf5_file(rhou2_B0, name0);
ops_fetch_dat_hdf5_file(rhoE_B0, name0);
ops_fetch_dat_hdf5_file(x0_B0, name0);
ops_fetch_dat_hdf5_file(x1_B0, name0);
ops_fetch_dat_hdf5_file(x2_B0, name0);
ops_fetch_dat_hdf5_file(D11_B0, name0);
ops_fetch_dat_hdf5_file(p_B0, name0);
ops_exit();
//Main program end 
}