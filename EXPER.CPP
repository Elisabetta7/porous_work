#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <dos.h>
//#include <values.h>
#include "defglob.h"
#include "matr.h"
#include "inpar.h"
#include "reoeq1.h"
#include "taskporous.h"
#include "diaproc.h"
#define MAX_STEP_DECREASE 0 // было 8
//#define MALBUF 16 - in diaproc.h
//extern unsigned _stklen = 30000; - defined in Module definition file "exper.def"

Real Excitation_displacement(TASK tsk ,Real tim, int k_harmonic, Real *ampl, Real *angular_frequency)
{
  Real res=0.0;

  if(tsk != vibr_trig_polynom) return 0.0;//no excitation
  else {//if excitation is specified by a trigonometric polynom
	  res=0.0;
	  for(int i=0;i<k_harmonic;i++){
		res+=ampl[i]*sin(angular_frequency[i]*tim);
	  }
  }
  return res;
}

Real Excitation_velocity(TASK tsk ,Real tim, int k_harmonic, Real *ampl, Real *angular_frequency)
{
  Real res=0.0;
  if(tsk != vibr_trig_polynom) return 0.0;//no excitation
  else {//if excitation is specified by a trigonometric polynom
	  res=0.0;
	  for(int i=0;i<k_harmonic;i++){
		Real f=angular_frequency[i];
		res+=ampl[i]*f*cos(f*tim);
	  }
  }
  return res;
}

Real Excitation_acceleration(TASK tsk ,Real tim, int k_harmonic, Real *ampl, Real *angular_frequency)
{
  Real res=0.0;
  if(tsk != vibr_trig_polynom) return 0.0;//no excitation
  else {//if excitation is specified by a trigonometric polynom
	  res=0.0;
	  for(int i=0;i<k_harmonic;i++){
		Real f=angular_frequency[i];
		res-=ampl[i]*f*f*sin(f*tim);
	  }
  }
  return res;
}

void print_step_header(FILE* fcyc){
  fprintf(fcyc,"\n###1time/s  2tem/K      3e33/pct   4s33/MPa       5g32/pct   6s32/MPa   7faza      8Ro/MPa   9iz   10bf_max 11kdmg 12Bmax_pl 13Bmean_pl 14Fmax_pl 15Fmean_pl");
}

void print_step(FILE* fttesf ,int print_step_flag,  int na_ekran, int istep, double tim, double tem, 
                 Matr33& sig, double force, Matr33& eps, InternalParameters& XX){
   if(print_step_flag){
       Real faza = XX.Phase;
       int k_dmg=XX.kzer_dmg;
       Real Ro = XX.RoGrainNorma*1e-6;
       Real s33=sig(2,2), s32=sig(2,1);
       Real e33=eps(2,2)*1e2, g32=eps(2,1)*1e2;
	   Real Bmax_pl=XX.Bmax_pl_dmg;
	   Real Bmean_pl=XX.Bmean_pl_dmg;
	   Real Fmax_pl=XX.Fmax_pl_dmg;
	   Real Fmean_pl=XX.Fmean_pl_dmg;
       int iz=XX.izer_bf_max; //numb. of grain with maximum value of damage criterion fabs(b)*f
       Real bf=XX.bf_dmg[iz]; //maximum abs(b)*f
     //printing:
       //printf("\r%4i %7.4lf %3.0lfK F=%2.0lf s/e%5.1lf %5.2lf t/g%5.1lf %5.2lf %3i %3i %4.2lf r=%2.0lf"
       //       ,istep,tim   ,tem    ,faza*1e2    ,s33   ,e33       ,s32   ,g32 ,iz,k_dmg,bf   ,Ro);
       if(na_ekran==1){
		printf("\r%6.3lf %3.0lfK %3.0lfK F=%2.0lf s/e%5.1lf %5.2lf t/g%5.1lf %5.2lf %3i %3i"         
			,istep,tim   ,tem	  ,faza*1e2  ,s33   ,e33     ,s32   ,g32    ,iz ,k_dmg );
       }
       fprintf(fttesf,"\n%10.3lg %5.1lf %5.1lf % 11.6lf % 8.4lf % 11.6lf % 8.4lf  %7.5lf   %5i %12.6lf %5i %le %le %le %le"
			 ,tim     ,tem 	,e33     ,s33     ,g32     ,s32    ,faza   ,Ro     ,iz ,bf    ,k_dmg  ,Bmax_pl ,Bmean_pl ,Fmax_pl ,Fmean_pl);
   }
 }

//void print_step_header2elems(FILE* f){
//  fprintf(f,"\n###1time/s  2tem/K    3g32/pct  4s32/MPa      5s32/MPa   6faza      7Ro/MPa   8iz   9bf_max 10kdmg 11Bmax_pl 12Bmean_pl 13Fmax_pl 14Fmean_pl");
//  fprintf(f,"\n");
//}
//void print_step2elems(FILE* fttesf ,int print_step_flag,  int na_ekran, int istep, double tim, double tem,
//                 Matr33& sig, InternalParameters& XX1, InternalParameters& XX2){
//   if(print_step_flag){
//       Real faza = XX2.Phase;
//       int k_dmg=XX2.kzer_dmg;
//       Real Ro = XX2.RoGrainNorma*1e-6;
//       Real s32=sig(2,1);
//	   Real g32=XX1.TotalStrain.C(2,1);
//	   Real Bmax_pl=XX2.Bmax_pl_dmg;
//	   Real Bmean_pl=XX2.Bmean_pl_dmg;
//	   Real Fmax_pl=XX2.Fmax_pl_dmg;
//	   Real Fmean_pl=XX2.Fmean_pl_dmg;
//       int iz=XX2.izer_bf_max; //numb. of grain with maximum value of damage criterion fabs(b)*f
//       Real bf=XX2.bf_dmg[iz]; //maximum abs(b)*f
//     //printing:
//       //printf("\r%4i %7.4lf %3.0lfK F=%2.0lf s/e=%5.1lf %5.2lf t/g=%5.1lf %5.2lf %3i %3i %4.2lf r=%2.0lf"
//       //       ,istep,tim   ,tem    ,faza*1e2    ,s33   ,e33       ,s32   ,g32 ,iz,k_dmg,bf  ,Ro);
//       if(na_ekran==1){
//		printf("\r%4i %7.4lf %3.0lfK F=%2.0lf t/g=1%5.1lf %5.2lf  %3i %3i"                                     //....добавлено
//				,istep,tim   ,tem ,faza*1e2    ,s32 ,g32     ,iz ,k_dmg  );
//       }
//       fprintf(fttesf,"\n%10.6lg  %6.1lf  % 12.6lf % 10.6lf %7.5lf  %12.6lf %5i %12.6lf   %5i   %le      %le     %le      %le"
//		           ,tim     ,tem      ,g32    ,s32    ,faza     ,Ro   ,iz   ,bf   ,k_dmg,Bmax_pl,Bmean_pl,Fmax_pl,Fmean_pl);
//   } 
// }///////// колонки

void print_step_header_2slotted_elems(FILE* f){ 
  fprintf(f,"\n###1time/s   2tem/K     3e33_1/%%  4s33_1/MPa 5e33_2/%% 6s33_2/MPa 7e33_3/%% 8s33_3/MPa 9e33_4/%%  10s33_4/MPa  11lenght_1/m 12s22_1/MPa  13e22_1/%%  14faza_1   15faza_2  16faza_3	17faza_4  18faza_5 19force1/N  20E_macro/%%  21Sigma_macro/MPa  22e33_5/%%  23e_6/%%");
  fprintf(f,"\n");
}
//print_step_2slotted_elems(fttesf ,1 ,0, 0, tim, tem, sig, XXnew1, XXnew2);
void print_step_2slotted_elems(FILE* fttesf ,int print_step_flag,  int na_ekran, int istep, double tim, double tem,
	double force_z, int Kuz, Circular_beam* Cb, double displ_1, double displ_2, double displ_3, double displ_4, double e_macro, double stress_macro, double force_new_eff){
   if(print_step_flag){
		int iuz_max=Kuz-1;
		Real faza_1 = Cb[0].XXnew.Phase;
		Real faza_2 = Cb[1].XXnew.Phase;    
		Real faza_3 = Cb[2].XXnew.Phase;
		Real faza_4 = Cb[3].XXnew.Phase;
		Real faza_5 = Cb[4].XXnew.Phase;
		//Real faza2 = Cb2[iuz_max].XXnew.Phase;
		int k_dmg=Cb[iuz_max].XXnew.kzer_dmg;
		Real Ro = Cb[iuz_max].XXnew.RoGrainNorma*1e-6;
		//Real s33=Cb[iuz_max].K_sig * force;    
		Real s33_1 =/*Cb[0].K_sig * force_z*/Cb[0].sigma33*1e-6;
		Real s33_2 =/*Cb[1].K_sig * force_z*/Cb[1].sigma33*1e-6;
		Real s33_3 = /*Cb[2].K_sig*force_z*/Cb[2].sigma33*1e-6;
		Real s33_4 = /*Cb[3].K_sig*force_z*/Cb[3].sigma33*1e-6;
		//real s33_2=Cb2[iuz_max].K_sig * force2*1e-6;
		Real s22_1 = Cb[0].sig_to_control*1e-6;
		Real e22_1 = Cb[0].XXnew.TotalStrain.C(1, 1)*1e2;
		Real e33_1 = Cb[0].XXnew.TotalStrain.C(2,2)*1e2;
		Real e33_2 = Cb[1].XXnew.TotalStrain.C(2,2)*1e2;
		Real e33_3 = Cb[2].XXnew.TotalStrain.C(2, 2)*1e2;
		Real e33_4 = Cb[3].XXnew.TotalStrain.C(2, 2)*1e2;
		Real e33_5 = Cb[4].XXnew.TotalStrain.C(2, 2)*1e2;
		Real e33_6 = Cb[5].XXnew.TotalStrain.C(2, 2)*1e2;
		Real length_1 = Cb[0].L_macro_new;
		//real e33_2=Cb2[iuz_max].XXnew.TotalStrain.C(2,2)*1e2;
	   //Real Bmax_pl=XX2.Bmax_pl_dmg;
	   //Real Bmean_pl=XX2.Bmean_pl_dmg;
	   //Real Fmax_pl=XX2.Fmax_pl_dmg;
	   //Real Fmean_pl=XX2.Fmean_pl_dmg;
       int iz=Cb[iuz_max].XXnew.izer_bf_max; //numb. of grain with maximum value of damage criterion fabs(b)*f
       Real bf=Cb[iuz_max].XXnew.bf_dmg[iz]; //maximum abs(b)*f
	   //printing:
       //printf("\r%4i %7.4lf %3.0lfK   F=%2.0lf   s/e%5.1lf %5.2lf     t/g%5.1lf  %5.2lf   %3i  %3i    %4.2lf  r=%2.0lf"
       //       ,istep ,tim    ,tem    ,faza*1e2    ,s33      ,e33       ,s32       ,g32    ,iz ,k_dmg,   bf     ,Ro);
       if(na_ekran==1 && istep%20==0){
		printf("\r%5i   time=%6.4lg   T=%5.8lg  F_1=%5.8lg  F_2=%5.8lg  e33_1=%6.8lg  e33_2=%6.8lg   s33_1=%6.8lg  s33_2=%6.8lg"
				,istep   ,tim           ,tem	 ,faza_1     ,faza_2      ,e33_1        ,e33_2        ,s33_1          ,s33_2);
       }
       fprintf(fttesf,"\n%10.6lg  %7.8lg  %7.8lg   %7.8lg   %7.8lg   %7.8lg   %7.8lg  %7.8lg   %7.8lg  %7.8lg    %7.8lg     %7.8lg     %7.8lg   %7.8lg     %7.8lg     %7.8lg   %7.8lg      %7.8lg      %7.8lg          %7.8lg             %7.8lg      %7.8lg  %7.8lg"
		                   , tim,  tem   , e33_1   , s33_1  , e33_2 , s33_2 , e33_3  , s33_3   , e33_4, s33_4  , length_1,  s22_1  ,e22_1 ,  faza_1   , faza_2,  faza_3, faza_4, faza_5,	force_z, e_macro * (-100), stress_macro*1e-6,  e33_5,  e33_6);//////колонки
	   //we print forces applied to the shaft, they are opposite to forces applied to RSS (round slotted spring)
   } 
 }


void print_cycle_header(FILE* fcyc){
  fprintf(fcyc,"\n######1cycle#   2time/s    3tem/K   4e33/pct   5s33/MPa     6g32/pct   7s32/MPa   8faza        9Ro/MPa     10iz  11bf_max   12kdmg 13Bmax_pl 14Bmean_pl 15Fmax_pl 16Fmean_pl");
}

void print_etap(FILE* fcyc, int sufficient_condition, CYCLE c, double tim, double tem,
		 Matr33& sig, Matr33& eps, InternalParameters& XX){
   long j=c.i+1L; //current iteration of the cycle
   if(j<=10 || j<=100L&& j%10L==0L || j<=1000L&& j%100L==0L || j<=10000L&& j%1000L==0L || j<=100000L&& j%10000L==0L
	    || j<=1000000L&& j%100000L==0L || j<=10000000L&& j%1000000L==0L || j<=100000000L&& j%10000000L==0L
	    || j<=1000000000L&& j%100000000L==0L
            || sufficient_condition){
      Real faza = XX.Phase;
      int k_dmg=XX.kzer_dmg;
      Real Ro = XX.RoGrainNorma*1e-6;
      Real s33=sig(2,2), s32=sig(2,1);
      Real e33=eps(2,2)*1e2, g32=eps(2,1)*2e2;
      int iz=XX.izer_bf_max; //numb. of grain with maximum abs(b)*f
      Real bf=XX.bf_dmg[iz]; //maximum abs(b)*f
	   Real Bmax_pl=XX.Bmax_pl_dmg;
	   Real Bmean_pl=XX.Bmean_pl_dmg;
	   Real Fmax_pl=XX.Fmax_pl_dmg;
	   Real Fmean_pl=XX.Fmean_pl_dmg;
    //printing:
      fprintf(fcyc,"\n %10li %6.3lg %6.3lg %6.1lf % 12.6lf % 10.6lf % 12.6lf % 10.6lf  %7.5lf   %5i %12.6lf %5i %le %le %le %le"
		       ,j    ,tim    ,tem  ,e33    ,s33     ,g32     ,s32      ,faza    ,Ro     ,iz   ,bf  ,k_dmg  ,Bmax_pl ,Bmean_pl ,Fmax_pl ,Fmean_pl);
    }
 }

int main(int argc, char **argv){
// struct time t;
 int retcode=0;
 int istep;//for counting steps of the experiment
 int na_ekran; //to print or not each step on the screen
 int print_step_flag; // this flag indicates printing to file ttesgtf.dat
 int pause_factor=0;
 pause(0);

int operation_flag=0; //this flag is used when vibration or other impact is applied
int operation_flag_prev; //this flag is used when vibration or other impact is applied
int load_flag;
int preload_flag; //preload_flag=1 means independent loading and heating of RSS
int add_magnetic_force_flag;//this flag is used when we introduce magnetic force
int vibr_trig_polynom_flag;//this flag is used when excitation ia a trigonometric polynom
int idle_flag; // this flag is used when there is no force excitation
				// we may wait or change temperatures of the 1-st or of the 2nd RSS or both

 int control;   //0 - free vibrations, 1- controlled vibrations
 char impulse_shape[MALBUF]; //it may == "1exp", "2exp", "3exp", "calculated"
 char wait_mode[MALBUF]; //it may == "seconds", "fractions"
 char tem_peaks_mode[MALBUF];//it may == "equal" or name-of-file for example: "tem_peak.d"
 //char excitation[MALBUF]; 
 int heatnum;//number of heating impulses
 int heating;//this flag==1 when heating is on (heating may not be on the first step)
// int k_elements; //k_elements may be 1 or 2 (2 - for vibration isolation by two elements). Default: k_elements=1.
 const int Kuz_max=10; //maximum number of SMA points (representative volumes)
 int Kuz, iuz_max;
 Real dtem_heat;//used only for controlled vibrations
 Real *tem_peaks=NULL;//for array of peak temperatures 
 int *macro_failure_pointer;

//variables needed for vibr32ation and excitation
 int k_harmonic;//number of terms (harmonics) in expansion of displacement or acceleration
// accel = ampl_1*cos((2*pi*freq_1)*tim) +...+ ampl_k*cos((2*pi*freq_k)*tim)
 Real *ampl=NULL   // amplitude of i-th harmonic
      ,*angular_frequency=NULL // frequency of i-th harmonic
      ,ground_velocity //velocity on the ground
      ,hground_velocity //increm. of velocity on the ground
      ,ground_displ   //displacement on the ground (excitation)
      ,hground_displ   //increm. of ground_displ
      ,mass_displ  //displacement of the mass
      ,mass_velocity
	  ,mass_accel
      ,hmass_displ  //increm. of the displacemenbt of the mass
      ,hmass_velocity
	  ,force
	  ,hforce
	  //,force2
	  //,hforce2
	  ,hmf_factor_preset
	;
	  //Synonims:
 //Real &displacement = mass_displ ,&velocity = mass_velocity;
 //Real &excit_displ = ground_displ;

 int VibrContinue;
 int /*semiperiod=0,*/iheat=0;
 Real amplit=MAXDOUBLE, amplit_prev=MAXDOUBLE, decrem=MAXDOUBLE;
 Real tim_rev=MAXDOUBLE, tim_rev_prev=MAXDOUBLE, tim_period=MAXDOUBLE;
 Real tim_wait_s=MAXDOUBLE; //wait time in seconds
 Real eps_peak=MAXDOUBLE, eps_peak_prev=MAXDOUBLE;

 Real siga_us33=MAXDOUBLE, epsa_us33=MAXDOUBLE, freq_us33=MAXDOUBLE;
 Real s_incr_us, e_incr_us;
 Real /*semiperiod_us33,*/ omega_us
     ,tim_etap=MAXDOUBLE
     ,s33_us ,hs33_us, e33_us, he33_us; // s33 and e33 due to ultrasound 

 double f_step,
      tim,          tim_ke,     tim_prev, 
      tim_wait,     tim_heat,   tim_cool,
      htim,         htim_preset,  
	  tim_print,    htim_print,
      
	  tem,         tem_z,     htem,
	  tem_fon,    tem_peak,   tem_ke,

      cool_exp_A,   cool_exp_AM,   cool_exp_M,
      fluence,      hfluence,      fluence_z,

      faza,
      zhest32,      zhest33,
      press,        hpress;

 //Real &tem = tem1, &htem =htem1;

 Real tim_incr=10.0, tem_incr= 1.0,  s_incr=1.0,  e_incr=1.0e-4, faz_incr=0.02; //max allowable increments
 Real htim_tem, htim_s, htim_e;//max time increments imposed by tem_incr, s_incr and e_incr

 FILE *fwzd=NULL, *fttesf=NULL, *fdec=NULL;
 FILE *fvibr; //file for registering vibrograms
 //FILE *fmagforce;//file for registering the process of switching on the magnetic force
 //FILE *fmag_u;
 FILE *ftem_peaks=NULL, *faccel=NULL,  *filrd=NULL; //for any input file

 //Declaration of FILE* and other variables for etap_ending_dat_files
 int k_cycle_files;
 FILE* fcycle[5];
 int ifil;
 int etap_print_numbers[5];
 //End of Declaration of FILE* and other variables for etap_ending_dat_files

 //Begin Declarations for numbers of cycle iterations, whose all steps are to be printed
 int k_icyc_print;
 int icyc;
 long printed_cycles[10];
 //End of Declarations for numbers of cycle iterations, whose all steps are to be printed
 
 Matr33 sig ,hsig  ,sig_prev
	 , eps, heps, deps, eps_prev, heps_prev, sig_elastic, hsig_elastic;
 //Matr33 &sig=sig1, &hsig=hsig1, &eps=eps1;

//---------------------------------------------------------------------------
 printf("\nExperiment started");


	//open files (1) for registering vibrograms and (2) magnetic forces: 
	fvibr = fopen("vibr.dat","w");
	fprintf(fvibr  ,"#tim_etap   Tem      mass_displ    mass_accel    ground_displ    ground_accel   force");
	fprintf(fvibr,"\n#  s         K          mm           m/s**2          mm             m/s**2        N");
	
	//fmagforce = fopen("add_mag_force.dat","w");
	//fprintf(fmagforce  ,"#tim_etap  Tem1  Tem2  mag_force_factor mass_displ   mass_accel");
	//fprintf(fmagforce,"\n#  s        K     K                        mm          m/s**2  ");


 //if argc==3 then last argv is argv[2]
 fttesf = argc>=3 ? fopen(argv[2],"w") : fopen("ttesgtf.dat","w");
 ////что-то исправить в колонках
 //fdec = fopen("decrem.dat","w"); //file for registering decrements

//printf("\n30:");
//pause();

MaterialConstants MC; 
MC.Init();
Real q0_C_Ro = MC.Q0/MC.SpecificHeatM/MC.DensityM;

printf("\n80: MC.Init passed successfully");
//pause();

//---------------- Assign initial values ---------------------------
sig = 0.0; hsig = 0.0; sig_elastic = 0.0; hsig_elastic = 0.0; eps = 0.0; heps = 0.0; deps = 0.0; faza = 0.0; tim = 0.0;// tem=MAXDOUBLE; //tem must be reread
 press=0.0; fluence=0.0;
 force=0.0;
 //if argc==2 then last argv is argv[1]
 fwzd = argc>=2 ? fopen(argv[1],"r") : fopen("wzdi.wzd","r");
 if(fwzd==0) { fprintf(stderr,"\nERROR: wzd-file absent\n"); pause(); return(1); }

 DATA_SECTION ini_section;
 //looking for section [ini_state]
 if( ffindsection(fwzd,"ini_state",ini_section) ){
   fprintf(stderr,"\n\nERROR: Cannot find section [ini_state] in wzd-file");
   fprintf(stderr,"\n Program is terminated.");
   pause();
   exit(1);
 }

	//printf("\n90:");
     //pause(pause_factor);

 fkeyread_double_in_section(fwzd, ini_section, "tim", &tim);
 if(fkeyread_double_in_section(fwzd, ini_section, "tem", &tem)){
   fprintf(stderr,"\n\nERROR: Cannot find key tem= in [ini_state] section of wzd-file");
   fprintf(stderr,"\n Program is terminated.");
   pause();
   exit(1);
 }

  /*if(fkeyread_double_in_section(fwzd, ini_section, "tem2", &tem2)){
   fprintf(stderr,"\n\nERROR: Cannot find key tem2= in [ini_state] section of wzd-file");
   fprintf(stderr,"\n Program is terminated.");
   pause();
   exit(1);
 }*/
     //printf("\n94:");
     //pause(pause_factor);

 if(tem<MC.Ms){
   fprintf(stderr,"\n\nERROR: in [ini_state] temperature must be >= Ms = %lg ", MC.Ms);
   fprintf(stderr,"\n Program is terminated.");
   pause();
   exit(2);
 }

 
  /*if(tem2<MC.Ms){
   fprintf(stderr,"\n\nERROR: in [ini_state] temperature2 must be >= Ms = %lg ", MC.Ms);
   fprintf(stderr,"\n Program is terminated.");
   pause();
   exit(2);
 }*/

// Read allowable icrements of temperature, stress and strain.
// If there are no keys for them in wzd-file, assign default values. 
 if( fkeyread_double_in_section(fwzd, ini_section, "tim_incr", &tim_incr) ) tim_incr=10.0;
 if( fkeyread_double_in_section(fwzd, ini_section, "tem_incr", &tem_incr) ) tem_incr=1.0;
 if( fkeyread_double_in_section(fwzd, ini_section, "s_incr", &s_incr) ) s_incr=1e6; else s_incr*=1e6;
 if( fkeyread_double_in_section(fwzd, ini_section, "e_incr", &e_incr) ) e_incr=1.0e-3; else e_incr*=1e-2;
 if( fkeyread_double_in_section(fwzd, ini_section, "faz_incr", &faz_incr) ) faz_incr=0.02;
 if( fkeyread_int_in_section(fwzd, ini_section, "pause_factor", &pause_factor) ) pause_factor=0;
 if( fkeyread_int_in_section(fwzd, ini_section, "print_na_ekran", &na_ekran) ) na_ekran=1;
// if( fkeyread_int_in_section(fwzd, ini_section, "k_elements", &k_elements) ) k_elements=1;
//	if(k_elements<0 || k_elements>2) Exit_error("Error in wzd-file. Invalid value k_elements=%i. Must be 0, 1 or 2",k_elements );
 //BEGIN Open etap_ending_dat_files
 //FILE* fcycle[5]; declared in the beginning of main() 
 if( fkeyread_int_in_section(fwzd, ini_section, "how_many_etap_ending_dat_files", &k_cycle_files) ){
   k_cycle_files=0;
   //default value is assigned when fkeyread_int_.. returns 1 i.e. key is not found
 }
 else{
   fprintf(stderr,"\nk_cycle_files=%i",k_cycle_files);
   //pause();
   if(k_cycle_files>5) Exit_error("Error in wzd-file. Maximum 5 etap_ending_dat_files are allowed");
   if(ffindkey_in_section(fwzd, ini_section, "etap_numbers"))
     Exit_error("Key  etap_numbers=  not found after key  how_many_etap_ending_dat_files= ");
   for(ifil=0; ifil<k_cycle_files ; ifil++) fscanf( fwzd,"%i", etap_print_numbers+ifil);

   fprintf(stderr,"\netap_print_numbers are:");
   for(ifil=0; ifil<k_cycle_files ; ifil++) fprintf(stderr,"%i ", etap_print_numbers[ifil]);
   //pause();

   if(ffindkey_in_section(fwzd, ini_section, "etap_data_filenames"))
     Exit_error("Key  etap_data_filenames=  not found after key  how_many_etap_ending_dat_files= ");

   fprintf(stderr,"\ncorresponding filenames are:");
   for(ifil=0; ifil<k_cycle_files ; ifil++){
     char tname[FILENAME_LENGTH];
     fread_word(fwzd, tname, FILENAME_LENGTH);
     fprintf(stderr,"  %s",tname);
     fcycle[ifil] = fopen_next_dat(tname,"w");
     print_cycle_header(fcycle[ifil]);
     //pause("Begin printing");
     //fprintf(fcycle[ifil],
     //"\n %10li  %10.6lf %7.1lf % 12.6lf % 10.6lf % 12.6lf % 10.6lf  %8.6lf  % 12.6lf %5i %12.6lf %5i"
     //  ,1L  ,2.22   ,3.33    ,4.44   ,5.55    ,6.66    ,7.77    ,8.88    ,9.99   ,10  ,11.1  ,12);
   //cycle#   tim(s)    tem(K) e22/pct s22MPa   g32/pct  s32MPa   faza      RoMPa   iz  bf_max kzer_dmg);

     //close files to check their opening
     //fclose(fcycle[ifil]);
   }
   //pause("check opening files for printing results of etaps");

  }//else
 //END Open etap_ending_dat_files

 //BEGIN define what iterations of [cycle 1] to print

 //For cyclic loading (when [cycle 1] is activated ) results of all steps are
 //printed only for iterations 1,2,...,9,10,20,...,90,100,200,...,900,1000,2000,...
 //and for iterations defined by user
 //This rule is applied only for [cycle 1] not for inner embedded cycles.

 if( fkeyread_int_in_section(fwzd, ini_section, "how_many_cycles_to_print", &k_icyc_print) ){
   k_icyc_print=0;
   //default value is assigned when fkeyread_int_.. returns 1 i.e. key is not found
 }
 else{//i.e. when k_icyc_print!=0
   if(k_icyc_print>10) Exit_error("in wzd-file  how_many_cycles_to_print  must be <=10");
   if(ffindkey_in_section(fwzd, ini_section, "what_numbers"))
     Exit_error("key  what_numbers  not found after key  how_many_cycles_to_print ");
   //if(!numeral_follows(fwzd)) Exit_error("some numbers must be after key  what_numbers  in wzd-file");
   for(icyc=0; icyc<k_icyc_print ; icyc++) fscanf( fwzd,"%li",printed_cycles+icyc);

   printf("\nhow_many_cycles_to_print= %i",k_icyc_print);
   printf("\nNumbers of cycles to print are:");
   for(icyc=0; icyc<k_icyc_print ; icyc++) printf("  %li",printed_cycles[icyc]);
   pause(10);
 }
 //END define what iterations of [cycle 1] to print

//BEGIN declare and initialise different variables

   // printf("\n140:");
   // pause(10);

 VariantsOrientations Vrnts;
 Vrnts.Init(MC);

     //printf("\n150:");
     //pause(pause_factor);

 GrainsParameters GrPar;
 GrPar.Init();

 printf("\n160: VariantsOrientations and GrainsParameters Init passed successfully ");
 //pause();

//Isolation Helix; 
InternalParameters XXold, XXnew;
//InternalParameters& XXold=XXold1, XXnew=XXnew1; //when there is 1 element we use notation XXold, XXnew

//Isolation by round slotted springs
//Declaration of the vibroisolation arret unit and its initialisation
Porous_exVIA porous; 
porous.Init(tem,force,MC,GrPar);
macro_failure_pointer = &porous.macro_failure;

printf("\n160: porous.Init passed successfully ");
//pause();

	////file fmag_u.dat describes the dependence of magnetic force on the displacement
	//fmag_u = fopen("fmag_u.dat","w");
	//fprintf(fmagforce  ,"#file fmag_u.dat describes the dependence of magnetic force on the displacement");
	//fprintf(fmagforce  ,"\n\n#displ   fmag ");
	//fprintf(fmagforce    ,"\n#  mm      N   ");
	//{	Real u=-9.9e-3;
	//	while(u<10.0e-3){
	//		fprintf(fmag_u,"\n %+6.3lf  %+6.3lf",u*1e3 ,porous.F_mag(u,1.0));
	//		u+=1.0e-4;
	//	}
	//}
	//fclose(fmag_u);


//------------------- print of the initial value in the output file "ttesgtf.dat" ------------------------------
//------------------------------------------------------------------------------------------------------------
	print_step_header_2slotted_elems(fttesf);
				//print info on strain and stress in (iuz_max)-th SMA circular_beam
	//void print_step_2slotted_elems(FILE* fttesf ,int print_step_flag,  int na_ekran, int istep, double tim, double tem,
	//                 real force1, real force2, int Kuz ,Circular_beam* Cb1, Circular_beam* Cb2){
	print_step_2slotted_elems(fttesf ,1 ,1, 0, tim, porous.tem_old, porous.force_old, porous.kc_beams 
									,porous.c_beam ,porous.displ_new_1, porous.displ_new_2 ,porous.displ_new_3, porous.displ_new_4, porous.e_macro ,porous.stress_macro ,porous.force_new_eff);
	//fclose(fttesf); // if we want to look at file ttesfgtf.dat we must close it and exit
	//pause(); 
	//exit(0);

printf("\n\nInitial state13: tim=%lg tem=%lg  \n", tim, tem);
//pause();
 
 TASK cur_task=load, prev_task=load;//mode of task to be defined here
 // enum TASK {...};

//---------------------------------------------------------------------------------
//--------------------------------- OPEN CYCLE ON ETAPS ---------------------------
//------------|||||||--------------------|||||||-----------------||||||------------
//------------vvvvvvv--------------------vvvvvvv-----------------vvvvvv------------
//
 
 CYCLE cyc[4];//CYCLE is structure defined in "diaproc.h". We allow 4 embedded cycles.
 for(int jcyc=0;jcyc<4;jcyc++) cyc[jcyc].activated=0;//assign initial flag that cycle was not read
 fseek(fwzd,ini_section.pos1,SEEK_SET); // position to the end of ini_section
 
  //   printf("\nExper_95: ");
	 //pause();

 //Start reading sections in wzd-file
 DATA_SECTION cur_section;
 while(freadnextsection(fwzd,cur_section)!=EOF && *macro_failure_pointer==0)
 { // calculate while wzd-file is not exhausted && no macro-fracture happens
   //cur_section is defined inside the condition of while-cycle
 
   ETAP cur_etap;
     //enum ETAP_TYPE {etap,cycle,endcycle};
     //struct ETAP { ETAP_TYPE type; int num; };
     //enum ETAP_TYPE {etap,cycle,endcycle};

   //Recognize cur_etap (cur_section was read in the condition of while-cycle)
   if (recognize_etap(cur_etap, cur_section) ){ //conditional expr == 0 <=> etap is well recognized 
     fprintf(stderr,"\n\nERROR: Unknown etap %s in wzd-file",cur_section.name);
     fprintf(stderr,"\n Program is terminated.");
	 pause();
     exit(1);
   }

  //   printf("\nExper_100: ");
	 //pause();

     switch (cur_etap.type) {
       case 0:   printf("etap "); break;
       case 1:   printf("cycle "); break;
       case 2:   printf("endcycle "); break;
       default:  Exit_error("Unknown etap");
     }
     printf("%i",cur_etap.num);
     //pause();

   if(cur_etap.type==cycle){
     int i_cycnum=cur_etap.num-1;//cycles in wzd-file are numbered from 1 to 4 and in arrays from 0 to 3
     if(i_cycnum>3 || i_cycnum<0){
       fprintf(stderr,"\n\nERROR: invalid cycle number (it must be integer from 1 to 4(maximum))");
       fprintf(stderr,"\n Program is terminated.");
       pause();
	   exit(2);
     }
     //check that cycle with index i_cycnum is not activated and previous is activated
     if(cyc[i_cycnum].activated==1)Exit_error("cycle number %i already used",i_cycnum+1);
     if(i_cycnum>0) if(cyc[i_cycnum-1].activated==0)Exit_error("cycle %i must be opened first",i_cycnum);
     cyc[i_cycnum].activated=1;//activate cycle

     //here cur_etap.type==cycle and the only data in this section is k_cycle
	 if( fkeyread_long_in_section(fwzd, cur_section, "k_cycle", &(cyc[i_cycnum].k)) ){
		printf("\n\n Number of cycle iterations must be specified after  k_cycle=");
		Exit_error("Cannot find key k_cycle= in section [cycle %i] in wzd-file",i_cycnum+1);
     }
     cyc[i_cycnum].pos0 = ftell(fwzd);
     cyc[i_cycnum].i = 0L;
     continue;//go to next etap
   }

   if(cur_etap.type==endcycle){
     int i_cycnum=cur_etap.num-1;//cycles in wzd-file are numbered from 1 to 4
     if(i_cycnum>3||i_cycnum<0){
       fprintf(stderr,"\n\nERROR: Invalid number in section [endcycle %i]",i_cycnum+1);
       fprintf(stderr,"\n Program is terminated.");
       pause();
	   exit(2);
     }
     //check that cycle with index i_cycnum is activated and with index i_cycnum+1 is not activated
     if(cyc[i_cycnum].activated==0)Exit_error("invalid number in [endcycle %i]",i_cycnum+1);
     if(i_cycnum<3) if(cyc[i_cycnum+1].activated==1)Exit_error("cycle %i must be closed first",i_cycnum+2);

     //print the results of the most outer cycle 
     //if(i_cycnum==0) print_etap(fcyc1, cyc[i_cycnum] ,tim ,tem ,sig ,eps ,XXnew);

     //if the number of iterations < cyc[i_cycnum].k then change the position to the beginning of the cycle:
     if(++cyc[i_cycnum].i < cyc[i_cycnum].k) fseek(fwzd, cyc[i_cycnum].pos0, SEEK_SET);
     else {cyc[i_cycnum].activated=0; //deactivate cycle
     //a reset to 0 of the iteration count is done in the beginning of the cycle
	  fseek(fwzd, cur_section.pos1, SEEK_SET);
     }
     continue;//go to next etap (next iteration of the while-cycle):
   }

   //Now cur_etap.type==etap as it cannot be otherwise
   if(!cyc[0].activated) print_step_flag = 1; 
   //if there are no opened (activated) cycles we shall print results for all steps
   else //i.e. if there are opened cycles 
   // we print steps only when the number of the most outer cycle satisfies conditions 1 and 2
   { long j=cyc[0].i+1L; //current iteration of the cycle
     int condition1 = j<10L || j==10L || j==20L || j==50L
                           || j==100L || j==200L || j==500L
                           || j==1000L || j==2000L || j==5000L
                           || j==10000L || j==20000L || j==50000L
                           || j==100000L || j==200000L || j==500000L
                           || j==1000000L || j==2000000L || j==5000000L
                           || j==10000000L || j==20000000L || j==50000000L
                           || j==100000000L || j==200000000L || j==500000000L
                           || j==1000000000L || j==2000000000L;
     int condition2=0;
     for(icyc=0; icyc<k_icyc_print ; icyc++) if(j==printed_cycles[icyc]){ condition2=1; break;}
     print_step_flag = condition1 || condition2 ? 1 : 0 ;
   }
    //admissible_keys=(admissible values):
   //tim=   htim=   tem=    htem=   dtem=  press=  hpress=   dpress=
   //s33=   hs33=   ds33=   s32=    hs32=  s32=    
   //e33=   he33=   de33=   e32=    he32=  de32=   zhest33=  zhest32=
   //fluence=  hfluence=
   //k_step=  k_print=
   //ultrasound33=(0/1)
   //vibr32=(0/1)
   int k_print;//number of printings for current etap
   int k_step;//number of steps for current etap
   //in "diaproc.h": enum TASK {...};
   double tim_z, tem_z, force_z, s33_z, s32_z, e33_z, e32_z, press_z;//finally they == incr.for etap
    double hsig_33, heps_33; // я дописала
   //press_z is merely added to all diagonal components of stress tensor
   //flags for the presence of corresponding keys in wzd-file:
   int tim_flag, htim_flag, any_tim_flag, 
	   tem_flag,  htem_flag, dtem_flag, any_tem_flag,
	   //tem2_flag,  htem2_flag, dtem2_flag, any_tem2_flag,
       s33_flag, hs33_flag, ds33_flag, any_s33_flag,  s32_flag, hs32_flag, ds32_flag, any_s32_flag,
       e33_flag, he33_flag, de33_flag, any_e33_flag,  e32_flag, he32_flag, de32_flag, any_e32_flag,
       press_flag,   hpress_flag,   dpress_flag,   any_press_flag,   zhest33_flag, zhest32_flag,
       fluence_flag, hfluence_flag, dfluence_flag, any_fluence_flag,
       k_step_flag, ultrasound33_flag ,vibr32_flag ,quake32_flag ,vibr33_flag //,quake33_flag
       ,siga_us33_flag, epsa_us33_flag;
   //note: !fkeyread...(...)==1 on success and 0 - when no key was found
   //if key tim= is found we calculate tim_z (incr.for etap) by substraction:

     //printf("\n160:");
     //pause();

   tim_z=100.0;//default value if both keys tim= and htim= are absent
   if( (tim_flag=!fkeyread_double_in_section(fwzd,cur_section,"tim",&tim_z))!=0) tim_z -= tim;

     //printf("\n162:tim_flag=%i tim_z=%lg" ,tim_flag ,tim_z);
     //pause();

   htim_flag = !fkeyread_double_in_section(fwzd, cur_section, "htim", &tim_z);
   any_tim_flag = tim_flag || htim_flag;

     printf("\n180:tim_flag=%i htim_flag=%i any_tim_flag=%i tim_z=%lg" ,tim_flag ,htim_flag ,any_tim_flag ,tim_z);
     //pause(pause_factor);

   if (tim_flag && htim_flag) Exit_error("conflict of keys tim= and htim= in [etap %i] in wzd-file",cur_etap.num);
   if(tim_z<=0) Exit_error("Invalid time in [etap %i] in wzd-file (time increment<=)",cur_etap.num);

     printf("\n200:htim_flag=%i tim_z=%lg",htim_flag,tim_z);
     //pause(pause_factor);

   k_step_flag = !fkeyread_int_in_section(fwzd, cur_section,"k_step",&k_step);
   if(!k_step_flag) k_step=50;//If there is no key k_step it will be calculated. This is its minimum value
   if(fkeyread_int_in_section(fwzd, cur_section,"k_print",&k_print)) k_print=0;//default - print result for each step

    if(fkeyread_int_in_section(fwzd, cur_section,"load_flag",&load_flag)) load_flag=0;
	if(fkeyread_int_in_section(fwzd, cur_section,"preload_flag",&preload_flag)) preload_flag=0;//default value is 0
	//preload_flag=1 means independent loading and heating of RSS   
	if(fkeyread_int_in_section(fwzd, cur_section,"add_magnetic_force_flag",&add_magnetic_force_flag))add_magnetic_force_flag=0;//default-0
	//add_magnetic_force_flag=1 is used when we introduce magnetic force (default - 0)
	if(fkeyread_int_in_section(fwzd, cur_section,"vibr_trig_polynom_flag",&vibr_trig_polynom_flag)) vibr_trig_polynom_flag=0;
	//vibr_trig_polynom_flag=1 is used when excitation is a trigonometric polynom (default - 0)
	if(fkeyread_int_in_section(fwzd, cur_section,"idle_flag",&idle_flag)) idle_flag=0;
	//vibr_trig_polynom_flag=1 is used when excitation is a trigonometric polynom (default - 0)
	operation_flag_prev=operation_flag;
	operation_flag = vibr_trig_polynom_flag || idle_flag;
	
	if(				 load_flag && (add_magnetic_force_flag || vibr_trig_polynom_flag || idle_flag)|| preload_flag && (add_magnetic_force_flag || vibr_trig_polynom_flag || idle_flag)
	    || add_magnetic_force_flag && (vibr_trig_polynom_flag || idle_flag)
		|| vibr_trig_polynom_flag && idle_flag
	  )//one and only one key of these must be present in the current stage (etap) of the wzd-file
	{
	  fprintf(stderr,"\n\nThere must be one and only one of the keys: ");
	  fprintf(stderr,"\nload_flag= ,preload_flag= ,vibr_trig_polynom_flag= ,add_magnetic_force_flag= ,idle_flag= in [etap%i] in wzd-file", cur_etap.num);
	  fprintf(stderr,"\nin [etap%i] in wzd-file", cur_etap.num);
	  fprintf(stderr,"\n\nProgram will be terminated");
	  pause();
	  exit(1);
	}
	
	//if( fkeyread_word_in_section(fwzd ,cur_section ,"excitation" ,excitation ,MALBUF) ) strcpy(excitation,"no_excitation");// - default value
	// default - excitation == "no_excitation". Other possible value "trig_polynom"
   if( (tem_flag=!fkeyread_double_in_section(fwzd, cur_section,"tem",&tem_z))!=0 ) tem_z -= tem;
   if((dtem_flag=!fkeyread_double_in_section(fwzd, cur_section, "dtem", &tem_z))!=0)tem_z *=tim_z;
   htem_flag = !fkeyread_double_in_section(fwzd, cur_section, "htem", &tem_z);
   any_tem_flag = tem_flag || dtem_flag || htem_flag;
   if(!any_tem_flag) tem_z=0.0;//default value

   if(tem_flag && htem_flag || tem_flag && dtem_flag || htem_flag && dtem_flag)
     Exit_error("Conflict of keys tem=, htem=, dtem in [etap%i] in wzd-file", cur_etap.num);

     printf("\n300:any_tem_flag=%i tem_z=%lg",any_tem_flag,tem_z);
     //pause(pause_factor);

   //if( (tem2_flag=!fkeyread_double_in_section(fwzd, cur_section,"tem2",&tem2_z))!=0 ) tem2_z -= tem2;
   //if((dtem2_flag=!fkeyread_double_in_section(fwzd, cur_section, "dtem2", &tem2_z))!=0)tem2_z *=tim_z;
   //htem2_flag = !fkeyread_double_in_section(fwzd, cur_section, "htem2", &tem2_z);
   //any_tem2_flag = tem2_flag || dtem2_flag || htem2_flag;
   //if(!any_tem2_flag) tem2_z=0.0;//default value

   //if(tem2_flag && htem2_flag || tem2_flag && dtem2_flag || htem2_flag && dtem2_flag)
   //  Exit_error("Conflict of keys tem2=, htem=2, dtem 2in [etap%i] in wzd-file", cur_etap.num);

   //  printf("\n300:any_tem2_flag=%i tem2_z=%lg",any_tem2_flag,tem2_z);
     //pause(pause_factor);

//==============================================================================================================
//==============================================================================================================
//										Determine cur_task: 
//==============================================================================================================
//==============================================================================================================
	prev_task=cur_task;//memorize task on the previous etep
	
	if(load_flag) cur_task=load;
	else if(preload_flag) cur_task=preload;
	else if(add_magnetic_force_flag) cur_task= add_magnetic_force;
	else if(vibr_trig_polynom_flag) cur_task= vibr_trig_polynom;
	else if(idle_flag) cur_task= idle;
	else Exit_error("Exper180: Cannot recognize the task for [etap %i] in wzd-file", cur_etap.num);
	
	printf("\n\nExper195: determined cur_task=%i",cur_task);
	//pause();


  tim_ke = tim + tim_z;//time when etap (stage) will end
  Real tim_ke_0 = tim_ke - 1e-10;//almost time at which a stage ends
 
 // tim_etap=0.0;
    s33_us=0.0;
    e33_us=0.0;

	

  printf("\n\nStarted etap %i: cur_task= %i", cur_etap.num ,cur_task);
  printf("\ntim=%lg tim_z=%lg tim_ke=%lg tem=%lg  tem_z=%lg",tim, tim_z ,tim_ke, tem, tem_z);
  pause(0);

	//printf("\n\nExper197: before switch cur_task=%i",cur_task);
	//pause();
  	XXold.Init(MC ,GrPar);
	XXnew.Init(MC ,GrPar);


  switch(cur_task){

   case load:
	   
	   int k_print;
	   Real hforce, force_z;
	   Real hdispl;
	   Real heps_33_1;
	   Real hsig_33_1;
	   Real heps_33_2;
	   Real hsig_33_2;
	   Real hfaz_1;
	   Real hfaz_2;
	   Real heps_33_3;
	   Real hsig_33_3;
	   Real hfaz_3;
	   Real heps_33_4;
	   Real hsig_33_4;
	   Real hfaz_4;
	   Real heps_33_5;
	   Real hsig_33_5;
	   Real hfaz_5;
	   Real heps_33_6;
	   Real hsig_33_6;
	   Real hfaz_6;
	   Real htem;
	   Real hfluence;


    printf("\nLOADING");
	if(fkeyread_double_in_section(fwzd, cur_section, "force", &force_z)) force_z=force;
	else hforce=0.0;
	if(k_step_flag) f_step = 1.0/(Real)k_step;
	else if(force_z!=0.0) f_step = 0.2/fabs(force_z); //there was "0.002N (200mg) per step", but then the cycle is too long
	else if(tem_z>1e-12) {f_step = 1.0/fabs(tem_z);
	                     }
	else f_step=5e-3;//divide etap in 200 steps
    pause(pause_factor);
	

    for(istep=0; tim < tim_ke && *macro_failure_pointer==0; istep++){
     Real tem_ke=tem+tem_z;  
       int good_step=0, iincrease_step=0, idecrease_step=0 ;


	   if (tim > 6.31579)
		   good_step = good_step;


       //---------------------------
			//int last_step=0; 
	while( !good_step ) {
	   if(tim+htim>tim_ke){htim=tim_ke-tim; f_step=htim/tim_z;}
	        htim = tim_z * f_step;
			//f_step = tim_z * htim;
			porous.htem = tem_z * f_step;
			porous.tem_new = porous.tem_old + porous.htem;
			//pause();
/////////////////////////////////////////////////////////////////////		    hfluence = fluence_z * f_step;
			hforce = force_z * f_step; 
			porous.hforce = hforce;
			porous.force_new = porous.force_old+porous.hforce;
			porous.displ_new = porous.displ_old + porous.hdispl;
            porous.displ_new_1 = porous.displ_old_1 + porous.hdispl_1;
            porous.displ_new_2 = porous.displ_old_2 + porous.hdispl_2;
			porous.displ_new_3 = porous.displ_old_3 + porous.hdispl_3;
			porous.displ_new_4 = porous.displ_old_4 + porous.hdispl_4;
			
		//printf("\n900:i_step=%i tim=%lg htim=%lg",istep ,tim ,htim);
		//pause();

			//sigma = sig;                      hsigma = hsig;
			//sigma.Add_to_diagonal(-press);    hsigma.Add_to_diagonal(-hpress);// 
			//sigma(2,2)+=s33_us;               hsigma(2,2)+=hs33_us;
			
			//correction of stress due to the decrease of the cross-section 
			// !! Valid only if the applied force is constant !!
			//	hsigma(2,2) = sigma(2,2) * heps(2,2)/(1.0 + eps(2,2));
			//end of the correction of stress due to the decrease of the cross-section 
			printf("\nsig(1,1)=%lg, \ne22_1=%lg", /*porous.c_beam[0].K_sig * porous.force_old*/porous.c_beam[0].sig_to_control*1e-6, porous.c_beam[0].XXnew.TotalStrain.C(1, 1)*1e2);
			
			
		//printf("\nBefore Task_sig:i_step=%i tim=%lg",istep,tim);
		//pause();

			//Task_VIA_operation(istep, htim , hdispl, porous	,MC	,Vrnts	, GrPar);

			Task_force_tem(istep, htim, tem, sig ,hsig ,sig_elastic ,hsig_elastic ,eps ,heps
				,MC , porous, Vrnts ,GrPar );
		/*	printf ("\n sig(2,2)=%lg, hsig(2,2)=%lg"  ,sig(2,2), hsig(2,2)  );
		pause();*/
		//printf("\nAfter Task_sig:i_step=%i tim=%lg",istep,tim);
		//pause();
			printf("\nsig(2,2)=%lg, \ne33_1=%lg", /*porous.c_beam[0].K_sig * porous.force_old*/porous.c_beam[0].sigma33*1e-6, porous.c_beam[0].XXnew.TotalStrain.C(2, 2)*1e2);
			//if(k_step_flag || ultrasound33_flag || last_step) break;//steps will not be adjusted// commentttt
		//---step adjusting ----
			//hsig=0.0; hsig(2,2) = porous.c_beam_max->K_sig * (+hforce); - коммент, так как то же в ТР
			//Real hfaz=XXnew.Phase-XXold.Phase;
			//hfaz=fabs(porous.c_beam_max->XXnew.Phase-porous.c_beam_max->XXold.Phase);
			htem=fabs(porous.htem);
			hfaz_1 = fabs(porous.c_beam[0].XXnew.Phase - porous.c_beam[0].XXold.Phase);
			hfaz_2 = fabs(porous.c_beam[1].XXnew.Phase - porous.c_beam[1].XXold.Phase);
			hfaz_3 = fabs(porous.c_beam[2].XXnew.Phase - porous.c_beam[2].XXold.Phase);
			hfaz_4 = fabs(porous.c_beam[3].XXnew.Phase - porous.c_beam[3].XXold.Phase);
			hfaz_5 = fabs(porous.c_beam[4].XXnew.Phase - porous.c_beam[4].XXold.Phase);
			hfaz_6 = fabs(porous.c_beam[5].XXnew.Phase - porous.c_beam[5].XXold.Phase);
			hsig_33_1 = fabs(porous.c_beam[0].K_sig * hforce); 
			hsig_33_2 = fabs(porous.c_beam[1].K_sig * hforce);
			hsig_33_3 = fabs(porous.c_beam[2].K_sig * hforce);
			hsig_33_4 = fabs(porous.c_beam[3].K_sig * hforce);
			hsig_33_5 = fabs(porous.c_beam[4].K_sig * hforce);
			hsig_33_6 = fabs(porous.c_beam[5].K_sig * hforce);
			heps_33_1 = fabs(porous.c_beam[0].XXnew.TotalStrain.C(2,2) - porous.c_beam[0].XXold.TotalStrain.C(2,2));
			heps_33_2 = fabs(porous.c_beam[1].XXnew.TotalStrain.C(2,2) - porous.c_beam[1].XXold.TotalStrain.C(2,2));
			hsig_33_3 = fabs(porous.c_beam[2].K_sig * hforce);
			heps_33_3 = fabs(porous.c_beam[2].XXnew.TotalStrain.C(2, 2) - porous.c_beam[2].XXold.TotalStrain.C(2, 2));
			hsig_33_4 = fabs(porous.c_beam[3].K_sig * hforce);
			heps_33_4 = fabs(porous.c_beam[3].XXnew.TotalStrain.C(2, 2) - porous.c_beam[3].XXold.TotalStrain.C(2, 2));
			hsig_33_5 = fabs(porous.c_beam[4].K_sig * hforce);
			heps_33_5 = fabs(porous.c_beam[4].XXnew.TotalStrain.C(2, 2) - porous.c_beam[4].XXold.TotalStrain.C(2, 2));
			hsig_33_6 = fabs(porous.c_beam[5].K_sig * hforce);
			heps_33_6 = fabs(porous.c_beam[5].XXnew.TotalStrain.C(2, 2) - porous.c_beam[5].XXold.TotalStrain.C(2, 2));
            good_step=1;
			if ((hsig_33_1>s_incr || hsig_33_2>s_incr || hsig_33_3>s_incr || hsig_33_4>s_incr || hsig_33_5>s_incr || hsig_33_6>s_incr) || (heps_33_1 >e_incr || heps_33_2>e_incr || heps_33_3>e_incr || heps_33_4>e_incr || heps_33_5>e_incr || heps_33_6>e_incr) || (hfaz_1>faz_incr && hfaz_2>faz_incr && hfaz_3>faz_incr && hfaz_4>faz_incr && hfaz_5>faz_incr && hfaz_6>faz_incr) || htim>tim_incr || htem>tem_incr)
			{ if(idecrease_step++<MAX_STEP_DECREASE){ htim /= 2.0; good_step=0;} }
			else if ((hsig_33_1<0.2*s_incr || hsig_33_2<0.2*s_incr || hsig_33_3<0.2*s_incr || hsig_33_4<0.2*s_incr || hsig_33_5<0.2*s_incr || hsig_33_6<0.2*s_incr) && (heps_33_1<0.2*e_incr || heps_33_2<0.2*e_incr || heps_33_3<0.2*e_incr || heps_33_4<0.2*e_incr || heps_33_5<0.2*e_incr || heps_33_6<0.2*e_incr) && (hfaz_1<0.2*faz_incr || hfaz_2<0.2*faz_incr || hfaz_3<0.2*faz_incr || hfaz_4<0.2*faz_incr || hfaz_5<0.2*faz_incr || hfaz_6<0.2*faz_incr) && htim<0.2*tim_incr && htem<0.2*tem_incr)
			{ if(iincrease_step++<0){ htim *= 2.0; good_step=0;} }//increase step was 2, not 0!!!


			/*if(hsig_33_2>s_incr && hsig_33_2>s_incr||heps_33_2>e_incr||hfaz>faz_incr||htim>tim_incr||htem>tem_incr)
			{ if(idecrease_step++<MAX_STEP_DECREASE){ htim /= 2.0; good_step=0;} }
		    else if(hsig_33_2<0.2*s_incr && heps_33_2<0.2*e_incr && hfaz<0.2*faz_incr && htim<0.2*tim_incr && htem<0.2*tem_incr)
			{ if(iincrease_step++<2){ htim *= 2.0; good_step=0;} }*/
		//	hsig=0.0; hsig(2,2) = porous.c_beam_max->K_sig * (+hforce);//the 1st RSS is pulled to the right 
			//hsig2=0.0; hsig2(2,2) = porous.c_beam2_max->K_sig * (-hforce);//the 2nd RSS is pulled to the left 
			printf("\nhetap_force=%lg,   \ndispl_new2=%lg,   \nhdissspl_2=%lg,   \ndispl_old2=%lg,           \nsig(2,2)=%lg,                               \ne33_1=%lg,                                    \ne33_2=%lg,          \ne33_3=%lg,                          \ne33_4=%lg,      \ne33_5=%lg, \ne33_6=%lg,                    \neps_macro=%lg,   \nsigma_macro=%lg   \nforce_real=%lg"
				, force_z, porous.displ_new_2, porous.hdispl_2, porous.displ_old_2, porous.c_beam[1].K_sig * porous.force_old*1e-6, porous.c_beam[0].XXnew.TotalStrain.C(2, 2)*1e2, porous.c_beam[1].XXnew.TotalStrain.C(2, 2)*1e2, porous.c_beam[2].XXnew.TotalStrain.C(2, 2)*1e2, porous.c_beam[3].XXnew.TotalStrain.C(2, 2)*1e2, porous.c_beam[4].XXnew.TotalStrain.C(2, 2)*1e2, porous.c_beam[5].XXnew.TotalStrain.C(2, 2)*1e2, porous.e_macro*1e2, porous.stress_macro*1e-6, porous.force_new_eff);
			
			//porous.hforce = +hforce;
			//porous.force_new=porous.force_old+porous.hforce;
			//htem=fabs(porous.htem);// ??
			//Real hfaz=fabs(porous.c_beam_max->XXnew.Phase - porous.c_beam_max->XXold.Phase);
			//hsig_33=fabs(porous.c_beam_max->K_sig * (+hforce));
			//heps_33 = fabs( porous.c_beam_max->XXnew.TotalStrain.C(2,2) - porous.c_beam_max->XXold.TotalStrain.C(2,2) );
			//good_step=1;
		
       }//for int good_step=0  //this cycle will not be done more than 8 times: 
        //htim may be increased only in <= 2 times and decreased <= 2**7 times
       //----------------------                      
       //Update parameters:
       tim += htim; tem=tem_ke;; force = force_z; sig+=hsig; eps+=heps; fluence += hfluence; // the force has been changed. (instead of sig)
       tim_etap+=htim;// s33_us += hs33_us;
       //XXold = XXnew;
       porous.c_beam_max->XXold = porous.c_beam_max->XXnew;
	   porous.Update();
	   real displ_1 = porous.displ_new_1;
	   real displ_2 = porous.displ_new_2;
	   real displ_3 = porous.displ_new_3;
	   real displ_4 = porous.displ_new_4;

      ////////////////////////////////////////////////////////////////// if(ultrasound33_flag && !aliquot(tim_etap,semiperiod_us33,1e-2) ) continue;
       //Function aliquot(x,y,eps)(x,y,eps>0) returns 1 if exists N : fabs(x-N*y)< eps*y
       ////////////////////////////////////////if(!k_print || tim>=tim_print){//print results (if k_print==0 - print all steps)
        /////////////////////////////////// tim_print+=htim_print;//find new time, when to print results
	 //print_step(fttesf ,print_step_flag ,na_ekran, istep, tim, tem, sig, force, eps, porous.c_beam_max->XXnew); // force<-sig
	 print_step_2slotted_elems(fttesf ,print_step_flag ,0 ,istep ,tim ,porous.tem_old 
		 , porous.force_old, porous.kc_beams, porous.c_beam, porous.displ_new_1, porous.displ_new_2, porous.displ_new_3, porous.displ_new_4, porous.e_macro, porous.stress_macro, porous.force_new_eff);
       ///////////////////////////////////////////////////////////}//if() printing results
    }// for int istep while tim < tim_ke
    break;//from switch case load:

	case preload: 
{
	//printf("\n\nExper215: started case preload");
	//pause();

	int k_print=0;//now we print results for all steps
	int elongation_flag=1;
	int control_of_elongation;
	//Real helongation_of_elem_z;
	//Real helongation_of_elem_cur;
	Real hforce, force_z;
//	Real tim_etap;
	Real hdispl;
	//Matr33 heps1, heps2;
	Real heps_33_1;
	Real hsig_33_1;
	Real heps_33_2;
    Real hsig_33_2;
	Real hfaz_1;
	Real hfaz_2;
	Real htem;
	//if(fkeyread_double_in_section(fwzd, cur_section, "max_helongation_of_elem", &helongation_of_elem_z)){
	//	helongation_of_elem_z=0.0;//default
	//	control_of_elongation=0;//no control of elongation
	//}
	//else control_of_elongation=1;//if any elongation (zero or non-zero) is specified, steps will stop when it is reached
	//helongation_of_elem_z*=1e-3; //mm -> m
	if(fkeyread_double_in_section(fwzd, cur_section, "hforce", &force_z)) force_z=0.0;//default
	force_z*=1e-2; // g->N //already in Newtons

	//Real tem_z = fabs(tem_z); 
	if(k_step_flag) f_step = 1.0/(Real)k_step;
	else if(tem_z>1e-12) {
		f_step = 1.0/fabs(tem_z);
	}
	else if(force_z!=0.0) f_step = 0.002/force_z; //0.002N (200mg) per step
	else f_step=5e-3;//divide etap in 200 steps

	Real tem_ke=tem+tem_z;
	//Real tem2_ke=tem2+tem2_z;
	//tim_ke=tim+tim_z; this was done before for all cases
	printf("\n\nExper780: case preload");
	printf("\ntim=%lg tim_z=%lg tim_ke=%lg",tim,tim_z,tim_ke);
	printf("\ntem=%lg tem_z=%lg tem_ke=%lg",tem,tem_z,tem_ke);
	//printf("\ntem2=%lg tem2_z=%lg tem2_ke=%lg",tem2,tem2_z,tem2_ke);
	printf("\nf_step=%lg", f_step);
	printf("\nhetap_force=%lg",force_z);
	//pause(); //pause(pause_factor);
	//if(control_of_elongation) elongation_flag = helongation_of_elem_cur < helongation_of_elem_z;
	
//================================================================================
//========================Open cycle on steps=====================================
//================================================================================
	tim_etap = 0.0;
	//helongation_of_elem_cur=0.0;
	//if(control_of_elongation) elongation_flag = helongation_of_elem_cur < helongation_of_elem_z;
	for(istep=0; tim<tim_ke && *macro_failure_pointer==0 ; istep++){
		//if(control_of_elongation) elongation_flag = helongation_of_elem_cur < helongation_of_elem_z;
		//else elongation_flag=1;
		//---step adjusting cycle----

	//printf("\n800:i_step=%i htim=%lg",istep,htim);
	//pause();
		int good_step=0, iincrease_step=0, idecrease_step=0;
		while( !good_step ) {
			if(tim+htim>tim_ke){ htim=tim_ke-tim; f_step=htim/tim_z; }
			
			htim = tim_z * f_step;
			porous.htem = tem_z * f_step;
			//porous.htem2 = tem2_z * f_step;
			porous.tem_new = porous.tem_old + porous.htem;
			//porous.tem2_new = porous.tem2_old + porous.htem2;
			
			hforce=force_z*f_step;
			porous.hforce = hforce;
			//porous.hforce2 = -hforce;//the 2-nd slotted element is pulled to the left
			porous.force_new=porous.force_old+porous.hforce;
			//porous.force2_new=porous.force2_old+porous.hforce2;
			
			//printf("\n\nExper850:i_step=%4i f_step=%lg",istep,f_step);
			//printf("\ntim=%lg tim_z=%lg tim_ke=%lg htim=%lg ",tim,tim_z,tim_ke,htim,htim);
			//printf("\ntem=%lg tem_z=%lg tem_ke=%lg htem=%lg ",tem,tem_z,tem_ke,htem);
			//printf("\ntim_z*f_step=%lg tem_z*f_step=%lg",tim_z*f_step ,tem_z*f_step);
			//pause();
			/*Task_force_2_slotted_elements
				(
					htim ,hdispl	,porous	,MC	,Vrnts	,GrPar	
				);*/

			Task_force_tem(istep, htim, tem, sig, hsig, sig_elastic, hsig_elastic, eps, heps
				,MC , porous, Vrnts ,GrPar);
			//---step adjusting ----
			//Real hfaz1=c_beam1[iuz_max].XXnew.Phase - c_beam1[iuz_max].XXold.Phase;
			
			hsig=0.0; hsig(2,2) = porous.c_beam_max->K_sig * hforce;//the 1st RSS is pulled to the right  -------------------коммент??это в тп тоже
			//hsig2=0.0; hsig2(2,2) = porous.c_beam2_max->K_sig * (-hforce);//the 2nd RSS is pulled to the left 
			
			htem=fabs(porous.htem);
			//htem2=fabs(porous.htem2);
			hfaz_1 = fabs(porous.c_beam[0].XXnew.Phase - porous.c_beam[0].XXold.Phase);
			hfaz_2 = fabs(porous.c_beam[1].XXnew.Phase - porous.c_beam[1].XXold.Phase);
			hsig_33_1 = fabs(porous.c_beam[0].K_sig * (+hforce));
			//hsig2_33=fabs(porous.c_beam2_max->K_sig * (-hforce));
			heps_33_1 = fabs( porous.c_beam[0].XXnew.TotalStrain.C(2,2) - porous.c_beam[0].XXold.TotalStrain.C(2,2) );
			//heps2_33 = fabs( porous.c_beam2_max->XXnew.TotalStrain.C(2,2) - porous.c_beam2_max->XXold.TotalStrain.C(2,2) );
			hsig_33_2 = fabs(porous.c_beam[1].K_sig * (+hforce));
			//hsig2_33=fabs(porous.c_beam2_max->K_sig * (-hforce));
			heps_33_2 = fabs( porous.c_beam[1].XXnew.TotalStrain.C(2,2) - porous.c_beam[1].XXold.TotalStrain.C(2,2) );
			good_step=1;

			/*if((hsig_33_1>s_incr || hsig_33_2>s_incr)||(heps_33_1 >e_incr || heps_33_2>e_incr)||(hfaz_1>faz_incr || hfaz_2>faz_incr)||htim>tim_incr||htem>tem_incr)
			{ if(idecrease_step++<MAX_STEP_DECREASE){ htim /= 2.0; good_step=0;} }
		    else if((hsig_33_1<0.2*s_incr || hsig_33_2<0.2*s_incr) && (heps_33_1<0.2*e_incr || heps_33_2<0.2*e_incr) && (hfaz_1<0.2*faz_incr || hfaz_2<0.2*faz_incr) && htim<0.2*tim_incr && htem<0.2*tem_incr)
			{ if(iincrease_step++<2){ htim *= 2.0; good_step=0;} }*/

			if(	   hsig_33_1>s_incr||heps_33_1>e_incr||hfaz_1>faz_incr||htem>tem_incr||htim>tim_incr )
				{ if(idecrease_step++<MAX_STEP_DECREASE){ f_step/=2.0; good_step=0;} }
			else if(hsig_33_1<0.1*s_incr && heps_33_1<0.1*e_incr && hfaz_1<0.1*faz_incr && htem<0.1*tem_incr && htim<0.1*tim_incr )
				{ if(iincrease_step++<2){ f_step*=2.0; good_step=0;} }

			if(	   hsig_33_2>s_incr||heps_33_2>e_incr||hfaz_2>faz_incr||htem>tem_incr||htim>tim_incr )
				{ if(idecrease_step++<MAX_STEP_DECREASE){ f_step/=2.0; good_step=0;} }
			else if(hsig_33_2<0.1*s_incr && heps_33_2<0.1*e_incr && hfaz_2<0.1*faz_incr && htem<0.1*tem_incr && htim<0.1*tim_incr )
				{ if(iincrease_step++<2){ f_step*=2.0; good_step=0;} }

		}//while good_step=0  //this cycle will not be done more than 8 times:
		//htim may be increased only in <= 2 times and decreased <= 2**7 times
		//----------------------                      
		//Update parameters:
		//helongation_of_elem_cur += fabs(hdispl); 
		tim += htim; tim_etap+=htim;
		tem=tem_ke;
        porous.Update();
	
		if(istep%20==0) printf("\rE235: istep=%i hdispl=%lg mm force=%lg N hforce=%lg N "
									   ,istep   ,hdispl*1e3 ,porous.force_new ,hforce );
		//pause();

		//void print_step_2slotted_elems(FILE* fttesf ,int print_step_flag, int na_ekran, int istep, double tim, double tem,
		//                 real force1, real force2, int Kuz ,Circular_beam* Cb1, Circular_beam* Cb2){

		print_step_2slotted_elems(fttesf ,print_step_flag ,0 ,istep ,tim ,porous.tem_old 
			, porous.force_old, porous.kc_beams, porous.c_beam, porous.displ_new_1, porous.displ_new_2, porous.displ_new_3, porous.displ_new_4, porous.e_macro, porous.stress_macro, porous.force_new_eff);
														//na_ekran	
	}// for istep 
	//printf("\nE940:At the end of etap %3i: cur_task#%i hetap_elongation=%lg mm", cur_etap.num ,cur_task,helongation_of_2_elems_cur*1e3);
	//printf("\n force1, force2 = %lg %lg N", porous.force1_old,porous.force2_old);
	//pause();
	//if(control_of_elongation && elongation_flag){//helongation_of_2_elems_cur < helongation_of_2_elems_z i.e.required elongation is not reached
	//	printf("\nExper960: required elongation %lg is not reached." ,helongation_of_elem_z*1e3);
	//	printf("\nTo reach desired elongation set bigger hforce.");
	//	pause();
	//}
}
	break;//from switch case preload:

//================================================================================================================

	case add_magnetic_force : 
	case vibr_trig_polynom :
	case idle :
  {//begin of unconditional block surrounding this case
	int i_harmc;
	Real mass;
	Real Dx; // temp. variable for reading
	Real external_damping;
	Real hmagnetic_factor;
	printf("\n\nExper 712: case operation started");
	  
	//files for registering vibrograms and magnetic forces (these files were opened in the beginning of Exper)
	//fvibr = fopen("vibr.dat","w");
	//fprintf(fvibr  ,"#tim_etap   Tem1   Tem2    mass_displ    mass_accel    ground_displ    ground_accel");
	//fprintf(fvibr,"\n#  s         K      K         mm           m/s**2          mm             m/s**2  ");
	//
	//fmagforce = fopen("add_mag_force.dat","w");
	//fprintf(fmagforce  ,"#tim_etap  Tem1  Tem2  mag_force_factor ");
	//fprintf(fmagforce,"\n#  s        K     K                     ");
	
	if(fkeyread_double_in_section(fwzd, cur_section,"htim_preset", &htim_preset) ) htim_preset=0.0001;
	if(fkeyread_double_in_section(fwzd, cur_section,"hmf_factor_preset", &hmf_factor_preset) ) hmf_factor_preset=0.0001;

	//if(cur_etap.num==5){ 
	//	printf("\nExper 713"); 
	//	pause();
	//}

	if(cur_task == vibr_trig_polynom){
	//\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
     //VVVVVVV  Begin reading and printing the type and parameters of the exitation force
       //reading coefficients and frequencies of excitation from cur_section of the wzd-file
		fkeyread_int_in_section(fwzd, cur_section,"Number_of_terms",&k_harmonic);

		ampl = new Real [k_harmonic];
		ffindkey_in_section(fwzd ,cur_section ,"Excitation_displacement_amplitudes");

		for(i_harmc=0;i_harmc<k_harmonic;i_harmc++){
			fscanf(fwzd,"%lg",&Dx); ampl[i_harmc] = Dx*1e-3; //mm -> m
		}

		angular_frequency = new Real [k_harmonic];
		ffindkey_in_section(fwzd ,cur_section ,"Frequencies_of_harmonics");

		for(i_harmc=0;i_harmc<k_harmonic;i_harmc++){
			fscanf(fwzd,"%lg",&Dx); angular_frequency[i_harmc]=Dx*two_pi; //Hz -> rad/s
		}

		printf("\n htim_preset==%lg", htim_preset);
		printf("\n k_harmonic=%d", k_harmonic);
		printf("\n excitation amplitudes(mm)=");
		for(i_harmc=0;i_harmc<k_harmonic;i_harmc++){
			printf("\t%8.3lg mm",ampl[i_harmc]*1e3); //m -> mm
		}
		printf("\n excitation_frequencies(Hz)=");
		for(i_harmc=0;i_harmc<k_harmonic;i_harmc++){
			printf("\t%8.3lg Hz",angular_frequency[i_harmc]/two_pi); // rad/s -> Hz
		}
	//	pause();
     //^^^^^^^  End of reading and printing the type and parameters of the exitation
     ////////////////////////////////////////////////////////////////////////////////////////
	}//endif(cur_task == vibr_trig_polynom)

	//Assign initial values before opening cycle on steps:

	//definition: operation_flag = vibr_trig_polynom_flag || idle_flag;
	if(operation_flag==1 && operation_flag_prev==1){}//do nothing
	else{ //i.e. if prev_task==preload || prev_task==add_magnetic_force 
		tim_etap=0.0;// time measured from the beginning of the vibrations
		ground_velocity=0.0; //velocity on the ground
		hground_velocity=0.0; //increm. of velocity on the ground
		ground_displ=0.0;   //displacement on the ground
		hground_displ=0.0;   //increm. of ground_displ
	}//only for the 1st operational etap we assign initial zero values

	
	htim = htim_preset; //Default is htim_preset = 0.0001s
	if(add_magnetic_force_flag){
		porous.hmf_factor = hmf_factor_preset;
		tim_z = htim_preset / hmf_factor_preset;//we set the duration of this etap such as needed to reach porous.hmf_factor==1
		tim_ke = tim + tim_z; 
		printf("\nExper729 Adding magnetic force is on");
	}
	else porous.hmf_factor = 0.0;

	if(idle_flag){
		if(k_step_flag) htim = tim_z / (Real)k_step;
	}

	printf("\nExper739 tim=%10.4lf tim_z=%10.5lf "
			   ,tim		   ,tim_z   );
	printf("\n         tem=%10.4lf tem_z=%10.5lf"
						,porous.tem_old  ,tem_z);
	//printf("\ntau1=%10.4lf tau2=%10.5lf",sig1(2,1)*1e-6 ,sig2(2,1)*1e-6);
	//pause(); //pause(pause_factor);

	Real acc0=0.0, acc1=0.0;
	//Main cycle on steps
	for(istep=0; tim<=tim_ke && *macro_failure_pointer==0; istep++){
		//cycle will be breaked also when file accel.d comes to end
		//step adjusting is not realized as htim may be specified in file accel.d 
		Real tim_etap_new = tim_etap + htim;
		acc0 = Excitation_acceleration(cur_task,tim_etap, k_harmonic, ampl, angular_frequency);
		acc1 = Excitation_acceleration(cur_task,tim_etap_new, k_harmonic, ampl, angular_frequency);
		hground_displ = Excitation_displacement(cur_task,tim_etap_new, k_harmonic, ampl, angular_frequency)
							-Excitation_displacement(cur_task,tim_etap, k_harmonic, ampl, angular_frequency);
		hground_velocity= Excitation_velocity(cur_task,tim_etap_new, k_harmonic, ampl, angular_frequency)
							-Excitation_velocity(cur_task,tim_etap, k_harmonic, ampl, angular_frequency);
		
		//if(cur_etap.num==5){
		//	printf("\n\nExper737 tim_etap=%lg hground_displ=%lg hground_velocity=%lg",tim_etap,hground_displ,hground_velocity);
		//	pause();
		//}
		
		   //------------------------------------------------------------------
	   //begin finding htim and htem
		f_step= htim / tim_z;
		porous.htem = tem_z * f_step;
		//porous.htem2 = tem2_z * f_step;
		porous.tem_new = porous.tem_old + porous.htem;
		//porous.tem2_new = porous.tem2_old + porous.htem2;
		porous.mf_factor_new = porous.mf_factor_old + porous.hmf_factor;
		//end finding htim and htem
		if(istep==0 && operation_flag_prev==0) VibrContinue=0; else VibrContinue=1;
						
		if( (retcode = 
						Task_force_tem
						( 
						istep, htim, tem, sig, hsig, sig_elastic, hsig_elastic, eps, heps, MC, porous, Vrnts, GrPar
						)
			)!=0 

			) return retcode;

	//Update parameters:
	ground_velocity+=hground_velocity;
	ground_displ+=hground_displ;
	tim += htim; tim_etap+=htim; 	//TaskVibr32_Isolation has calculated hmass_displ, hmass_velocity, heps and hsigma32
	porous.Update();
	mass_accel = porous.hvelocity / htim;

	//tim_prev=tim; sig_prev=sig; eps_prev=eps; heps_prev=heps;

	//printf("\nE239: istep=%i mass_displ=%lg force1=%lg  force2=%lg"
	//				,istep  ,mass_displ*1e3 ,force1*1e2 ,force2*1e2);//force N->g

	if(vibr_trig_polynom_flag){
		//fprintf(fvibr  ,"#tim_etap   Tem1   Tem2    mass_displ    mass_accel    ground_displ    ground_accel");
		//fprintf(fvibr,"\n#  s         K      K         mm           m/s**2          mm             m/s**2  ");
		fprintf(fvibr,"\n %e   %5.1lf  %5.1lf    %+e  %+e  %+6.3lf  %+6.3lf"
			,tim_etap , porous.tem_new , porous.displ_new*1e3 ,mass_accel ,ground_displ*1e3 ,acc1 , porous.force_new);
	}
	//if(add_magnetic_force_flag){
	//	//fprintf(fmagforce  ,"#tim_etap  Tem1  Tem2  mag_force_factor mass_displ   mass_accel");
	//	//fprintf(fmagforce,"\n#  s        K     K                        mm          m/s**2  ");
	//	fprintf(fmagforce,"\n  %e  %7.3lf   %6.4lf  %e  %e"
	//		,tim_etap ,porous.tem_new ,porous.mf_factor_new ,porous.displ_new*1e3 ,mass_accel);
	//}

{ //block replacing if(...)
	//if(!k_print || tim>=tim_print){//print results
	//	tim_print+=htim_print;//find new time, when to print results
		//sig1=0.0; sig1(2,2) = c_beam1[iuz_max].K_sig * force1;
		//sig2=0.0; sig2(2,2) = c_beam2[iuz_max].K_sig * force2;
		int na_ekran=1;
		print_step_2slotted_elems(fttesf ,print_step_flag ,na_ekran ,istep ,tim ,porous.tem_old  
										,porous.force_old ,porous.kc_beams ,porous.c_beam ,porous.displ_new_1 ,porous.displ_new_2, porous.displ_new_3, porous.displ_new_4, porous.e_macro ,porous.stress_macro ,porous.force_new_eff);
}//end block							//old values here are the same as new
	
	//pause();

  }//for istep

if (faccel!=NULL) fclose(faccel);
if (faccel!=NULL) fclose(faccel);

delete[] ampl;
delete[] angular_frequency;
//if (tem_peaks!=NULL) delete[] tem_peaks;
printf("\n\nExper 791: case operation ended");
//if(add_magnetic_force_flag==1){
//	Real Fmag=porous.F_mag(porous.displ_old,porous.mf_factor_old);
//	printf("\nExper793: Magnetic force has been added");
//	printf("\nExper795 Magnetic force =%lgN displ=%lgmm",Fmag ,porous.displ_old*1e-3);
//}

}//end of unconditional block surrounding case operation
	break;


/*	case vibr_2_elements: //SMA-helix as isolator from vibrations
      //Acceleration of the ground: accel =SUM_i{ accel_ampl[i]*cos((2*pi/accel_period[i])*tim)}
      // k_vibr32_harmonic - number of terms (harmonics)
  {//begin of unconditional block surrounding this case
	//Majatnik Maj1;     
	
	int i_harmc;
	printf("\n\nExper: case vibr_2_elements");
	printf("\ntim=%10.4lf tim_z=%10.5lf tem1=%10.4lf tem1_z=%10.5lf tem2=%10.4lf tem2_z=%10.5lf"
			,tim		,tim_z			,tem1		,tem1_z			,tem2		,tem2_z);
	printf("\ntau1=%10.4lf tau2=%10.5lf",sig1(2,1)*1e-6 ,sig2(2,1)*1e-6);
	pause(0); //pause(pause_factor);
	  
	  //Maj1.Init();
	//open file for registering vibrograms:
	fvibr = fopen("vibr.dat","w");
	//fprintf(fvibr,"#    time   temperature ground_displ mass_displ ground_accel beam_accel");
	fprintf(fvibr  ,"#tim_etap      temperature mass_displ    mass_accel    ground_displ    ground_accel");
	fprintf(fvibr,"\n#s                  K         mm           m/s**2          mm             m/s**2  ");
     //-------------------------------------------------------------------------------------
     //VVVVVVV  Begin reading and printing the type and parameters of the exitation32
     if(strcmp(excitation, "trig_polynom")==0){
       //open file for reading coefficients and frequencies of quake excitation
       if( (filrd = fopen("excit_trig_poly.d","r")) == 0 ) Exit_error("excit_trig_poly.d not found\n");
       fkeyread_int(filrd, "Number_of_terms", &k_harmonic);

       Real Dx; // temp. variable for reading
		ampl = new Real [k_harmonic];
        ffindkey(filrd, "Excitation_displacement_amplitudes");
       for(i_harmc=0;i_harmc<k_harmonic;i_harmc++){
         fscanf(filrd,"%lg",&Dx); ampl[i_harmc] = Dx*1e-3; //mm -> m
       }

       angular_frequency = new Real [k_harmonic];
       ffindkey(filrd, "Frequencies_of_harmonics");
       for(i_harmc=0;i_harmc<k_harmonic;i_harmc++){
         fscanf(filrd,"%lg",&Dx); angular_frequency[i_harmc]=Dx*two_pi; //Hz -> rad/s
       }

       printf("\n\n k_harmonic=%d", k_harmonic);
       printf("\n excitation amplitudes(mm)=");
       for(i_harmc=0;i_harmc<k_harmonic;i_harmc++){
         printf("\t%8.3lg mm",ampl[i_harmc]*1e3); //m -> mm
       }
       printf("\n excitation_frequencies(Hz)=");
       for(i_harmc=0;i_harmc<k_harmonic;i_harmc++){
         printf("\t%8.3lg Hz",angular_frequency[i_harmc]/two_pi); // rad/s -> Hz
       }
	   pause(0);
       fclose(filrd);
	 }// if excitation == "trig_polynom"
     //open file accel.d 
     else if(strcmp(excitation, "quake_record")==0){
       if( !(faccel = fopen("accel.d","r")) ) Exit_error("file accel.d not found\n");
     }
	 else Exit_error("\n\nInvalid type of excitation");
     //^^^^^^^  End of reading and printing the type and parameters of the exitation32
     //--------------------------------------------------------------------------------------
     ////////////////////////////////////////////////////////////////////////////////////////
	 //VVVVVVV      |||| Specification of control ||||
     //             VVVV                          VVVV
 
     // Inherited from previous etap are:
     // tim_rev, tim_rev_prev, tim_period, tim_wait_s,
     // eps_peak, eps_peak_prev (initially all of them = MAXDOUBLE)
	 if(!control){ 
		printf("\nPassive control of vibration32 with isolation; end time=%lg\n ", tim+tim_z);
		pause(pause_factor);
	 }
     else{
		Exit_error("\n\nvibration32 with control temporarily are not envisaged ");
		//printf("\nvibration32 with control impulse_shape= %s; end time=%e ",impulse_shape, tim+tim_z);
		 // Declare and enter pendulum parameters
		 //if( !(tem_peaks = new Real [heatnum]) ) Exit_error("Cannot allocate array tem_peaks\n");
		//if(strcmp(impulse_shape,"exp1")==0) printf("\ncool_exp_A=%7.3lf",cool_exp_A);
		//if(strcmp(impulse_shape,"exp2")==0 || strcmp(impulse_shape,"calculated")==0)
		//	printf("\ncool_exp_A=%7.3lf cool_exp_M=%7.3lf",cool_exp_A ,cool_exp_M);
		//	if(strcmp(impulse_shape,"exp3")==0)
		//		printf("\ncool_exp_A=%7.3lf cool_exp_M=%7.3lf cool_exp_AM=%7.3lf"
		//				,cool_exp_A ,cool_exp_M ,cool_exp_AM);
		//	printf("\n");
		//	if(strcmp(wait_mode,"seconds")==0) printf("tim_wait=%7.3lf seconds",tim_wait);
		//	else printf("tim_wait=%7.3lf part of previous period",tim_wait);

		//	printf("\ntim_heat=%7.3lf tem_fon=%5.1lf" ,tim_heat ,tem_fon);

		//	if(strcmp(tem_peaks_mode,"equal")==0) {for(int i=0;i<heatnum;i++) tem_peaks[i]=tem_peak; }
		//	else {
		//	  ftem_peaks = fopen(tem_peaks_mode,"r");
		//	  if(ftem_peaks==NULL) Exit_error("File %s not found", tem_peaks_mode);
		//	  int heatnum1; fkeyread_int(ftem_peaks, "heatnum", &heatnum1); //number of heating impulses
		//	  if(heatnum1!=heatnum) Exit_error("\nheatnum in file %s is not equal to heatnum in wzd-file"
		//										, tem_peaks_mode);
		//	  //heatnum=heatnum1;// altering the value of heatnum
		//	  ffindkey(ftem_peaks,"tem_peak");
		//	  for(int i=0;i<heatnum;i++) fscanf(ftem_peaks,"%lg", &tem_peaks[i]);
		//	  fclose(ftem_peaks);
		//	}//end of else - when temperature peaks are different
		//	{ int i;
		//	  printf("\nHeating impulses:\n");
		//	  for(i=0;i<heatnum;i++) printf("%5d\t",i);
		//	  printf("\n");
		//	  for(i=0;i<heatnum;i++) printf("%5.1lf\t",tem_peaks[i]);
		//	}//unconditional block
	 }//else i.e. if control==1

     //        ^^^^                          ^^^^
     //        |||| Specification of control ||||

     //Assign initial values before opening cycle on steps:
	Real tim_etap=0.0; // time measured from the beginning of the stage (etap)
	ground_velocity=0.0; //velocity on the ground
	hground_velocity=0.0; //increm. of velocity on the ground
	ground_displ=0.0;   //displacement on the ground
	hground_displ=0.0;   //increm. of ground_displ
	//mass_displ = element_displ - ground_displ;
	mass_displ = 0.0;
	mass_velocity = 0.0;
	mass_accel = 0.0;


	htim = htim_preset; //Default for htim_preset is 0.0001s
     //tim_wait_s=1e100;   //very big value, such that condition tim>=tim_rev+tim_wait_s == FALSE
     //tim_rev_prev=1e100; //then no heating could be prior to 2 semiperiods
     //tim_rev=1e100;
     //tim_prev=tim;//to avoid warning for assignment tim_rev=(tim_prev+tim)/2.0;
     //dtem_heat = MAXDOUBLE;
     //semiperiod=0;
     //iheat=0;
     //heps_prev=0.0;//with this value condition of revers (heps(2,1)<0.0)&&(heps_prev(2,1)>0.0)==false
     //heating=0;//this flag==1 when heating is on (heating may not be on the first step)

     Real acc0=0.0, acc1=0.0;
     //Main cycle on steps
	 for(istep=0; tim<=tim_ke && *macro_failure_pointer==0; istep++){
     //cycle will be breaked also when file accel.d comes to end
     //step adjusting is not realized as htim may be specified in file accel.d 
	Real tim_etap_new = tim_etap + htim;
	if(strcmp(excitation, "trig_polynom")==0){
		acc0 = Excitation_acceleration(excitation,tim_etap, k_harmonic, ampl, angular_frequency);
		acc1 = Excitation_acceleration(excitation,tim_etap_new, k_harmonic, ampl, angular_frequency);
		hground_displ = Excitation_displacement(excitation,tim_etap_new, k_harmonic, ampl, angular_frequency)
							-Excitation_displacement(excitation,tim_etap, k_harmonic, ampl, angular_frequency);
		hground_velocity= Excitation_velocity(excitation,tim_etap_new, k_harmonic, ampl, angular_frequency)
							-Excitation_velocity(excitation,tim_etap, k_harmonic, ampl, angular_frequency);
	}
	else if(strcmp(excitation, "quake_record")==0){
		acc0 = acc1;
		Real tim_x;
		if(fscanf(faccel,"%lg %lg", &tim_x, &acc1)==EOF) break;//from cycle on istep
		htim = tim_x - tim;
		if(htim<=0.0) Exit_error("invalid time in file accel.d");
		//this scheme is similar to Euler-Cauchy method for eqn. y'' = a(t)
		hground_velocity = htim * (acc0+acc1)/2.0;
		hground_displ = htim*(ground_velocity+hground_velocity/2.0);
	}
       //------------------------------------------------------------------
       //start finding htim and htem
	if(control){
		//if(heating){
		//	if(dtem_heat==MAXDOUBLE) Exit_error("Impossible: dtem_heat == MAXDOUBLE");
		//	htem = tem_incr; htim = htem/dtem_heat; //heating can be only after cooling so that 
		//	if(tem>=tem_peak) heating=0;             //dtem_heat is known
		//}
		//else{//i.e. if heating==0
		//	dtem_heat = (tem_peaks[iheat]-tem)/tim_heat;
		//	if(iheat<heatnum && tim>=tim_rev+tim_wait_s ){ //start heating
		//		iheat++; htem = tem_incr; htim = htem/dtem_heat;
		//		heating=1;
		//	}
		//	if(strcmp(impulse_shape,"calculated")==0 || strcmp(impulse_shape,"exp1")==0) tim_cool = cool_exp_A;
		//	if(strcmp(impulse_shape,"exp2")==0) tim_cool = XXold.Phase < 0.0001 ? cool_exp_A : cool_exp_M;
		//	if(strcmp(impulse_shape,"exp3")==0)
		//	if(XXold.Phase < 0.0001) tim_cool = cool_exp_A;
		//	else if(XXold.Phase > 0.9999) tim_cool = cool_exp_M;
		//	else                          tim_cool = cool_exp_AM;
		//	if(tem>tem_fon){
		//		Real tem_incr_at_cooling = XXold.dPhase>1e-12 ? -0.1 : -1.0;
		//		Real dT=(tem_fon-tem)/tim_cool;

		//		if(strcmp(impulse_shape,"calculated")==0 && XXold.dPhase>0.0)
		//		//dT/=1+(XXnew.dPhase*htim*MC.Q0/MC.SpecificHeatM/MC.DensityM/htem);
		//		dT += q0_C_Ro * XXold.dPhase;

		//		htim=fabs(tem_incr_at_cooling/dT);
		//		if(htim > htim_preset) htim=htim_preset;

		//		htem=dT*htim;
		//		if(tem+htem<=tem_fon) htem=tem_fon-tem;
		//	}//if tem>tem_fon
  //         else{htem=0.0; htim=htim_preset;}
		//}//else (if heating==0)
	}//if control
	else{//i.e. if control==0
		htim = htim_preset;
		f_step= htim / tim_z;
		htem1 = tem1_z * f_step;
		htem2 = tem2_z * f_step;
	}// end of else (case when control==0)
	//end finding htim and htem
	if(istep==0 && prev_task!=vibr32_quake) VibrContinue=0; else VibrContinue=1;
	//Matr33 sigma, hsigma;
	//f_step = htim / tim_z;
	//hsig(2,2) = s33_z * f_step;
	//press = press_z * f_step;
	//hfluence = fluence_z * f_step;
	//hs33_us= ultrasound33_flag ? siga_us33*sin(omega_us*(tim_etap+htim)) - s33_us : 0.0;
	//sigma = sig;                      hsigma = hsig;
	//sigma.Add_to_diagonal(-press);    hsigma.Add_to_diagonal(-hpress);
	//if(ultrasound33_flag){
	//	hs33_us=siga_us33*sin(omega_us*(tim_etap+htim)) - s33_us; 
	//	sigma(2,2)+=s33_us; hsigma(2,2)+=hs33_us; 
	//}
	if( (retcode =
				TaskVibr32_2elems(VibrContinue	,htim	,tem	,htem	,fluence	,hfluence
							,sig1	,hsig1	,sig2	,hsig2
							,hground_displ	,hground_velocity	
							,mass_displ	,hmass_displ
							,mass_velocity ,hmass_velocity
							,MC ,Vrnts ,GrPar ,Helix
							,XXold1 ,XXnew1 ,XXold2 ,XXnew2 )

		)!=0 ) {return retcode;}

//>>>>>>>>>>>>>>>>
	//if revers
	//if( ((heps(2,1)<0.0)&&(heps_prev(2,1)>0.0)) ||   //if revers
	//		((heps(2,1)>0.0)&&(heps_prev(2,1)<0.0)) ) {
	//	if(tim_rev<MAXDOUBLE) tim_rev_prev=tim_rev;
	//	tim_rev=(tim_prev+tim)/2.0;
	//	if(tim_rev_prev<MAXDOUBLE){
	//		tim_period = 2.0*(tim_rev-tim_rev_prev); 
	//		if(strcmp(wait_mode,"fractions")==0) tim_wait_s = tim_period * tim_wait;
	//		else tim_wait_s=tim_wait;
	//	}
	//	if(eps_peak<MAXDOUBLE) eps_peak_prev=eps_peak;
	//	eps_peak=(eps_prev(2,1)+eps(2,1))/2.0;
	//	if(eps_peak_prev<MAXDOUBLE){
	//		if(amplit<MAXDOUBLE) amplit_prev=amplit;
	//		amplit=fabs(eps_peak-eps_peak_prev);
	//		//decrem= semiperiod>2 ? log(amplit_prev_prev/amplit) : 0.0;
	//		if(amplit_prev<MAXDOUBLE){
	//			decrem= 2.0*log(amplit_prev/amplit);
	//			fprintf(fdec,"%5u %12.6lf %12.6lf %12.6lf\n",semiperiod, tim_prev, amplit*2e2, decrem);
	//		}
	//	}
	//	semiperiod++;
	//}// if revers

	//Update parameters:
	ground_velocity+=hground_velocity;
	ground_displ+=hground_displ;
	//TaskVibr32_Isolation has calculated hmass_displ, hmass_velocity, heps and hsigma32
	mass_displ += hmass_displ;
	mass_velocity += hmass_velocity;
	mass_accel = hmass_velocity / htim;

	//tim_prev=tim; sig_prev=sig; eps_prev=eps; heps_prev=heps;
	tim += htim; tim_etap+=htim; tem += htem; sig1 += hsig1; sig2 += hsig2; 
	XXold1 = XXnew1; XXold2 = XXnew2;

	//Print vibrogram for each step
	//fprintf(fvibr  ,"#  tim_etap      temperature mass_displ    mass_accel    ground_displ    ground_accel");
	//fprintf(fvibr,"\n#     s               K          mm           m/s**2          mm             m/s**2  ");
	fprintf(fvibr,"\n %e   %8.3lf  %e  %e  %e  %e"
					,tim_etap ,tem ,mass_displ*1e3 ,mass_accel ,ground_displ*1e3 ,acc1);

	Matr33 eps1=XXnew1.TotalStrain;
	if(!k_print || tim>=tim_print){//print results
		tim_print+=htim_print;//find new time, when to print results
		print_step2elems(fttesf ,print_step_flag ,na_ekran, istep, tim, tem, sig1, sig2, XXnew1, XXnew2);
	}//if() print results
//<<<<<<<<<<<<<<<<
  }//for istep

if (faccel!=NULL) fclose(faccel);

delete[] ampl;
delete[] angular_frequency;
//if (tem_peaks!=NULL) delete[] tem_peaks;

}//end of unconditional block surrounding case vibr32_quake
	break;
*/
	default:
		Exit_error("899: unknown mode of boundary condition");
		break;
     //this case can never happen unless some errror with memory distribution
  }//switch
//-----------------------------------------------------------------------  
//--------------- end of switch of cases --------------------------------
//-----------------------------------------------------------------------
  
  int sufficient_print_condition=*macro_failure_pointer; //print to "cycle.dat" if macro failure occurs
  for(ifil=0; ifil<k_cycle_files ; ifil++){
    if(etap_print_numbers[ifil]==cur_etap.num)
      print_etap(fcycle[ifil],sufficient_print_condition, cyc[0] ,tim ,tem,sig ,eps ,XXnew);
  }
  printf("\n\nEnded etap %3i: cur_task#%i", cur_etap.num ,cur_task);
  printf("\ttim=%lg tem=%8.3lf ",tim, tem);
  pause(pause_factor);

 }//while fread_next_section && XXnew.PI_macro==0
//----------------------------------------------------------------------- 

 if(*macro_failure_pointer){
   printf("\nMacro-fracture occurred. %3i of %3i grains destroyed", XXnew.kzer_dmg, MC.kzer);
 }
//if(ampl!=NULL) delete[] ampl;
//if(angular_frequency!=NULL) delete[] angular_frequency;
//if(tem_peaks!=NULL) delete[] tem_peaks;
fclose(fvibr);
//fclose(fmagforce);

_fcloseall();
// gettime(&t);
// printf("\nEnd of the experiment in %2d:%02d:%02d\n", t.ti_hour ,t.ti_min ,t.ti_sec);
 printf("\nEnd of the experiment");
 
 pause(); 
 
 exit(0);
 return 0;
} /* main */