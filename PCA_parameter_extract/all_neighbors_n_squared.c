/* 6.10.19 
 * having some interface problems with this code.  
 * lets dig in and see what is going wrong
 *
 * 1.23.19
 *
 * ...eventually get this to qsort() your data...
 *
 *making an n^2 algorithm to find all neighbor distances.
 *using cluster code since it is good at reading in smtest data
 *
 *
 * Revisions log:
 * 1.06.19 printing out r_nn, x_nn, y_nn and calculating mean, <r_nn>
 *          \delta r_nn' = r_nn - <r_nn>
 *	   sorting from small to large(?)
 *
 * 7.22.17 use cell N log N algorithm to find near neighbors, 
 *          measure the distances, 
 *          look at near neighbors and angles
 *          also consider \psi_6
 *
 *	   also clean up the code to remove the older variables
 *	   that are no longer used
 *
 * 5.13.16 default input name "smtest" 
 */

#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<string.h>
#include <sys/stat.h>

#define PI 3.14159265359
#define DEBUG 0
#define VERBOSE 0

FILE *f_debug;

struct vortex{
  int id;
  int color;
  double x;
  double y;
  double radius;

  //DM measurements
  int neighbor_count;
  //double neighbor_distance_list[100];


  /*trying to emulate this:
    for(j=0;j<lookupdata.ncellsy;j++){
    for(i=0;i<lookupdata.ncellsx;i++){
    lookuptable[i][j].num=0;
    lookuptable[i][j].list=malloc(maxcell*sizeof(int));
    }
    }

   */
  
};

struct syssize{
  double SX;
  double SY;
  double SX2;
  double SY2;
};

struct mylookup{
  int num;
  int *list;
};

struct lookupdata{
  int ncellsx;
  int ncellsy;
  double xscale;
  double yscale;
};

struct parameters{

  int nV;
  int nP;
  //int nV1, nV2;
  int maxnum;
  
  //int maxnum; //_small;
  //int maxnum_large;
  
  double radius; //_small;
  //double radius_large;
  double runforce;
  int runtime;
  //
  double density;
  //double phi2_small;
  //double phi1_big;
  double pdensity;
  //
  double dt;
  int maxtime;
  int starttime;           //useful for restarting simulation from config file
  int writemovietime;
  double potential_radius; //change to pin_radius
  double potential_mag;    //change to pin_max_force;
  double kspring;
  double drive_mag;
  double drive_frq;
  int decifactor;

  //change the driving force at a given time by a given amount.
  int drive_step_time;
  double drive_step_force;

  int restart;
};



void distance(double *dr, double *dx, double *dy,
	      double x1, double y1,
	      double x2, double y2,
	      struct syssize syssize);

//DM read in Pa0 
void get_parameters_file(struct parameters *parameters,
			 struct syssize *syssize);

void pair_interact(int id1,int id2,struct vortex *vortex,
		   struct syssize syssize,
		   struct parameters parameters,
		   FILE *out);




main(int argc,char *argv[])
{
  //file names/pointers
  FILE *in,*out; //

  char ascii_file[120]; //="velocity_data/XV_data_t=";
  char str_time[10];

  //if directory doesn't exist, make directory
  //make and/or open directory snapshot_files
  struct stat st = {0};       //null directory
  if (stat("neighbor_data", &st) == -1) {
    mkdir("neighbor_data",0700);
  }
  
  //parameters
  double density;
  struct syssize syssize;
  int runtime;
  int maxtime;

  //contains everything from Pa0, easy to pass around
  struct parameters parameters;
  
  //
  int completeframe;
  struct vortex *vortex;
  int maxnum;
  
  char filename[120]="smtest";
  int nV, time, nVold;
  int i,j,k;

  int frame;
  
  //grab everything from Pa0
  get_parameters_file(&parameters,&syssize);

  //reassign some variables locally - mainly for historic reasons
  maxnum = parameters.maxnum;
  maxtime = parameters.maxtime;
  runtime = parameters.runtime;
  density = parameters.density;


  //more system specific.  vortex struct only contains maxnum
  //and we create room for a neighbor list at the end.
  vortex=malloc(maxnum*sizeof(struct vortex));
  // + (maxnum+1)*sizeof(double) );

  //make an array to calculate all of the averages
  //within every frame.
  //so every frame we zero all the values
  //and then do a running total within the frame
  double rN_average[maxnum];

  //add the information from smtest
  if(argc<2){
    // printf("Enter name of movie: ");
    // scanf("%s",filename);
  }
  else
    sscanf(argv[1],"%s",filename);

  if((in=fopen(filename,"r"))==NULL){
    printf("Error opening file %s\n",filename);
    exit(-1);
  }

  //create a list of probe particle centered neighbor distances
  //each particle has its own list of neighbors,
  //sorted from small to large
  //and each value is normalized across particles
  //such that nn=1 is averaged and subtracted, likewise nn=2, etc
  
  // out=fopen("delta_rN.dat","w");

  if(DEBUG){
    f_debug = fopen("debug.dat","w");
  }
  
  //TODO-----------------------------------------
  //total movie stats - update for new measures
  frame=0;

  //for each new frame, the average neighbor distances are reset to zero
  for(i=0;i<nV-1;i++){
    rN_average[i] = 0.0;
  }
  
  //read in all frames of smtest and analyze individually
  while(!feof(in)){
    frame++;
    nVold=nV;
    completeframe=read_frame(in,&nV,&time,vortex);

    //cut early frames
    if(!completeframe) break;
    if(frame==0) continue;

    //open file to write data
    strcpy(ascii_file,"neighbor_data/frame_");
    sprintf(str_time,"%08d",time); //convert current to a string   
    strcat(ascii_file,str_time);
    if(DEBUG) printf("printing frame data to %s\n",ascii_file);

    if((out=fopen(ascii_file,"w"))==NULL){
      printf("Error opening file %s\n",ascii_file);
      exit(-1);
    }
    //------------------------------------

    for(i=0;i<nV;i++){
      fprintf(out,"#frame %d, particle %d\n",frame,i);
      for(j=0;j<nV;j++){

	if(i!=j){
	  pair_interact(i,j,vortex,syssize,parameters,out);
	}
	
      }      
    }
    
    
    fclose(out);

  }//do the same for a new frame (i.e. end while(feof) loop)

  

  //--------------------------------------------------------

  //  fclose(out);

  if(DEBUG){
    fclose(f_debug);
  }

}

//------------------------------------------------------------------
int read_frame(FILE *in,int *nV,int *time,struct vortex *vortex)
{
  int i;
  float xin,yin,zin;

  fread(nV,sizeof(int),1,in);
  fread(time,sizeof(int),1,in);
  if(feof(in)) return 0;
  for(i=0;i<*nV;i++){
    fread(&(vortex[i].color),sizeof(int),1,in);
    fread(&(vortex[i].id),sizeof(int),1,in);
    fread(&xin,sizeof(float),1,in);
    fread(&yin,sizeof(float),1,in);
    fread(&zin,sizeof(float),1,in);
    vortex[i].x=(double)xin;
    vortex[i].y=(double)yin;
    vortex[i].radius=(double)zin;

    vortex[i].neighbor_count = 0;

    //
    
  }
  return 1;
  
}

//---------------------------------------------------------
//---------------------------------------------------------
void pair_interact(int id1,int id2,
		   struct vortex *vortex,
		   struct syssize syssize,
		   struct parameters parameters,
		   FILE *out)
{
  double dr,dx,dy;

  if(DEBUG){
    fprintf(f_debug,"%d %d \n",id1,id2);
    fflush(stdout);

    fprintf(f_debug,"%lf %lf \n",vortex[id1].x,vortex[id2].x);
    fprintf(f_debug,"%lf %lf \n",vortex[id1].y,vortex[id2].y);
    fflush(stdout);
  }
  
  //center to center distance with PBC, returns to dr
  distance(&dr,&dx,&dy,
	   vortex[id1].x,vortex[id1].y,
	   vortex[id2].x,vortex[id2].y,
	   syssize);

  //just print dr since that is the only analysis we are using.  save data.
  //fprintf(out,"%d\t %d\t %f\t %f\t %f \n",id1, id2, dr,dx,dy);    
  fprintf(out,"%f\n",dr);    

  //Assign the measured distance to the neighbor arrays in each vortex struct.
  //vortex[id1].neighbor_distance_list[vortex[id1].neighbor_count] = dr;
  //vortex[id1].neighbor_count++;
  
  //if(vortex[id1].neighbor_count > parameters.maxnum){
  //  printf("Exceeded max number of total distances!\n");
  //  fflush(stdout);
  //  exit(-1);
  //}

  //vortex[id2].neighbor_distance_list[vortex[id2].neighbor_count] = dr;
  //vortex[id2].neighbor_count++;

  //if(vortex[id2].neighbor_count > parameters.maxnum){
  //  printf("Exceeded max number of total distances!\n");
  //  fflush(stdout);	
  //  exit(-1);
  //}
	

  return;
  
}//end subroutine.  

//-------------------------------------------------------------


//--------------------------------------------------------------------
//-----------------------------------------------------------------
void get_parameters_file(struct parameters *parameters,
			 struct syssize *syssize)
{
  FILE *in;
  char trash[120];
  double cellsize,length_scale;
  double resolution;

  
  resolution=1e-6;

  if((in=fopen("Pa0","r"))==NULL){
    printf("Input file Pa0 not found\n");
    exit(-1);
  }
  else if (VERBOSE) printf("reading Pa0.\n");
  fflush(stdout);

  fscanf(in,"%s %lf\n",trash,&((*parameters).density));
  //fscanf(in,"%s %lf\n",trash,&((*parameters).phi2_small));
  //fscanf(in,"%s %lf\n",trash,&((*parameters).phi1_big));
  fscanf(in,"%s %lf\n",trash,&((*parameters).pdensity));

  //----------------------------------------------------------------------
  //check whether the user gave you sensible choices for the two densities
  //----------------------------------------------------------------------
  //double test1 = (*parameters).phi2_small + (*parameters).phi1_big;
  //if( fabs( ((*parameters).density-test1) >0.01) ){
  //  printf("Inconsistency in requested densities.\n");
  //  printf("Total density: %f\n", (*parameters).density);
  //  printf("Small disk density: %f\n", (*parameters).phi2_small);
  //  printf("Large disk density: %f\n", (*parameters).phi1_big);
  //  exit(-1);
  //}
  //---------------------------------------------------------------------
  
  if (VERBOSE) printf("reading parameters density %f, pin density %f.\n", (*parameters).density, (*parameters).pdensity);
  fflush(stdout);

  //-
  fscanf(in,"%s %lf\n",trash,&((*syssize).SX));
  (*syssize).SX2=(*syssize).SX*0.5;
  //--
  fscanf(in,"%s %lf\n",trash,&((*syssize).SY));
  (*syssize).SY2=(*syssize).SY*0.5;
  
  fscanf(in,"%s %lf\n",trash,&((*parameters).radius)); //_small));
  //fscanf(in,"%s %lf\n",trash,&((*parameters).radius_large));

  //--where it is appropriate to accurately calculate maxnum
  //--based on particle sizes

  //double prefix = (*syssize).SX*(*syssize).SY/PI;
  //double phiS_term = (*parameters).phi2_small/( (*parameters).radius_small*(*parameters).radius_small) ;
  //double phiB_term = (*parameters).phi1_big/( (*parameters).radius_large*(*parameters).radius_large);
   
  //(*parameters).maxnum= (int) floor(prefix * (phiS_term + phiB_term));
  (*parameters).maxnum= (int) floor( (*syssize).SX*(*syssize).SY/PI  *  (*parameters).density /( (*parameters).radius*(*parameters).radius) );
  
  if(VERBOSE) printf("\n maxnum is: %d\n", (*parameters).maxnum);
  //--
 
  fscanf(in,"%s %d\n",trash,&((*parameters).runtime));
  fscanf(in,"%s %lf\n",trash,&((*parameters).runforce));
  fscanf(in,"%s %lf\n",trash,&((*parameters).dt));
  fscanf(in,"%s %d\n",trash,&((*parameters).maxtime));
  if(VERBOSE) printf("\n maxtime is: %d\n", (*parameters).maxtime);

    
  fscanf(in,"%s %d\n",trash,&((*parameters).writemovietime));
  fscanf(in,"%s %lf\n",trash,&((*parameters).kspring));
  
  fscanf(in,"%s %lf\n",trash,&cellsize); // size of lookup cell
  fscanf(in,"%s %lf\n",trash,&((*parameters).potential_radius));
  fscanf(in,"%s %lf\n",trash,&((*parameters).potential_mag));
  
  fscanf(in,"%s %lf\n",trash,&length_scale);
  fscanf(in,"%s %lf\n",trash,&((*parameters).drive_mag));
  fscanf(in,"%s %lf\n",trash,&((*parameters).drive_frq));
  fscanf(in,"%s %d\n",trash,&((*parameters).decifactor));
  
  //starting from a file?
  fscanf(in,"%s %d\n",trash,&((*parameters).restart));

  //ramping the drive rate?
  fscanf(in,"%s %d\n",trash,&((*parameters).drive_step_time));
  fscanf(in,"%s %lf\n",trash,&((*parameters).drive_step_force));
  
  fclose(in);
}

//-----------------------------------------------------------------
void distance(double *dr,double *dx,double *dy,double x1,double y1,
	      double x2,double y2,struct syssize syssize)
{
  double locdx,locdy;

  locdx=x1-x2;
  if(locdx>syssize.SX2) locdx-=syssize.SX;
  if(locdx<=-syssize.SX2) locdx+=syssize.SX;
  locdy=y1-y2;
  if(locdy>syssize.SY2) locdy-=syssize.SY;
  if(locdy<=-syssize.SY2) locdy+=syssize.SY;
  *dr=sqrt(locdx*locdx+locdy*locdy);
  *dx=locdx;
  *dy=locdy;
  
  return;
}


