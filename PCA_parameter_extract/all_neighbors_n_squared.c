/* 6.19.19 
 * Revised to measure not only dr, but sort dr before printing.
 * also calculating the average
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

  //you could put dr here and save... but I don't think it would help.
  
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

double pair_interact(int id1,int id2,struct vortex *vortex,
		   struct syssize syssize,
		   struct parameters parameters,
		   FILE *out);

int compare (double *num1, double *num2);

void name_outfile(int time, char *ascii_file, char *filename);

int read_frame(FILE *in,int *nV,int *time,struct vortex *vortex);


//----------------------------------------------------------------//


int main(int argc,char *argv[])
{
  //file names/pointers
  FILE *in,*out,*feature_out; //

  //for naming files and reading in with fscanf
  char ascii_file[120]; 
  char feature_file[120]; 

  char str_time[10];
  int read_items;
  char trash[120];
  int trash_int;
  float dr;
  
  //if directory doesn't exist, make directory
  //make and/or open directory snapshot_files
  struct stat st = {0};       //null directory
  if (stat("neighbor_data", &st) == -1) {
    mkdir("neighbor_data",0700);
  }
  
  //simulation parameters
  double density;
  struct syssize syssize;
  int runtime;
  int maxtime;
  int writemovietime;
  
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
  writemovietime = parameters.writemovietime;

  nV = maxnum;
  
  //more system specific.  vortex struct only contains maxnum
  //and we create room for a neighbor list at the end.
  vortex=malloc(maxnum*sizeof(struct vortex));

  //make an array to calculate all of the averages
  //within every frame.
  //so every frame we zero all the values
  //and then do a running total within the frame
  double rN_average[maxnum-1];
  double dr_probe[maxnum-1];

  int arraylen = sizeof(dr_probe) / sizeof(double);

  if(DEBUG)
    printf("%d %d\n",nV,arraylen);
  
  //for each new frame, the average neighbor distances are reset to zero
  for(i=0;i<(nV-1);i++){
    
    //do a += continuously as you accumulate each frame -
    //normalize at end of smtest
    rN_average[i] = 0.0;

    //this is a particle by particle measure,
    //so it gets filled up and rewritten for every particle
    //so you need to write it somewhere - either an array or file
    //and then reload it to subtract off normalized values later
    dr_probe[i] = 0.0;
  }
    
  //open the smtest file
  if((in=fopen(filename,"r"))==NULL){
    printf("Error opening file %s\n",filename);
    exit(-1);
  }

  //create a list of probe particle centered neighbor distances
  //each particle has its own list of neighbors,
  //sorted from small to large
  //and each value is normalized across particles
  //such that nn=1 is averaged and subtracted, likewise nn=2, etc
  

  if(DEBUG){
    f_debug = fopen("debug.dat","w");
  }
  
  //TODO-----------------------------------------
  //total movie stats - update for new measures
  frame=0;
  
  //read in all frames of smtest and analyze individually
  while(!feof(in)){
    frame++;  
    nVold=nV;
    completeframe=read_frame(in,&nV,&time,vortex);
    
    //cut early frames
    if(!completeframe) break;

    //often this is a square array, we don't want stats on
    if(time==0) continue;

    //name file to write data - we'll do in subroutine to clean up
    name_outfile(time,ascii_file,"neighbor_data/frame_");    

    if((out=fopen(ascii_file,"rw"))==NULL){
      printf("Error opening file %s\n",ascii_file);
      exit(-1);
    }

    //zero all rN data before new particle analysis double loop
    //for each new frame, the average neighbor distances are reset to zero
    for(i=0;i<(nV-1);i++){
    
      //do a += continuously as you accumulate each frame -
      //normalize at end of FRAME D
      rN_average[i] = 0.0;
    }

    //------------------------------------
    //loop through all particles
    for(i=0;i<nV;i++){

      //zero the neighbor counter for every new probe particle i
      k=0;
      
      //fprintf(out,"#particle %d\n",i);

      //and all their neighbors -
      //each particle is a distinct probe, so we don't do j<i
      for(j=0;j<nV;j++){

	if(i!=j){
	  //calculate the dr for each pair
	  dr_probe[k] = pair_interact(i,j,vortex,syssize,parameters,out);
	  k++;
	}
	
      }//exit the j portion of the double loop, the neighbors of i

      //make some i probe specific calculations
      
      //qsort the dr values
      qsort (dr_probe,arraylen,sizeof(double),
	     (int (*)(const void *, const void *)) compare);

      //before we write over the dr_probe values for the new j
      //write the values to the array to calculate their average 
      for(j=0;j<(nV-1);j++){
	fprintf(out,"%lf\n",dr_probe[j]);
	rN_average[j] += dr_probe[j];
      }

    }//exit the i loop over all nV probe particles

    //-------------------------------------------------
    //end the data accumlation - frame level analysis
    //-------------------------------------------------
    
    //normalize the average for the frame D (not over the simulation S).
    for(i=0;i<(nV-1);i++){
      rN_average[i] /= (nV-1);
    }

    //read the data back in from the written frame_data and write feature_data
    //REWIND!
    rewind(out);

    //open a new file for writing the feature data -
    //frame by frame for downsampling
    name_outfile(time,feature_file,"neighbor_data/feature_");    

    if((feature_out=fopen(feature_file,"w"))==NULL){
      printf("Error opening file %s\n",feature_file);
      exit(-1);
    }

      
    //FILE STRUCTURE - COMMENT LINE FOR EACH nV particles
    //followed by nV-1 data lines
    //i.e. nV lines per particle, of which there are nV, AKA (nV^2)

    for(i=0;i<nV;i++){
      
      //read_items = fscanf(out,"%s %d",trash,&trash_int);
      //printf("%s %d\n",trash,trash_int);
      //fflush(stdout);
      //exit(0);
      
      fprintf(feature_out,"#particle %d\n",i);
      
      for(j=0;j<(nV-1);j++){
	read_items = fscanf(out,"%f",&dr);
	fprintf(feature_out,"%f\n",dr - rN_average[j]);
      }//end j neighbor loop
      
    }//end i probe loop
    
    fclose(out);
    fclose(feature_out);
    
  }//do the same for a new frame (i.e. end while(feof) loop)


  //NORMALIZE THE AVERAGES
  //why oh why is it frame-2?
  //time 0 = frame 1, time 1 = frame 2. etc, so frames-1?
  //nope!  get another increment from the while loop going one last time.
  
  //for(i=0;i<(nV-1);i++){
  //  rN_average[i] /= ((frame-2)*(nV-1)) ;
  //}

  //-----------------------------------------------------------
  //reopen the written data and subtract off the averages...
  //---------------------------------------------------------------

  //for(time=writemovietime; time<maxtime; time+=writemovietime){
  //fclose(out);
  //}//end time loop
  
  //--------------------------------------------------------

  //  fclose(out);

  if(DEBUG){
    fclose(f_debug);
  }

  return 0;
}
//--------------------------------------------------------------------


//------------------------------------------------------------------
int read_frame(FILE *in,int *nV,int *time,struct vortex *vortex)
{
  int i;
  float xin,yin,zin;

  int read_items;

  read_items = fread(nV,sizeof(int),1,in);
  read_items = fread(time,sizeof(int),1,in);
  if(feof(in)) return 0;
  for(i=0;i<*nV;i++){
    read_items = fread(&(vortex[i].color),sizeof(int),1,in);
    read_items = fread(&(vortex[i].id),sizeof(int),1,in);
    read_items = fread(&xin,sizeof(float),1,in);
    read_items = fread(&yin,sizeof(float),1,in);
    read_items = fread(&zin,sizeof(float),1,in);
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
double pair_interact(int id1,int id2,
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
  //fprintf(out,"%f\n",dr);    
	

  return dr;
  
}//end subroutine.  return dr by value

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

  //to eliminate all the damn warnings
  int read_int;
  
  resolution=1e-6;

  if((in=fopen("Pa0","r"))==NULL){
    printf("Input file Pa0 not found\n");
    exit(-1);
  }
  else if (VERBOSE) printf("reading Pa0.\n");
  fflush(stdout);

  read_int = fscanf(in,"%s %lf\n",trash,&((*parameters).density));
  //fscanf(in,"%s %lf\n",trash,&((*parameters).phi2_small));
  //fscanf(in,"%s %lf\n",trash,&((*parameters).phi1_big));
  read_int =fscanf(in,"%s %lf\n",trash,&((*parameters).pdensity));

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
  read_int = fscanf(in,"%s %lf\n",trash,&((*syssize).SX));
  (*syssize).SX2=(*syssize).SX*0.5;
  //--
  read_int = fscanf(in,"%s %lf\n",trash,&((*syssize).SY));
  (*syssize).SY2=(*syssize).SY*0.5;
  
  read_int = fscanf(in,"%s %lf\n",trash,&((*parameters).radius)); //_small));
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
 
  read_int = fscanf(in,"%s %d\n",trash,&((*parameters).runtime));
  read_int = fscanf(in,"%s %lf\n",trash,&((*parameters).runforce));
  read_int = fscanf(in,"%s %lf\n",trash,&((*parameters).dt));
  read_int = fscanf(in,"%s %d\n",trash,&((*parameters).maxtime));
  if(VERBOSE) printf("\n maxtime is: %d\n", (*parameters).maxtime);

    
  read_int = fscanf(in,"%s %d\n",trash,&((*parameters).writemovietime));
  read_int = fscanf(in,"%s %lf\n",trash,&((*parameters).kspring));
  
  read_int = fscanf(in,"%s %lf\n",trash,&cellsize); // size of lookup cell
  read_int = fscanf(in,"%s %lf\n",trash,&((*parameters).potential_radius));
  read_int = fscanf(in,"%s %lf\n",trash,&((*parameters).potential_mag));
  
  read_int = fscanf(in,"%s %lf\n",trash,&length_scale);
  read_int = fscanf(in,"%s %lf\n",trash,&((*parameters).drive_mag));
  read_int = fscanf(in,"%s %lf\n",trash,&((*parameters).drive_frq));
  read_int = fscanf(in,"%s %d\n",trash,&((*parameters).decifactor));
  
  //starting from a file?
  read_int = fscanf(in,"%s %d\n",trash,&((*parameters).restart));

  //ramping the drive rate?
  read_int =fscanf(in,"%s %d\n",trash,&((*parameters).drive_step_time));
  read_int = fscanf(in,"%s %lf\n",trash,&((*parameters).drive_step_force));
  
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


//------------------------------------------------------------
int compare (double *num1, double *num2)
{
  if (*num1 < *num2) return -1;
  else if (*num1 == *num2) return 0;
  return 1;
}

//--------------------------------------------------------------
void name_outfile(int time, char *ascii_file, char *filename){

  //char ascii_file[120]; 
  char str_time[10];
  //open file to write data
  strcpy(ascii_file,filename);
  sprintf(str_time,"%08d",time); //convert current to a string   
  strcat(ascii_file,str_time);
  
  if(DEBUG)
    printf("printing frame data to %s\n",ascii_file);

  //------------------------------------------------------------
  //CAUSES SEGFAULT - EASER TO DO IT IN MAIN FUNCTION
  //------------------------------------------------------------
  //if((out=fopen(ascii_file,"w"))==NULL){
  //  printf("Error opening file %s\n",ascii_file);
  //  exit(-1);
  //}

  return;
}
