/*1.23.19
 *making adding an n^2 algorithm to find all nn.
 */

/* Revisions log:
 * 1.06.19 printing out r_nn, x_nn, y_nn and calculating mean, <r_nn>
           \delta r_nn' = r_nn - <r_nn>
	   sorting from small to large(?)

 * 7.22.17 use cell N log N algorithm to find near neighbors, 
           measure the distances, 
           look at near neighbors and angles
           also consider \psi_6

	   also clean up the code to remove the older variables
	   that are no longer used


 * 5.13.16 default input name "smtest"
 * 7.17.15 Replacing ancient Stuart/Jared lookup table with my newly
           written version.
 * 7.15.15 Find overall largest cluster also.
 * 6.19.15 Writing a read_ascii subroutine to test float->double issues
 * 9.6.13  Well, that didn't work.  Next, I try using Hermann's method.
 * 8.29.13 Written */

#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<string.h>

#define PI 3.14159265359

//hardwired - make lots of space!
#define MAXPAR 10000       //max particles!
#define MAXSTACK 50000     //max NOT USED!

#define DEBUG 0
#define VERBOSE 0

FILE *f_debug;

struct vortex{
  int id;
  int color;
  double x;
  double y;
  double radius;
  int clusterid;

  //DM measurements
  int near_neighbor_id;
  double near_neighbor_dr;
  double near_neighbor_dx;
  double near_neighbor_dy;
  double near_neighbor_angle;
  
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
  int nV1, nV2;
  int maxnum;
  
  //note: improve the code using
  //      cynthia's old granular/active code had arrays
  //      to hold the many possible disk radii
  //      cleaner, but harder to read.  
  int maxnum_small;
  int maxnum_large;
  
  double radius_small;
  double radius_large;
  double runforce;
  int runtime;
  //
  double density;
  double phi2_small;
  double phi1_big;
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

void cell_interact(struct vortex *vortex,int nV,
		   struct mylookup **lookuptable,
		   struct lookupdata lookupdata,
		   struct syssize syssize,
		   int **distarray);

//for clusters distarray measures overlap (binary)
//here we want to look within a longer distance
//they don't actually have to be touching
//and then find the closest


void distance(double *dr, double *dx, double *dy,
	      double x1, double y1,
	      double x2, double y2,
	      struct syssize syssize);

//void get_parameters_file(double *density, struct syssize *syssize,
//			 int *runtime, int *maxtime,
//			 struct lookupdata *lookupdata,
//			 int *maxcell,int *maxnum);

//DM read in Pa0 
void get_parameters_file(struct parameters *parameters,
			 struct syssize *syssize,
			 struct lookupdata *lookupdata,
			 int *maxcell);

void map_particles_to_cells(struct vortex *vortex,int nV,
			    struct mylookup **lookuptable,
			    struct lookupdata lookupdata);

void pair_interact(int id1,int id2,struct vortex *vortex,
		   struct syssize syssize,int **distarray);




main(int argc,char *argv[])
{
  FILE *in,*out; //,*out2,*out3;
  double density;
  struct syssize syssize;
  int runtime;
  int maxtime;

  //particle cell data 
  int maxcell;
  struct mylookup **lookuptable;
  struct lookupdata lookupdata;

  struct parameters parameters;
  
  //
  int completeframe;
  struct vortex *vortex;
  int maxnum;
  
  char filename[120]="smtest";
  int nV, time, nVold;
  int i,j,k;
  //int *stack;
  int **distarray;
  int maxcluster1;
  int maxcluster2;
  int overallmax;
  int frame;
  int avgmax;
  int maxcount;

  //DM, calculate the average angle for each frame individually
  double avg_dr_nn_frame;
  double avg_dx_nn_frame;
  double avg_dy_nn_frame;
  float counter; //frame counter - make float since we'll be dividing with it
  
  int **cluster;
  int length[MAXPAR];
  int null;

  null=-1000;

  //generic malloc for some unknown number of particles (MAXPAR)
  //will zero and populate for every frame
  cluster=(int **)malloc((int)((MAXPAR)*sizeof(int *)));
  cluster[0]=(int *)malloc((int)((MAXPAR*MAXPAR)*sizeof(int)));
  for(i=1;i<=MAXPAR;i++){
    cluster[i]=cluster[i-1]+MAXPAR;
  }

  //generic malloc for some unknown number of particles (MAXPAR)
  //will zero and populate for every frame
  distarray=(int **)malloc((int)((MAXPAR)*sizeof(int *)));
  distarray[0]=(int *)malloc((int)((MAXPAR*MAXPAR)*sizeof(int)));
  for(i=1;i<=MAXPAR;i++){
    distarray[i]=distarray[i-1]+MAXPAR;
  }

  //be careful on this one...
  //get_parameters_file(&density,&syssize,&runtime,
  //&maxtime,&lookupdata,&maxcell,&maxnum);

  get_parameters_file(&parameters,&syssize,&lookupdata,&maxcell);

  maxnum = parameters.maxnum;
  maxtime = parameters.maxtime;
  runtime = parameters.runtime;
  density = parameters.density;
  
  //more system specific.  vortex struct only contains maxnum, not MAXPAR
  vortex=malloc(maxnum*sizeof(struct vortex));

  //likewise lookupdata comes from Pa0
  //make (malloc) the lookup cells
  lookuptable=(struct mylookup **)malloc(lookupdata.ncellsx
					 *sizeof(struct mylookup *));

  for(i=0;i<lookupdata.ncellsx;i++){
    lookuptable[i]=malloc(lookupdata.ncellsy*sizeof(struct mylookup));
  }
  
  for(j=0;j<lookupdata.ncellsy;j++){
    for(i=0;i<lookupdata.ncellsx;i++){
      lookuptable[i][j].num=0;
      lookuptable[i][j].list=malloc(maxcell*sizeof(int));
    }
  }

  //now that we've made all these structs in memory,
  //we can add the information from smtest
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

  out=fopen("delta_r_nn.dat","w");

  if(DEBUG){
    f_debug = fopen("debug.dat","w");
  }
  
  //TODO-----------------------------------------
  //total movie stats - update for new measures
  avgmax=0;
  maxcount=0;
  frame=0;
  overallmax=0;
  
  //read in all frames of smtest and analyze individually
  while(!feof(in)){
    frame++;
    nVold=nV;
    completeframe=read_frame(in,&nV,&time,vortex);

    if(!completeframe) break;
    if(frame==1) continue;
    
    // Prepare a distance array.
    for(i=0;i<nV;i++){
      for(j=0;j<nV;j++){
	distarray[j][i]=null;
      }
    }

    map_particles_to_cells(vortex,nV,lookuptable,lookupdata);
    cell_interact(vortex,nV,lookuptable,lookupdata,syssize,distarray);
	  
    // We've now marked all particles that are touching each other.
    // Next we find the near neighbor, by species
    // and characterize the NN distance and angle

    for(i=0;i<nV;i++){

      if(vortex[i].near_neighbor_dr == 1000.0){
	for(j=0;j<i;j++){
	  //printf("%d %d %f %f %f %f\n",i, j, vortex[i].x, vortex[i].y, vortex[j].x, vortex[j].y);
	    pair_interact(i,j,vortex,syssize,distarray);
	}
      }
    }
    
    

    avg_dr_nn_frame = 0.0;
    avg_dx_nn_frame = 0.0;
    avg_dy_nn_frame = 0.0;
    counter=0;
    
    for(i=0;i<nV;i++){

	
      //make sure you haven't counted any of those big initial values
      if(fabs(vortex[i].near_neighbor_dr-1000.0) > 0.1){
	avg_dr_nn_frame += fabs(vortex[i].near_neighbor_dr);
	avg_dx_nn_frame += fabs(vortex[i].near_neighbor_dx);
	avg_dy_nn_frame += fabs(vortex[i].near_neighbor_dy);
      }      
      counter += 1.0;
 
      
    }//end i loop

    avg_dr_nn_frame /= counter;
    avg_dx_nn_frame /= counter;
    avg_dy_nn_frame /= counter;
	  
    fprintf(out,"#time \t dr \t dx \t dy \n");
    fprintf(out,"#%d averages: \t %e \t %e \t %e\n",time, avg_dr_nn_frame, avg_dx_nn_frame, avg_dy_nn_frame);

    for(i=0;i<nV;i++){


      fprintf(out,"%d \t %e \t %e \t %e\n",	\
	      i, 
	      fabs(vortex[i].near_neighbor_dr)-avg_dr_nn_frame,
	      fabs(vortex[i].near_neighbor_dx)-avg_dx_nn_frame,
	      fabs(vortex[i].near_neighbor_dy)-avg_dy_nn_frame);

      // else{
      //fprintf(out,"%d \t %e \t %e \t %e\n", i, 0.0, 0.0, 0.0);
      //}
	
    }
     
    avg_dr_nn_frame = 0.0;
    avg_dx_nn_frame = 0.0;
    avg_dy_nn_frame = 0.0;
    counter = 0.0;
  }//do the same for a new frame (i.e. end while(feof) loop)

  

  //--------------------------------------------------------

  fclose(out);

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
    vortex[i].clusterid=-1;

    vortex[i].near_neighbor_id = i;
    
    //set these to arbitrarily large values intentionally
    vortex[i].near_neighbor_dr = 1000.0;
    vortex[i].near_neighbor_dx = 1000.0;
    vortex[i].near_neighbor_dy = 1000.0;
    vortex[i].near_neighbor_angle = 1000.0;
    
  }
  return 1;
  
}

//------------------------------------------------------------------------
void cell_interact(struct vortex *vortex,int nV,
		   struct mylookup **lookuptable,
		   struct lookupdata lookupdata,
		   struct syssize syssize,
		   int **distarray)
{
  int i,j,k,kk;
  int id1,id2;
  int loop;
  int indx,jndx;

  //loop through the cells (i,j)
  for(i=0;i<lookupdata.ncellsx;i++){
    for(j=0;j<lookupdata.ncellsy;j++){
      
      //within each cell, dig into the number of particles in the cell
      for(k=0;k<lookuptable[i][j].num;k++){

	//id1 is the universal particle number from smtest
	id1=lookuptable[i][j].list[k];

	//do each case below, since nn can be in bordering cells
	for(loop=0;loop<5;loop++){
	  switch(loop){
	  case 0:
	    // Same cell
	    indx=i;
	    jndx=j;
	    break;
	  case 1:
	    // Cell to the right
	    indx=i+1;
	    jndx=j;
	    if(indx>=lookupdata.ncellsx) indx-=lookupdata.ncellsx;
	    break;
	  case 2:
	    // Cell to diagonal upper right
	    indx=i+1;
	    jndx=j+1;
	    if(indx>=lookupdata.ncellsx) indx-=lookupdata.ncellsx;
	    if(jndx>=lookupdata.ncellsy) jndx-=lookupdata.ncellsy;
	    break;
	  case 3:
	    // Cell above
	    indx=i;
	    jndx=j+1;
	    if(jndx>=lookupdata.ncellsy) jndx-=lookupdata.ncellsy;
	    break;
	  case 4:
	    // Cell to diagonal upper left
	    indx=i-1;
	    jndx=j+1;
	    if(indx<0) indx+=lookupdata.ncellsx;
	    if(jndx>=lookupdata.ncellsy) jndx-=lookupdata.ncellsy;
	    break;
	  }

	  //loop through all the particles in the neighboring cell
	  for(kk=0;kk<lookuptable[indx][jndx].num;kk++){
	    // Do each interacting pair only once

	    //same cell, same particle
	    if((!loop)&&(kk<=k)) continue;
	    
	    //unique particle
	    id2=lookuptable[indx][jndx].list[kk];

	    //this is the place to test whether both
	    //types are the same radius or not
	    //in other words, what are the nature of the lanes?
	    //if(vortex[id1].color == vortex[id2].color){
	    //}
	    
	    //check the distances of the two unique particles
	    pair_interact(id1,id2,vortex,syssize,distarray);
	  }
	}
      }
    }
  }
}

//---------------------------------------------------------
//---------------------------------------------------------
void pair_interact(int id1,int id2,struct vortex *vortex,
		   struct syssize syssize,int **distarray)
{
  double dr,dx,dy;
  double deltar;
  double epsilon;


  //center to center distance with PBC, returns to dr
  distance(&dr,&dx,&dy,
	   vortex[id1].x,vortex[id1].y,
	   vortex[id2].x,vortex[id2].y,
	   syssize);
  
  //distance between particles
  deltar=dr-(vortex[id1].radius+vortex[id2].radius);

  distarray[id1][id2]=1;
  distarray[id2][id1]=1;

  //test if they are candidates to be nearest neighbors
  //if so, perform the next measure
  if(dr <= vortex[id1].near_neighbor_dr){

    vortex[id1].near_neighbor_id = id2;
    vortex[id1].near_neighbor_dr = dr;

    //take the absolutel value because we care about
    //first quadrant angle.
    vortex[id1].near_neighbor_dx = fabs(dx);
    vortex[id1].near_neighbor_dy = fabs(dy);

    //calculate angle with horizontal
    //choose the first qua
    vortex[id1].near_neighbor_angle = atan2(fabs(dy),fabs(dx));
      
    if(DEBUG){
      fprintf(f_debug,"%d %lf %lf %lf \n",id1,
	      dx,dy,vortex[id1].near_neighbor_angle);
    }

  } //end if
  else if(DEBUG){
    fprintf(f_debug,"%d %lf %lf %lf \n",id1,
	    dx,dy,dr);
  }
	    
    
  if(dr < vortex[id2].near_neighbor_dr){

    vortex[id2].near_neighbor_id = id1;
    vortex[id2].near_neighbor_dr = dr;
    vortex[id2].near_neighbor_dx = fabs(dx); //same as particle i
    vortex[id2].near_neighbor_dy = fabs(dy); //same as particle i

    //calculate angle with horizontal
    //probably just want the absolute value
    vortex[id2].near_neighbor_angle = atan2(fabs(dy),fabs(dx));

    if(DEBUG && fabs(dx)>0.0){
      fprintf(f_debug,"AACK! %d %lf %lf %lf %lf \n",id2,
	      -dx,-dy,dr,vortex[id2].near_neighbor_dr );
    }
	    
  }//end if
	

  return;
  
}//end subroutine.  

//-------------------------------------------------------------
//-------------------------------------------------------------
void map_particles_to_cells(struct vortex *vortex,int nV,
			    struct mylookup **lookuptable,
			    struct lookupdata lookupdata)
{
  int i,j,k;
  int xindex,yindex;

  //initialize the cells with no particles within
  for(j=0;j<lookupdata.ncellsy;j++){
    for(i=0;i<lookupdata.ncellsx;i++){
      lookuptable[i][j].num=0;
    }
  }

  //run through all of the particles,
  //adding them to the proper cell
  //in a discrete grid
  for(i=0;i<nV;i++){

    //here we map to the cell grid
    xindex=floor(vortex[i].x*lookupdata.xscale);
    yindex=floor(vortex[i].y*lookupdata.yscale);

    //here we put the nth element into the 0th
    if(xindex==lookupdata.ncellsx)
      xindex=0;
    if(yindex==lookupdata.ncellsy)
      yindex=0;

    //here we get the number of particles in the cell
    k=lookuptable[xindex][yindex].num;

    //the lookuptable has a list that the kth particle
    //has the following id from smtest
    lookuptable[xindex][yindex].list[k]=i;

    //now we increment the number of particles in the cell
    (lookuptable[xindex][yindex].num)++;
    
  }//end for(i) loop

  return;
}

//--------------------------------------------------------------------
//-----------------------------------------------------------------
void get_parameters_file(struct parameters *parameters,
			 struct syssize *syssize,
			 struct lookupdata *lookupdata,
			 int *maxcell)
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
  fscanf(in,"%s %lf\n",trash,&((*parameters).phi2_small));
  fscanf(in,"%s %lf\n",trash,&((*parameters).phi1_big));
  fscanf(in,"%s %lf\n",trash,&((*parameters).pdensity));

  //----------------------------------------------------------------------
  //check whether the user gave you sensible choices for the two densities
  //----------------------------------------------------------------------
  double test1 = (*parameters).phi2_small + (*parameters).phi1_big;
  if( fabs( ((*parameters).density-test1) >0.01) ){
    printf("Inconsistency in requested densities.\n");
    printf("Total density: %f\n", (*parameters).density);
    printf("Small disk density: %f\n", (*parameters).phi2_small);
    printf("Large disk density: %f\n", (*parameters).phi1_big);
    exit(-1);
  }
  //---------------------------------------------------------------------
  
  if (VERBOSE) printf("reading parameters density %f, pin density %f.\n", (*parameters).density, (*parameters).pdensity);
  fflush(stdout);

  //-
  fscanf(in,"%s %lf\n",trash,&((*syssize).SX));
  (*syssize).SX2=(*syssize).SX*0.5;
  //--
  fscanf(in,"%s %lf\n",trash,&((*syssize).SY));
  (*syssize).SY2=(*syssize).SY*0.5;
  
  fscanf(in,"%s %lf\n",trash,&((*parameters).radius_small));
  fscanf(in,"%s %lf\n",trash,&((*parameters).radius_large));

  //--where it is appropriate to accurately calculate maxnum
  //--based on particle sizes

  double prefix = (*syssize).SX*(*syssize).SY/PI;
  double phiS_term = (*parameters).phi2_small/( (*parameters).radius_small*(*parameters).radius_small) ;
  double phiB_term = (*parameters).phi1_big/( (*parameters).radius_large*(*parameters).radius_large);
   
  (*parameters).maxnum= (int) floor(prefix * (phiS_term + phiB_term));
  
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

  //Create the cells for ln(N) lookups
  (*lookupdata).ncellsx=(int)(*syssize).SX/cellsize;
  (*lookupdata).ncellsy=(int)(*syssize).SY/cellsize;
  
  if(fabs((*lookupdata).ncellsx*cellsize-(*syssize).SX)>resolution){
    printf("Error, SX mismatch with cell size\n");
    //printf("Error, int(SX) = %f\n",(int)(*syssize).SX/cellsize);
    fflush(stdout);

    //printf("Error, calc: %f\n",fabs((*lookupdata).ncellsx*cellsize-(*syssize).SX));
    exit(-1);
  }
  
  if(fabs((*lookupdata).ncellsy*cellsize-(*syssize).SY)>resolution){
    printf("Error, SY mismatch with cell size\n");
    exit(-1);
  }
  
  (*lookupdata).xscale=(double)(*lookupdata).ncellsx/(*syssize).SX;
  (*lookupdata).yscale=(double)(*lookupdata).ncellsy/(*syssize).SY;
  
  // Compute largest number of particles that could possibly fit
  // inside a cell of this size
  *maxcell=(int)2*cellsize*cellsize/(PI*((*parameters).radius_small)*((*parameters).radius_small));
  
  //printf("maxcell:%d\n",*maxcell); //YY debuging the size of lookuptable
  
  fclose(in);
}

/*
//--------------------------------------------------------------------
void get_parameters_file(double *density,struct syssize *syssize,
			 int *runtime,int *maxtime,struct lookupdata *lookupdata,
			 int *maxcell,int *maxnum)
{
  FILE *in;
  char trash[120];
  double pdensity;
  double cellsize;
  double resolution;
  double radius;
  double runforce;
  double dt;
  //int maxtime;  //DM I need this
  int writemovietime;
  double kspring;

  resolution=1e-6;

  if((in=fopen("Pa0","r"))==NULL){
    printf("Input file Pa0 not found\n");
    exit(-1);
  }
  fscanf(in,"%s %lf\n",trash,density);
  fscanf(in,"%s %lf\n",trash,&pdensity);
  fscanf(in,"%s %lf\n",trash,&((*syssize).SX));
  (*syssize).SX2=(*syssize).SX*0.5;
  fscanf(in,"%s %lf\n",trash,&((*syssize).SY));
  (*syssize).SY2=(*syssize).SY*0.5;
  *maxnum=*density*(*syssize).SX*(*syssize).SY;
  fscanf(in,"%s %lf\n",trash,&radius);
  fscanf(in,"%s %d\n",trash,runtime);
  fscanf(in,"%s %lf\n",trash,&runforce);
  fscanf(in,"%s %lf\n",trash,&dt);
  fscanf(in,"%s %d\n",trash,maxtime);
  fscanf(in,"%s %d\n",trash,&writemovietime);
  fscanf(in,"%s %lf\n",trash,&kspring);
  fscanf(in,"%s %lf\n",trash,&cellsize); // size of lookup cell
  (*lookupdata).ncellsx=(int)(*syssize).SX/cellsize;
  (*lookupdata).ncellsy=(int)(*syssize).SY/cellsize;
  if(fabs((*lookupdata).ncellsx*cellsize-(*syssize).SX)>resolution){
    printf("Error, SX mismatch with cell size\n");
    exit(-1);
  }
  if(fabs((*lookupdata).ncellsy*cellsize-(*syssize).SY)>resolution){
    printf("Error, SY mismatch with cell size\n");
    exit(-1);
  }
  (*lookupdata).xscale=(double)(*lookupdata).ncellsx/(*syssize).SX;
  (*lookupdata).yscale=(double)(*lookupdata).ncellsy/(*syssize).SY;
  // Compute largest number of particles that could possibly fit
  // inside a cell of this size
  *maxcell=2*(int)cellsize*cellsize/(PI*radius*radius);
  fclose(in);

  return;
}
*/

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

