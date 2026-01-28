#include <riegl/scanlib.hpp>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <exception>
#include <cmath>
#include <limits>
#include <memory>
#include <string>
#include <cstring>
#include <cstdlib>
/*#include <libLidarHDF.h>*/

#ifdef _MSC_VER 
//not #if defined(_WIN32) || defined(_WIN64) because we have strncasecmp in mingw
#define strncasecmp _strnicmp
#define strcasecmp _stricmp
#endif

using namespace scanlib;
using namespace std;

/*thresholds for recording points/beams*/
float maxZen,maxDev;
float minRefl, maxRefl;


/*#########################################################*/
/* A class that reads the .rxp files*/

class reader : public pointcloud  /*we have inherited pointcloud public members which in turn inherits "compressed_packets"*/
{
private:
  ostream& o;
  ostream& oAsc;
  unsigned long line;
  int doneThis;
  int donethat=1;
  float zen,az;
  float range[20],refl[20], refl_mb[20];
  float x_mb[20], y_mb[20], z_mb[20];
  float tic;
  float dev;         /*deviation of pulse shape*/
  float scanCent[3];
  double trans[4][4];      /*4*4 translation matrix for geolocation*/
  void transform();
  uint8_t nHits;
public:
  uint32_t shotN,totShot;         // total number of shots
  double coordOffset[3];
  double bounds[6];
  int ascOut,binOut;
  int nDec;
  double M_PI = 3.14159265358979323846;
  bool hdf;   // HDF5 writing switch


  // the constructor function
  reader(ostream& o_,ostream& oAsc_)
        : pointcloud(false) // set this to true if you need gps aligned timing. What happens in the poincoud constructor?
        , o(o_)             // initialises o as o_, which is what is passed to the class constructor.
        , oAsc(oAsc_)
  {
    o.precision(10);    /*sets the number of decimal places of output*/
    line=0;
    shotN=0;
    totShot=0;
    
  }/*creator function*/

  /*This call is invoked for every pulse, even if there is no retur, by importer, before on_echo_transformed()*/
  void on_shot()   // this seems to be an overloading of a pointcloud member function.
  {
    /*transform vector and origin to geolocate*/
    transform();
    int i = 0;

    //write bounds to file header
    if (donethat) {
    if (ascOut) {
        oAsc << std::setprecision(10) << bounds[0] << " " << bounds[1] << " " << bounds[2] << " " << bounds[3] << " " << bounds[4] << " " << bounds[5] << endl;
    
    }
    if (binOut) {
        for (i = 0; i < 6; i++)o.write((char*)&bounds[i], sizeof(double));
    }
    donethat = 0;
    }

    if(shotN==0){
      coordOffset[0]=beam_origin[0];
      coordOffset[1]=beam_origin[1];
      coordOffset[2]=beam_origin[2];
    }

    zen=(float)(atan2(sqrt(beam_direction[0]*beam_direction[0]+beam_direction[1]*beam_direction[1]),beam_direction[2])*180.0/M_PI);
    az=(float)(atan2(beam_direction[0],beam_direction[1])*180.0/M_PI);
    if((zen>361.0)||(zen<-361.0)||(az>361.0)||(az<-361.0))fprintf(stderr,"Angle error %f %f\n",zen,az);
    scanCent[0]=(float)(beam_origin[0]-coordOffset[0]);
    scanCent[1]=(float)(beam_origin[1]-coordOffset[1]);
    scanCent[2]=(float)(beam_origin[2]-coordOffset[2]);

    if (ascOut)writePulse();
    if (binOut)writePulse_binary();

    nHits=0;
    return;
  }/*on_shot*/

  /*Invoked for every return by importer, after on_shot()*/
  void on_echo_transformed(echo_type echo)   // this seems to be an overloading of a pointcloud member function.
  {
    target& t(targets[target_count-1]);

    /*check deviation before recording*/
    if(t.deviation<maxDev){
      if(nHits>=20){
        fprintf(stderr,"Too many hits in this beam\n");
        exit(1);
      }

      range[nHits]=(float)t.echo_range;
      refl[nHits]=pow(10.0,t.reflectance/10.0);
      refl_mb[nHits] = t.reflectance;
      x_mb[nHits] = t.vertex[0];
      y_mb[nHits] = t.vertex[1];
      z_mb[nHits] = t.vertex[2];

    


      /*write point cloud*/
      if(ascOut&&((shotN%nDec)==0)){
        dev=t.deviation;
        tic = t.time;
        writePoint();
      }

      if (binOut && ((shotN % nDec) == 0)) {
          dev = t.deviation;
          tic = t.time;
          writePoint_binary();
      }
      

      nHits++;
    }/*deviation checker*/
    return;
  }/*on_echo_transformed*/


  /*Check if we're done with a pulse*/
  void on_shot_end()   /*called at end of pulse*/
  {
    /*this is where we will call the voxelisation program*/
    if((zen<maxZen)&&((totShot%nDec)==0)){
      //MB took out if(binOut)writeBeam();   /*write data to binary file*/
      shotN++;
    }
    totShot++;
    return;
  }/*on_shot_end*/

  /*reset markers*/
  void reset(){doneThis=1;}

  /*read translation matrix*/
  void readTrans(char *,char *);

  /*write a beam to binary output*/
  //void writeBeam();

  /*write point to ASCII file*/
  void writePoint();

  void writePulse();

  void writePulse_binary();

  void writePoint_binary();

  /*end of class*/
};/*reader class*/


/*#######################################################################*/
/*write a binary beam*/

/*void reader::writeBeam()
{
  int i=0,len=25+8*nHits;
  int offset=0;
  //char data[len];
  char* data = new char[len];

  //out data into char array and write it
  memcpy(&(data[0]),&zen,4);
  memcpy(&(data[4]),&az,4);
  memcpy(&(data[8]),&(scanCent[0]),4);
  memcpy(&(data[12]),&(scanCent[1]),4);
  memcpy(&(data[16]),&(scanCent[2]),4);
  memcpy(&(data[20]),&shotN,4);
  memcpy(&(data[24]),&nHits,1);

  /*if((nHits==0)&&((zen>90.0)||(zen<=-90.0))){
    cout << zen<<" "<<az<< " " << scanCent[0] << " " << scanCent[1] << " " << scanCent[2] << endl;
  }*/
/*
  offset=25;
  for(i=0;i<nHits;i++){
    memcpy(&(data[offset]),&(range[i]),4);
    offset+=4;
    memcpy(&(data[offset]),&(refl[i]),4);
    offset+=4;
  }
  o.write(data,len);

  return;
}/*writeBeam*/

void reader::writePulse_binary()
{
    int i = 0;// , len = 25;// +8 * nHits;
    //int offset = 0;
    //char data[len];
    //char* data = new char[len];

    /*out data into char array and write it*/
    /*memcpy(&(data[0]), &i, 1);
    memcpy(&(data[1]), &(beam_origin[0]), 4);
    memcpy(&(data[5]), &(beam_origin[1]), 4);
    memcpy(&(data[9]), &(beam_origin[2]), 4);
    memcpy(&(data[13]), &(beam_direction[0]), 4);
    memcpy(&(data[17]), &(beam_direction[1]), 4);
    memcpy(&(data[21]), &(beam_direction[2]), 4);*/
    


   /* oAsc << std::setprecision(10) << "0 " << beam_origin[0] << " " << beam_origin[1] << " " << beam_origin[2] << " "
        << beam_direction[0] << " " << beam_direction[1] << " " << beam_direction[2] << endl;/*

    /*if((nHits==0)&&((zen>90.0)||(zen<=-90.0))){
      cout << zen<<" "<<az<< " " << scanCent[0] << " " << scanCent[1] << " " << scanCent[2] << endl;
    }*/

    /*offset = 25;
    for (i = 0; i < nHits; i++) {
        memcpy(&(data[offset]), &(range[i]), 4);
        offset += 4;
        memcpy(&(data[offset]), &(refl[i]), 4);
        offset += 4;
    }*/
    //fprintf(stdout, " data : %c beam_origin[0]: %f \n", data, beam_origin[0]);
    //o.write(data, len);
    o.write((char*)&i, sizeof(int));
    o.write((char*)&beam_origin[0], sizeof(double));
    o.write((char*)&beam_origin[1], sizeof(double));
    o.write((char*)&beam_origin[2], sizeof(double));
    o.write((char*)&beam_direction[0], sizeof(double));
    o.write((char*)&beam_direction[1], sizeof(double));
    o.write((char*)&beam_direction[2], sizeof(double));


    return;
}/*writePulse_binary*/

void reader::writePoint_binary()
{
    double x = 0, y = 0, z = 0;
    double zenRad = 0, azRad = 0;
    double reflect;
    int i = 2;// , len = 13;// +8 * nHits;
   //int offset = 0;
   //char data[len];
    //char* data = new char[len];

    zenRad = (double)zen * M_PI / 180.0;
    azRad = (double)az * M_PI / 180.0;

    /*is this correct?*/
    x = sin(zenRad) * sin(azRad) * (double)range[nHits] + (double)scanCent[0] + coordOffset[0];
    y = sin(zenRad) * cos(azRad) * (double)range[nHits] + (double)scanCent[1] + coordOffset[1];
    z = cos(zenRad) * (double)range[nHits] + (double)scanCent[2] + coordOffset[2];

    
    

    if ((x >= bounds[0]) && (y >= bounds[1]) && (z >= bounds[2]) && (x <= bounds[3]) && (y <= bounds[4]) && (z <= bounds[5]) && (nHits == 0)) {
        //oAsc << std::setprecision(10) << "1 " << x << " " << y << " " << z << endl;
        reflect = pow(10, (refl_mb[nHits] + 20) / 10); //for some reason this equation differs from the one used by Steven
        //MB: classifying point on the basis of return deviation and reflectance
        if (range[nHits] < 36){
            reflect = reflect * (0.7155 + 0.04562 * range[nHits] - 0.003741 * pow(range[nHits],2) + 0.00007449 * pow(range[nHits], 3));
        }
        if ((reflect> minRefl)&&(reflect < maxRefl) && (dev <25)){
            //the point is a leaf
            i = 1;
        }
         

    /*out data into char array and write it*/
    /*memcpy(&(data[0]), &i, 1);
    memcpy(&(data[1]), &(x), 4);
    memcpy(&(data[5]), &(y), 4);
    memcpy(&(data[9]), &(z), 4);
    

    /*if((nHits==0)&&((zen>90.0)||(zen<=-90.0))){
      cout << zen<<" "<<az<< " " << scanCent[0] << " " << scanCent[1] << " " << scanCent[2] << endl;
    }*/
        o.write((char*)&i, sizeof(int));
        o.write((char*)&x, sizeof(double));
        o.write((char*)&y, sizeof(double));
        o.write((char*)&z, sizeof(double));
    //o.write(data, len);
    }
    
    return;
}/*writePoint_binary*/

/*#######################################################################*/
/*write an ASCII pulse*/

void reader::writePulse()
{
    
    double zenRad = 0, azRad = 0;

    zenRad = (double)zen * M_PI / 180.0;
    azRad = (double)az * M_PI / 180.0;

    /*is this correct?*/
    

    
        oAsc << std::setprecision(10) << "0 " << beam_origin[0] << " " << beam_origin[1] << " " << beam_origin[2] << " " << beam_direction[0] << " " << beam_direction[1] << " " << beam_direction[2] << endl;
        //testing line: oAsc << std::setprecision(10) << "0 " << beam_origin[0] << " " << beam_origin[1] << " " << beam_origin[2] << " " << scanCent[0] << " " << scanCent[1] << " " << scanCent[2] << " " << coordOffset[0] << " " << coordOffset[1] << " " << coordOffset[2] << " " << beam_direction[0] << " " << beam_direction[1] << " " << beam_direction[2] << endl;

    return;
}/*writePulse*/


/*#######################################################################*/
/*write an ASCII point*/

void reader::writePoint()
{
  double x=0,y=0,z=0;
  double zenRad=0,azRad=0;
  double reflect;
  int i = 2;

  zenRad=(double)zen*M_PI/180.0;
  azRad=(double)az*M_PI/180.0;

  /*is this correct?*/
  x=sin(zenRad)*sin(azRad)*(double)range[nHits]+(double)scanCent[0]+coordOffset[0];
  y=sin(zenRad)*cos(azRad)*(double)range[nHits]+(double)scanCent[1]+coordOffset[1];
  z=cos(zenRad)*(double)range[nHits]+(double)scanCent[2]+coordOffset[2];

  if((x>=bounds[0])&&(y>=bounds[1])&&(z>=bounds[2])&&(x<=bounds[3])&&(y<=bounds[4])&&(z<=bounds[5])&& (nHits==0)){
      
          reflect = pow(10, (refl_mb[nHits] + 20) / 10); //for some reason this equation differs from the one used by Steven
        //MB: classifying point on the basis of return deviation and reflectance
          if (range[nHits] < 36) {
              reflect = reflect * (0.7155 + 0.04562 * range[nHits] - 0.003741 * pow(range[nHits], 2) + 0.00007449 * pow(range[nHits], 3));
          }
          if ((reflect > minRefl) && (reflect < maxRefl) && (dev < 25)) {
              //the point is a leaf
              i = 1;
          }
    //oAsc << std::setprecision(10) << i << " "  << x << " " << y << " " << z << " " <<reflect << " " << dev << endl;
          //testing line : oAsc << std::setprecision(10) << i << " " << x << " " << y << " " << z << " " << scanCent[0] << " " << scanCent[1] << " " << scanCent[2] << " " << coordOffset[0] << " " << coordOffset[1] << " " << coordOffset[2] << endl;
         oAsc << std::setprecision(10) << i << " " << x << " " << y << " " << z <<  endl;
  }
  
  //oAsc << std::setprecision(10) << x_mb[nHits] << " " << y_mb[nHits] << " " << z_mb[nHits] << " " << x << " " << y << " " << z << " " << refl[nHits] << " " << zen << " " << az << " " << range[nHits] << " hits " << (float)nHits << " " << dev << endl;
  //mb: above is testing if vertex are the same as the x, y ,z computed by Steven, they are, up to about 7th decimal place

  return;
}/*writePoint*/


/*#######################################################################*/
/*read transformation matrix*/

void reader::readTrans(char *namen,char *globNamen)
{
  int i=0,j=0,k=0;
  ifstream ipoo(namen); /*file is opened by constructor*/
  double tempTrans[4][4],newTrans[4][4];

  /*do we need to apply a local transformation?*/
  if(strncasecmp(namen,"none",4)){
    /*check that file is opened*/
    if(!ipoo.is_open()){
      cerr << "Input file not open " << namen << endl;
      exit(1);
    }

    /*read file*/
    i=0;
    while((!ipoo.eof())&&(i<4)){
      ipoo >> trans[i][0] >> trans[i][1] >> trans[i][2] >> trans[i][3];
      i++;
    }
    ipoo.close();
  }else{   /*unity matrix*/
    for(i=0;i<4;i++){
      for(j=0;j<4;j++){
        trans[i][j]=0.0;
      }
    }
    for(i=0;i<4;i++)trans[i][i]=1.0;
  }

  /*do we need to apply a global transformation too?*/
  if(strncasecmp(globNamen,"none",4)){
    ifstream gpoo(globNamen); /*file is opened by constructor*/
    /*read global ,matrix*/
    if(!gpoo.is_open()){
      cerr << "Input file not open " << globNamen << endl;
      exit(1);
    }
    i=0;
    while((!gpoo.eof())&&(i<4)){
      gpoo >> tempTrans[i][0] >> tempTrans[i][1] >> tempTrans[i][2] >> tempTrans[i][3];
      i++;
    }
    gpoo.close();

    /*apply global matrix*/
    for(i=0;i<4;i++){
      for(j=0;j<4;j++){
        newTrans[i][j]=0.0;
        for(k=0;k<4;k++){
          newTrans[i][j]+=tempTrans[i][k]*trans[k][j];
        }
      }
    }

    /*copy over*/
    for(i=0;i<4;i++){
      for(j=0;j<4;j++){
        trans[i][j]=newTrans[i][j];
      }
    }
  }/*global matrix if needed*/

  return;
}/*readTrans*/


/*#######################################################################*/
/*apply transformation*/

void reader::transform()
{
  double *tempDir=NULL,*tempCent=NULL;
  double *setBits(double *,double **);
  double *rotateVect(double[3],double[4][4]);

  /*rotate*/
  tempDir=rotateVect(beam_direction,trans);
  tempCent=rotateVect(beam_origin,trans);

  /*translate and put into original vectors*/
  beam_origin[0]=tempCent[0]+trans[0][3];
  beam_origin[1]=tempCent[1]+trans[1][3];
  beam_origin[2]=tempCent[2]+trans[2][3];
  beam_direction[0]=tempDir[0];
  beam_direction[1]=tempDir[1];
  beam_direction[2]=tempDir[2];

  /*free arrays*/
  if(tempDir)delete[] tempDir;
  if(tempCent)delete[] tempCent;

  return;
}/*transform*/


/*#######################################################################*/
/*rotate a vector*/

double *rotateVect(double vect[3],double matrix[4][4])
{
  int i=0,j=0;
  double *rotated=NULL;

  rotated=new double[3];

  for(i=0;i<3;i++){
    rotated[i]=0.0;
    for(j=0;j<3;j++)rotated[i]+=vect[j]*matrix[i][j];
  }

  return(rotated);
}/*rotateVect*/


/*#######################################################################*/
/*control structure*/

typedef struct {
  char inName[200];    // I'll need a list of these eventually.
  char transName[200]; // transformation file
  char globName[200];  // global transformation file
  char outNamen[200];  // output filename
  char ascNamen[200];  // ascii output filename
  char ascOut;         // ascii output switch
  char binOut;         // binary output switch
  int nDec;            // amount to decimate by
  bool hdf;            // HDF5 writing switch
  double bounds[6];    // bounds for ASCII output
}control;


/*#######################################################################*/
/*main*/

int main(int argc, char** argv)
{
  int i=0;
  shared_ptr<basic_rconnection> rc;   /*the file object? Exciting bits are in here*/
  control *dimage=NULL;
  control *readCommands(int,char **);
  ofstream oAsc;
  ofstream opoo;

  /*read the command line*/
  dimage=readCommands(argc,argv);

  if(dimage->binOut)opoo.open(dimage->outNamen,ios::binary);
  if(dimage->ascOut)oAsc.open(dimage->ascNamen);

  /*open the file*/
  cout << "# Reading " << dimage->inName << "\n";
  rc=basic_rconnection::create(dimage->inName);
  rc->open();

  /*decode the file*/
  decoder_rxpmarker dec(rc);     /*needs to be declared after rc is opened for constructor*/
  reader            imp(opoo,oAsc);   /*this is a declaration*/
  buffer            buf;

  /*copy switches to class for later*/
  imp.ascOut=dimage->ascOut;
  imp.binOut=dimage->binOut;
  imp.nDec=dimage->nDec;
  for(i=0;i<6;i++)imp.bounds[i]=dimage->bounds[i];
 

  /*read the transformation file*/
  imp.readTrans(dimage->transName,dimage->globName);

  /*loop over data packets. Functions are called within the inherited classes*/
  for(dec.get(buf);!dec.eoi();dec.get(buf)){ /*loops through the data packets*/
    imp.reset();   /*reset markers*/
    /*many shots are called in the packets, so calls need to be in there*/
    imp.dispatch(buf.begin(), buf.end());   /*dispatch is a member of basic_packets*/
  }

  /*close file*/
  if(dimage->binOut){
    /*write out total number of beams written*/
    /*opoo.write((char *)(&imp.coordOffset[0]),8);
    opoo.write((char *)(&imp.coordOffset[1]),8);
    opoo.write((char *)(&imp.coordOffset[2]),8);
    opoo.write((char *)(&imp.shotN),4);*/
    opoo.close();
    cout << "Binary to " << dimage->outNamen << endl;
  }
  if(dimage->ascOut){
    oAsc.close();
    cout << "ASCII to " << dimage->ascNamen << endl;
  }

  /*close file and tidy arrays*/
  rc->close();
  if(dimage){
    free(dimage);
    dimage=NULL;
  }
  return(0);
}/*main*/


/*#######################################################################*/
/*read the command line*/

control *readCommands(int argc,char **argv)
{
  int i=0,j=0;
  control *dimage=NULL;
  void checkArguments(int,int,int,string);

  if(!(dimage=(control *)calloc(1,sizeof(control)))){
    fprintf(stderr,"error in input filename structure.\n");
    exit(1);
  }

  /*defaults*/
  strcpy(dimage->inName,"/mnt/urban-bess/shancock_work/data/bess/ground_truth/riegl/luton/2014-08-05.001.riproject/ScanPos002/140805_094556.rxp");
  strcpy(dimage->transName,"none");
  strcpy(dimage->outNamen,"teast.bin");
  strcpy(dimage->ascNamen,"teast.pts");
  strcpy(dimage->globName,"none");
  dimage->ascOut=0;
  dimage->binOut=1;
  dimage->nDec=1;
  dimage->hdf=false;
  dimage->bounds[0]=dimage->bounds[1]=dimage->bounds[2]=-100000000.0;
  dimage->bounds[3]=dimage->bounds[4]=dimage->bounds[5]=100000000.0;
  maxZen=10000000.0;
  maxDev=10000000.0;
  minRefl = -10000000.0;
  maxRefl= 10000000.0;

  /*read the command line*/
  for (i=1;i<argc;i++){
    if (*argv[i]=='-'){
      if(!strncasecmp(argv[i],"-input",6)){
        checkArguments(1,i,argc,"-input");
        strcpy(dimage->inName,argv[++i]);
      }else if(!strncasecmp(argv[i],"-trans",6)){
        checkArguments(1,i,argc,"-trans");
        strcpy(dimage->transName,argv[++i]);
      }else if(!strncasecmp(argv[i],"-globMat",8)){
        checkArguments(1,i,argc,"-globMat");
        strcpy(dimage->globName,argv[++i]);
      }else if(!strncasecmp(argv[i],"-output",7)){
        checkArguments(1,i,argc,"-output");
        strcpy(dimage->outNamen,argv[++i]);
      }else if(!strncasecmp(argv[i],"-ascii",6)){
        checkArguments(1,i,argc,"-ascii");
        dimage->ascOut=1;
        dimage->binOut=0;
        strcpy(dimage->ascNamen,argv[++i]);
      }else if(!strncasecmp(argv[i],"-dec",4)){
        checkArguments(1,i,argc,"-dec");
        dimage->nDec=atoi(argv[++i]);
      }else if(!strncasecmp(argv[i],"-bounds",7)){
        checkArguments(6,i,argc,"-bounds");
        for(j=0;j<6;j++)dimage->bounds[j]=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-maxZen",7)){
        checkArguments(1,i,argc,"-maxZen");
        maxZen=atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-maxDev",7)){
        checkArguments(1,i,argc,"-maxDev");
        maxDev=atof(argv[++i]);
      }
      else if (!strncasecmp(argv[i], "-minRefl", 7)) {
          checkArguments(1, i, argc, "-minRefl");
          minRefl = atof(argv[++i]);
      }
      else if (!strncasecmp(argv[i], "-maxRefl", 7)) {
          checkArguments(1, i, argc, "-maxRefl");
          maxRefl = atof(argv[++i]);
      }else if(!strncasecmp(argv[i],"-noBin",6)){
        dimage->binOut=0;
      }else if(!strncasecmp(argv[i],"-hdf",4)){
        dimage->hdf=true;
      }else if(!strncasecmp(argv[i],"-help",5)){
        fprintf(stdout,"\n-input name;       input rxp filename\n-trans name;       transformation filename\n-globMat name;    global coordiunate matrix if needed. Unity otherwise\n-output name;      output filename\n-ascii name;       ascii output and filename\n-noBin;            don't output binary\n-maxZen zen;       max zen to trust. Degrees\n-maxDev dev;      maximum return shape deviation to accept\n-dec n;            rate to decimate pointcloud by\n-bounds minX minY minZ maxX maxY maxZ;    bounds for ASCII output\n-hdf;      write output in compressed HDF5 format\n\n");
        exit(1);
      }else{
        fprintf(stderr,"%s: unknown argument on command line: %s\nTry readRXP -help\n",argv[0],argv[i]);
        exit(1);
      }
    }
  }

  return(dimage);
}/*readCommands*/


/*#########################################################################*/
/*check the number of command line arguments*/

void checkArguments(int numberOfArguments,int thisarg,int argc,string option)
{
  int i=0;

  for(i=0;i<numberOfArguments;i++){
    if(thisarg+1+i>=argc){
      cerr << "error in number of arguments for " << option << " option: " << numberOfArguments << "required\n";
      exit(1);
    }
  }
  return;
}/*checkArguments*/

/*the end*/
/*#######################################################################*/

