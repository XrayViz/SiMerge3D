// postproc3d.cpp
//
// Postprocessing of 3D images made by SiMerge3D
//
// Copyright  2021 Deutsches Elektronen-Synchrotron DESY,
//                  a research centre of the Helmholtz Association.
//
// Author: Oleksandr Yefanov (oleksandr.yefanov@desy.de)

#include "stdio.h"
#include "string.h"
#include "ctype.h"
#include "stdlib.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#define SQR(A) ((A)*(A))
typedef float avec[3];
#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif


const float Threshold = 0;
const float PanelGapValue = -1e10;

#include "MyMath.h"

struct rvec
{ float u;
  float v;
  float w;
};

void read_params(char* pstr, const char* separ, float* params)
{
  char *pch = strtok (pstr,"=");
  int ci = 0;
  pch = strtok (NULL, separ);
  while (pch != NULL)
  { params[ci] = atof(pch);
    ci++;
    pch = strtok (NULL, separ);
  }
}

void read_params_int(char* pstr, const char* separ, int* params)
{
  char *pch = strtok (pstr,"=");
  int ci = 0;
  pch = strtok (NULL, separ);
  while (pch != NULL)
  { params[ci] = atoi(pch);
    ci++;
    pch = strtok (NULL, separ);
  }
}

int file_type(char* fName1)
{
  int ftype = -1;
  char str1[MaxStrLen];
  if (strrchr(fName1,'.')!=NULL)
  { strcpy(str1,strrchr(fName1,'.')+1);
    strLower1(str1);
    if (strncmp(str1,"tif",3)==0 || strncmp(str1,"mccd",4)==0) ftype = 2;
    else if (strncmp(str1,"mrc",3)==0) ftype = 6;
    else if (strncmp(str1,"h5",2)==0) ftype = 5;
    else ftype = 0;
  } else ftype = 0;  //!!! not good

  return ftype;
}

bool read_file(char* fName, size_t dims[3], float** ar3D)
{   float cell[6] = {0,0,0,0,0,0};
    char comment[MaxStrLen];
    bool readOK = false;
    int ftype = file_type(fName);
    int _dims[3] = {0,0,0};

    if (ftype==2) readOK = TIFReader3D(fName, ar3D, &_dims[0], &_dims[1], &_dims[2], comment);
    //!!! HDF5!
    else if (ftype==5) 
    { printf("HDF5 is not supported yet\n");
      readOK = false; //readOK = MRCReader(fName, ar3D, &_dims[0], &_dims[1], &_dims[2], cell, comment);
    }
    else if (ftype==6) readOK = MRCReader(fName, ar3D, &_dims[0], &_dims[1], &_dims[2], cell, comment);
    else if (ftype==0) 
    {
      size_t fulldims = (size_t)_dims[0]*(size_t)_dims[1]*(size_t)_dims[2];
      if (fulldims<=0)
      { printf("For raw files you must set reasonable dimentions dims=x,y,z!\n");
        readOK = false;
      }
      size_t numEl = fulldims;
      readOK = RAWFloatReader(fName, ar3D, &numEl);
      if (numEl != fulldims) printf("Strange, read %d elements instead of expected %d\n",numEl,fulldims);
    } 
    else 
      printf("Not supported file type. Supported extentions are: tif, mrc, h5, bin/raw\n");

    for (int i=0; i<3; i++) dims[i] = (size_t)_dims[i];
    return readOK;
}  


size_t read_file_list_3D(char* fName, char** nameslist, avec** rotangles, bool* angles, bool* eulerangles, bool* scalepatterns, float** scalePats, int** frameNums)
{

  FILE* stList;
  { stList = fopen(fName,"rt");
    if (stList==NULL)
    { printf("List of files %s not found!\n",fName);
      return 0;
    }
  }

  char str1[MaxStrLen];
  size_t numFi = 0;
  while (!feof(stList))
  { if (fgets(str1,MaxStrLen,stList) != NULL)
      if (strlen(str1) > 2)
        numFi++;
  }
  rewind(stList);

  char* filesList = new char[MaxStrLen*numFi];
  *nameslist = filesList;
  avec* rotanglel = new avec[numFi];
  *rotangles = rotanglel;
  float* scalePat = new float[numFi];
  *scalePats = scalePat;
  int* frameNum = new int[numFi];
  *frameNums = frameNum;
  for (int i=0; i<numFi; i++)
  { scalePat[i] = 1;
    frameNum[i] = 0;
  }
  *scalepatterns = false;
  bool framenumpresent = false;

  for (int i=0; i<numFi; i++)
  { fgets(str1,MaxStrLen,stList);
    TrimNoCh(str1);
    char *pch = strtok (str1," \t");
    sprintf(&filesList[i*MaxStrLen],"%s\0",pch);
    pch = strtok (NULL, " \t");
    if (pch != NULL) 
    { if (pch[0]=='/' && pch[1]=='/')
      { framenumpresent = true;
        frameNum[i] = atof(&pch[2]);
        pch = strtok (NULL, " \t");
        if (pch == NULL) continue;
      }  
      rotanglel[i][0] = atof(pch) * M_PI/180.;
      *angles = true;
      pch = strtok (NULL, " \t");
      if (pch != NULL) 
      { rotanglel[i][1] = atof(pch);
        pch = strtok (NULL, " \t");
        if (pch != NULL) 
        { rotanglel[i][1] *= M_PI/180.;
          rotanglel[i][2] = atof(pch) * M_PI/180.;
          *eulerangles = true;
           pch = strtok (NULL, " \t");
           if (pch != NULL) 
           { scalePat[i] = atof(pch);
             *scalepatterns = true;
           }
        } else
        { scalePat[i] = rotanglel[i][1];
          *scalepatterns = true;
        }
      }
    }
  }

  if (!*scalepatterns)
  { delete[] scalePat;
    *scalePats = NULL;
  }

  printf("Found %ld files in the list %s %s\n",numFi,fName,(eulerangles?"with Euler angles":"with rotation angles"));
  fclose(stList);

  return numFi;
}

int main(int argc, char* argv[])
{
  printf("Postprocessing of 3D volume\n");

// TO ADD:
// cut the resulting volume - newdims and shifts OR roi
// rotate the volume (rot and angle or euler angles) -
// Also rotate with respect to "normal" vector - from paraview
// for radialSub set the center (should be an output from simerge)
// extract a slice using the origin and normal (from paraview) 


  if (argc<2)
  { printf("   postproc3d file/list [raw,tif,mrc] [dims=x,y,z] [center=x,y,z] [badval3d=X] [radialsnr3d=X] [localbg3d=+-I] \n");
//???    printf("           [averagefiles=name.lst]\n");
//    printf("           [roi=fs1,ss1,sss1,fs2,ss2,sss2] [rot=1,0,0] [angle=X] [eulerangles=X,Y,Z] [newdims=x,y,z] [shifts=x,y,z] [scale=X]\n");

//    printf("           [interp=X] [ringsmooth] [ringdiff] [badvalin=X] [badvalout=X]\n");
    printf("Files supported: raw/bin, tif, mrc\n");
    printf("localbg3d=I - local background subtraction (I - integer): +I - subtracting median calculated in I*I window, -I - replacing with median\n");
    printf("radialsnr3d=X - radial background subtraction with some SNR (usually 3-8). By defauls - off\n");
    printf("badval3d=X - the value below which the input data is masked. Default -1e10\n");
    printf("\nEXAMPLE: postproc3d file3d.tif radialsnr3d=5\n");
    return 0;
  }


  char fName[MaxStrLen] = {0};
  int convertTo = 0; //0 - raw, 1 - tif, 2 - mrc
  size_t dims[3] = {0,0,0};
  long shifts[3] = {-1,-1,-1}; 
  size_t newdims[3] = {0,0,0};
  float center[3] = {0,0,0};
  int roi[6] = {0,0,0,0,0,0}; // fs_min, ss_min, fs_max, ss_max
  float scale = 1.;
  float rot[3] = {0.,0.,0.}; //From command line, needs normalization 
  int localBG3D = 0;
  float radialSNR3D = -1;
  float outlierSNR3D = -1;
  float outliersSNR3D = -1;  // function???
  float ringDiff = -1;  //???
  int ringSmooth = 0;
  float interpRad = -1;
  float badVal3D = -10000.;
  float badValOut = -10000.;
  float ADCthresh = 0; //??? parameter?

  strcpy(fName,argv[1]);
  // analysis of the command line parameters
  for (int i=2; i<argc; i++)
  { // First we convert all keywords into lower case
    int lens = strlen(argv[i]);
    if (strchr(argv[i], '=') != NULL) lens -= strlen(strchr(argv[i], '='));
    for (int j=0; j<lens; j++)
	    argv[i][j] = tolower(argv[i][j]);

    if (strncmp(argv[i],"raw",3)==0) convertTo = 0;
    else if (strncmp(argv[i],"tif",3)==0) convertTo = 1;
    else if (strncmp(argv[i],"mrc",3)==0) convertTo = 2;
    else if (strncmp(argv[i],"dims",4)==0)
    { int _dims[3] = {0,0,0};
      read_params_int(argv[i], ",", _dims);
      for (int i=0; i<3; i++) dims[i] = _dims[i];
      printf("dims=%ld,%ld,%ld\n",dims[0],dims[1],dims[2]);
    }
    else if (strncmp(argv[i],"shifts",6)==0)
    { int _shifts[3] = {0,0,0};
      read_params_int(argv[i], ",", _shifts);
      for (int i=0; i<3; i++) shifts[i] = _shifts[i];
      printf("shifts=%ld,%ld,%ld\n",shifts[0],shifts[1],shifts[2]);
    }
    else if (strncmp(argv[i],"center",6)==0)
    { read_params(argv[i], ",", center);
      printf("center=%0.2f,%0.2f,%0.2f\n",center[0],center[1],center[2]);
    }
    else if (strncmp(argv[i],"roi",3)==0)
    { read_params_int(argv[i], ",", roi);
      printf("roi=%d,%d,%d,%d\n",roi[0],roi[1],roi[2],roi[3]);
    }
    else if (strncmp(argv[i],"rot",3)==0)
    { read_params(argv[i], ",", rot);
      printf("rot=%0.2f,%0.2f,%0.2f\n",rot[0],rot[1],rot[2]);
    }
    else if (strncmp(argv[i],"localbg3d",9)==0) localBG3D = atoi(strchr(argv[i],'=')+1);       // Those should go to postproc3d
    else if (strncmp(argv[i],"radialsnr3d",10)==0) radialSNR3D = atof(strchr(argv[i],'=')+1);
    else if (strncmp(argv[i],"outliersnr3d",11)==0) outlierSNR3D = atof(strchr(argv[i],'=')+1);
    else if (strncmp(argv[i],"scale",5)==0) scale = atof(strchr(argv[i],'=')+1);
    else if (strncmp(argv[i],"ringsmooth",10)==0) ringSmooth = atoi(strchr(argv[i],'=')+1);
    else if (strncmp(argv[i],"ringdiff",8)==0) ringDiff = atof(strchr(argv[i],'=')+1);
    else if (strncmp(argv[i],"interp",6)==0) interpRad = atof(strchr(argv[i],'=')+1);
    else if (strncmp(argv[i],"badvalin",8)==0) badVal3D = atof(strchr(argv[i],'=')+1);
    else if (strncmp(argv[i],"badval3d",8)==0) badValOut = atof(strchr(argv[i],'=')+1);
    else printf("Unknown parameter %s ! Skipping...\n", argv[i]);
  }

  // reading text file with file names and rotation angles
  size_t numFi = 0;
  char* filesList = NULL;
  avec* rotangles = NULL;
  float* scalesPat = NULL;
  int* frameNums = NULL;
  bool anglesPresent = false;
  bool scalePatterns = false;
  numFi = read_file_list_3D(fName, &filesList, &rotangles, &anglesPresent, NULL, &scalePatterns, &scalesPat, &frameNums);



    // Reading the 3D file
    float *ar3D = NULL;

    bool readOK = read_file(fName, dims, &ar3D);

    size_t num3D = dims[0]*dims[1]*dims[2];
    if (num3D<=0) readOK = false;
  
    if (!readOK)
    { printf("File %s cannot be read, exiting...\n",fName);
      return 0;
    }
  
    if (center[0]<0 || center[1]<0 || center[2]<0)
    {
      for (int i=0; i<3; i++)
        center[i] = 0.5*(dims[i]-1);
      printf("Center not set, so it is calculated: %0.2f %0.2f %0.2f\n", center[0], center[1], center[2]);
    }
  
    // Postprocessing 3D
    // Radial 3D
    if (radialSNR3D > 0)
    { int* pix_r3 = NULL;
      int maxRad = 0;
      BuildSimpleRadialArray3D(dims, &pix_r3);  //!!! WRONG! has to be around center!
      SubtractBgLoop(ar3D, num3D, pix_r3, radialSNR3D, ADCthresh, badValOut, NULL); //??? maybe save the curve?
      if (pix_r3 != NULL) delete[] pix_r3;
    }  
  
    // Local 3D
    if (localBG3D != 0)
      if (localBG3D > 0)
        SubLocalBG3D(ar3D, (size_t)dims[0], (size_t)dims[1], (size_t)dims[2], localBG3D, badValOut, NULL);
      else
      { float* smAr = new float[num3D];
        SubLocalBG3D(ar3D, (size_t)dims[0], (size_t)dims[1], (size_t)dims[2], localBG3D, badValOut, smAr);
        delete[] ar3D;
        ar3D = smAr;
      }  
  
    // Mask outliers
    if (outlierSNR3D > 0)
      MaskOutliersAtRings3D(ar3D, dims, NULL, badValOut, outlierSNR3D);
  
    char maskName[MaxStrLen];
    strcpy(maskName, fName);
    ExtractFileNameOnly(maskName);
  
    //saving 3D
    if (convertTo == 0)
    { strcat(maskName,".bin");
      SaveFloat3DtoFile(maskName, ar3D, num3D);
    } else if (convertTo == 1)
    { strcat(maskName,".tif");
      TIFWriter3D(maskName, ar3D, NULL, dims[0], dims[1], dims[2], "SimpleMerge by O.Yefanov");
    } else if (convertTo == 2)
    { strcat(maskName,".mrc");
      MRCWriterF(maskName, ar3D, dims[0], dims[1], dims[2], NULL, NULL, "SimpleMerge by O.Yefanov");
    }
  
//  printf("Cleaning\n");
    if (ar3D != NULL) delete[] ar3D;
  
  return 1;
}
//---------------------------------------------------------------------------
