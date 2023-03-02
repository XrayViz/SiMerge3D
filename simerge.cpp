// simerge.cpp
//
// Simple merging using list of frames/events with orientation angles.
//   Either Euler for each pattern or 1 angle around axis that must be set. 
//   Also possible to set region of interest.
//
// Copyright  2020 Deutsches Elektronen-Synchrotron DESY,
//                  a research centre of the Helmholtz Association.
//
// Author: Oleksandr Yefanov (oleksandr.yefanov@desy.de)

#include "stdio.h"
#include "string.h"
#include "ctype.h"
//#include <hdf5.h>    
//#include <hdf5_hl.h>
#include "stdlib.h"
//?#include <stdint.h>

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
//const int MaxStrLen = 255;

//#include "fileformats.h"
#include "readgeom.cpp"

struct rvec
{ float u;
  float v;
  float w;
};

void get_q_from_xyz_old(float rx, float ry, float rz, float k, float q[3])
{
//	float q[3];
  double r = sqrt(rx*rx + ry*ry);
  double twotheta = atan2(r, rz);
  double az = atan2(ry, rx);

  q[0] = k * sin(twotheta)*cos(az);
  q[1] = k * sin(twotheta)*sin(az);
  q[2] = k * (cos(twotheta) - 1.0);

//	return q;
};

void get_q_from_xyz_new(float rx, float ry, float rz, float k, float q[3])
{
//	float q[3];
  double L_1 = k/sqrt(rx*rx + ry*ry + rz*rz);

  q[0] = rx*L_1;
  q[1] = ry*L_1;
  q[2] = rz*L_1 - k;

//	return q;
};


void recalc_pixel_to_3D(float xx, float yy, float zz, float kwp, float trmat[9], float cqm1[])
{
//  struct rvec cqui;
  float cqm[3];
//  get_q_from_xyz(xx, yy, zz, kwp, cqm);
  get_q_from_xyz_new(xx, yy, zz, kwp, cqm);
  MultMatrixF(trmat, 3, 3, cqm, 3, 1, cqm1);
};


int add_pixel_to_3D(float cqm1[], float val, int dims[3], int shifts[3], float* ar, int* arN)
{
  int xi = roundMy(cqm1[0])+shifts[0];     // SHIFTS in pixels             
  int yi = roundMy(cqm1[1])+shifts[1];
  int zi = roundMy(cqm1[2])+shifts[2];

  if (xi>=dims[0] || yi>=dims[1] || zi>=dims[2] || xi<0 || yi<0 || zi<0) return 0;
  size_t ind = (size_t)xi + ((size_t)yi + (size_t)(zi * dims[1])) * (size_t)dims[0];
  ar[ind] += val;
  arN[ind]++;
  return 1;
}

int add_pixel_to_TPoints(float cqm1[], float val, int dims[3], int shifts[3], TPoints** allPo)
{
  int xi = roundMy(cqm1[0])+shifts[0];     // SHIFTS in pixels             
  int yi = roundMy(cqm1[1])+shifts[1];
  int zi = roundMy(cqm1[2])+shifts[2];

  if (xi>=dims[0] || yi>=dims[1] || zi>=dims[2] || xi<0 || yi<0 || zi<0) return 0;
  size_t ind = (size_t)xi + ((size_t)yi + (size_t)(zi * dims[1])) * (size_t)dims[0];

  TPoints** cPo = &allPo[ind];
  while (*cPo != NULL)
    cPo = &(*cPo)->ne;

  *cPo = new TPoints;
  (*cPo)->ne = NULL;
  (*cPo)->z = cqm1[2]+shifts[2];
  (*cPo)->y = cqm1[1]+shifts[1];
  (*cPo)->x = cqm1[0]+shifts[0];
  (*cPo)->I = val;

  return 1;
}


// The main function that recalculates a pattern into 3D
int add_pattern_to_3D(float* pat, int numcomp, avec rotangle, float rot[3], bool eulerangles, float* x, float* y, float* z, int dims[3], int shifts[3], float kwp, float badVal, float* ar, int* arN, float globTrMat[9], TPoints** allPo)
{ // here loop over all pixels in a pattern
  // DIST - in px, k - in 1/px, and scaling into k!
  int numadded = 0;
  float trmat[9];
  float trmat1[9];

  // making transformation matrix
  if (eulerangles)
    RotEulerMatrix(rotangle[0], rotangle[1], rotangle[2], trmat1);
  else
    RotAroundAxisMatrix(-rotangle[0], rot[0], rot[1], rot[2], trmat1);

  MultMatrixF(trmat1, 3, 3, globTrMat, 3, 3, trmat);

  // processing all pixels of a pattern
  for (int i=0; i<numcomp; i++)
  { if (pat[i] > badVal+1)
    { float cqm[3];
      recalc_pixel_to_3D(x[i], y[i], z[i], kwp, trmat, cqm);
      if (allPo == NULL)
        numadded += add_pixel_to_3D(cqm, pat[i], dims, shifts, ar, arN); 
      else 
        numadded += add_pixel_to_TPoints(cqm, pat[i], dims, shifts, allPo); 
    }    
  }
  
  return numadded;
}

void merge_interpolated(float interpRad, int dims[3], float* ar, TPoints** allPo, float fillWith)
{
  printf("Start merging interpolated data\n");  
  int radius = roundMy(interpRad);
#ifdef _OPENMP
  printf("OMP is on. For merge using %d threads\n",omp_get_max_threads());
#endif
#pragma omp parallel for shared(dims, ar, allPo, radius, interpRad)
  for (int zi=0; zi<dims[2]; zi++)
  {
//#ifdef _OPENMP
//	  printf("Thread %d executes zi=%d of %d\n", omp_get_thread_num(),zi,dims[2]);
//#endif
    for (int yi=0; yi<dims[1]; yi++)
     for (int xi=0; xi<dims[0]; xi++)
    { 
      size_t ind = (size_t)xi + ((size_t)yi + (size_t)(zi * dims[1])) * (size_t)dims[0];
      float sumWi = 0.; // weighting factor
      float _cI = 0.;  // intensity in a pixel

      for (int xii=xi-radius; xii<=xi+radius; xii++)
       for (int yii=yi-radius; yii<=yi+radius; yii++)
        for (int zii=zi-radius; zii<=zi+radius; zii++)
        { if (xii<0 || yii<0 || zii<0 || xii>=dims[0] || yii>=dims[1] || zii>=dims[2]) continue;
          TPoints** cPo = &allPo[xii+dims[0]*(yii+dims[1]*zii)];
          while (*cPo != NULL)
          { float _dist = sqrt(SQR((*cPo)->x-xi)+SQR((*cPo)->y-yi)+SQR((*cPo)->z-zi));
            if (_dist > interpRad)
            { cPo = &(*cPo)->ne;
              continue;
            }
            if (_dist<MinVal) // a point at this location found!
            { sumWi = 1.;
              _cI = (*cPo)->I;
              xii=dims[0]+radius; yii=dims[1]+radius; zii=dims[2]+radius;
              break; // must be up to sumWi=0
            }

            float wi = SQR((interpRad-_dist)/(interpRad*_dist));
            sumWi += wi;
            _cI += wi * (*cPo)->I;
            cPo = &(*cPo)->ne;
          }
        }
      if (sumWi>MinVal) //???
        ar[ind] = _cI/sumWi;
      else 
        ar[ind] = fillWith;
    }
  }

  printf("End merging interpolated data\n");  
}

void update_boders(float cqm[3], float bordX[2], float bordY[2], float bordZ[2])
{     if (cqm[0] < bordX[0]) bordX[0] = cqm[0];
      if (cqm[0] > bordX[1]) bordX[1] = cqm[0];
      if (cqm[1] < bordY[0]) bordY[0] = cqm[1];
      if (cqm[1] > bordY[1]) bordY[1] = cqm[1];
      if (cqm[2] < bordZ[0]) bordZ[0] = cqm[2];
      if (cqm[2] > bordZ[1]) bordZ[1] = cqm[2];
}      

void resize_dataset (int dims[3], int shifts[3], float bordX[2], float bordY[2], float bordZ[2])
{ if (dims[0] < 1 || dims[1] < 0 || dims[2] < 0)
  { dims[0] = roundMy(bordX[1] - bordX[0]) + 1;
    dims[1] = roundMy(bordY[1] - bordY[0]) + 1;
    dims[2] = roundMy(bordZ[1] - bordZ[0]) + 1;
  }
  if (shifts[0] < 0 && shifts[1] < 0 && shifts[2] < 0)
  { shifts[0] = roundMy(0.5*(dims[0]-1) - 0.5*(bordX[1] + bordX[0]));
    shifts[1] = roundMy(0.5*(dims[1]-1) - 0.5*(bordY[1] + bordY[0]));
    shifts[2] = roundMy(0.5*(dims[2]-1) - 0.5*(bordZ[1] + bordZ[0]));
  }  
}  


void calculate_boders(float roi[4], float dist, int numFi, avec* rotangles, bool eulerangles, float rot[3], float kwp, int dims[3], int shifts[3], float globTrMat[9])
{ // here check max and min for x,y,z and calculate shifts
  // here I need to estimate the shifts in 3D - go through rotations with edges of the ROI  
  // using just roi corners
  float arxr[8],aryr[8],arzr[8]; // 4 corners of ROI x1,y1,x2,y2
  for (int i=0; i<8; i++)
    arzr[i] = dist;
  // corners
  arxr[0] = arxr[3] = roi[0];
  arxr[1] = arxr[2] = roi[2];
  aryr[0] = aryr[1] = roi[1];
  aryr[2] = aryr[3] = roi[3];
  // edges
  arxr[4] = roi[0];
  arxr[6] = roi[2];
  aryr[5] = roi[1];
  aryr[7] = roi[3];
  arxr[5] = arxr[7] = 0.5*(roi[0]+roi[2]);
  aryr[4] = aryr[6] = 0.5*(roi[1]+roi[3]);

  // now finding min-max after transformation
  float bordX[2] = {1e6,-1e6};
  float bordY[2] = {1e6,-1e6};
  float bordZ[2] = {1e6,-1e6};
  float trmat1[9];
  float trmat[9];
  for (int fi=0; fi<numFi; fi++)
  { if (eulerangles)
      RotEulerMatrix(rotangles[fi][0], rotangles[fi][1], rotangles[fi][2], trmat1);
    else
      RotAroundAxisMatrix(-rotangles[fi][0], rot[0], rot[1], rot[2], trmat1);

  MultMatrixF(trmat1, 3, 3, globTrMat, 3, 3, trmat);

    // calculating borders for all corners
    for (int i=0; i<8; i++)
    { float cqm[3];
      recalc_pixel_to_3D(arxr[i], aryr[i], arzr[i], kwp, trmat, cqm);
      update_boders(cqm, bordX, bordY, bordZ);
    }  
  }

  // checking center of diffr. pattern - needed for curved Ewald sphere
  for (int fi=0; fi<numFi; fi++)
  { if (eulerangles)
      RotEulerMatrix(rotangles[fi][0], rotangles[fi][1], rotangles[fi][2], trmat1);
    else
      RotAroundAxisMatrix(-rotangles[fi][0], rot[0], rot[1], rot[2], trmat1);

    MultMatrixF(trmat1, 3, 3, globTrMat, 3, 3, trmat);

    float cqm[3];
    float cqm1[3] = {0,0,-kwp};
    MultMatrixF(trmat, 3, 3, cqm1, 3, 1, cqm);
    cqm[2] += (cqm[2]<0 ? kwp : -kwp); 
    if ((bordX[0] < cqm[0]) && (bordX[1] > cqm[0]) && (bordY[0] < cqm[1]) && (bordY[1] > cqm[1]))
    { cqm[0] = bordX[0];
      cqm[1] = bordY[0];
      update_boders(cqm, bordX, bordY, bordZ);                                      
    }
  }

  printf("borders: %0.3f - %0.3f, %0.3f - %0.3f, %0.3f - %0.3f\n",bordX[0], bordX[1], bordY[0], bordY[1], bordZ[0], bordZ[1]);
  // resize to the resulting data size:
  resize_dataset (dims, shifts, bordX, bordY, bordZ);

}

void calculate_boders2(int* gp, float* x, float* y, float* z, int numcomp, int numFi, avec* rotangles, bool eulerangles, float rot[3], float kwp, int dims[3], int shifts[3], float globTrMat[9])
{ // here check ALL x,y,z and calculate shifts
  printf("Calculating exact borders for the resulting 3D. Might take long.\n");
  // now finding min-max after transformation
  float bordX[2] = {1e6,-1e6};
  float bordY[2] = {1e6,-1e6};
  float bordZ[2] = {1e6,-1e6};
  float trmat1[9];
  float trmat[9];
  for (int fi=0; fi<numFi; fi++)
  { if (eulerangles)
      RotEulerMatrix(rotangles[fi][0], rotangles[fi][1], rotangles[fi][2], trmat1);
    else
      RotAroundAxisMatrix(-rotangles[fi][0], rot[0], rot[1], rot[2], trmat1);
    
    MultMatrixF(trmat1, 3, 3, globTrMat, 3, 3, trmat);
    
    for (int i=0; i<numcomp; i++)
    { float cqm[3];
      recalc_pixel_to_3D(x[i], y[i], z[i], kwp, trmat, cqm);
      update_boders(cqm, bordX, bordY, bordZ);
    }  
  }

  // resize to the resulting data size:
  resize_dataset (dims, shifts, bordX, bordY, bordZ);

}


void applying_roi(float roixy[4], int roi[4], int numfs, int numss, float* arx, float* ary, int* gp)
{ int numcomp = numfs*numss;
  // Applying ROI
  if ((roi[2]-roi[0])>0 && (roi[3]-roi[1])>0)
    for (int ssi=0; ssi<numss; ssi++)
      for (int fsi=0; fsi<numfs; fsi++)
        if (!(fsi >= roi[0] && fsi <= roi[2] && ssi >= roi[1] && ssi <= roi[3]))
          gp[fsi + numfs*ssi] = -1;
  // Applying ROIXY
  if ((roixy[2]-roixy[0])>0 && (roixy[3]-roixy[1])>0)
    for (int i=0; i<numcomp; i++)
      if (!(arx[i] >= roixy[0] && arx[i] <= roixy[2] && ary[i] >= roixy[1] && ary[i] <= roixy[3]))
        gp[i] = -1;
}

void clean_good_pixels(int* numcomp, float** arxR, float** aryR, float** arzR, int** gpR)
{
  float* arx = *arxR;
  float* ary = *aryR;
  float* arz = *arzR;
  int* gp = *gpR;

    float* ary1 = new float[*numcomp];
    float* arx1 = new float[*numcomp];
    float* arz1 = new float[*numcomp];  
    int* gp1 = new int[*numcomp];  
    int nc1 = 0;
    for (int i=0; i<*numcomp; i++)
      if (gp[i]>=0)                    // to take into account pixels already disabled by mask
      { gp1[nc1] = gp[i];
        arx1[nc1] = arx[gp[i]];
        arz1[nc1] = arz[gp[i]];
        ary1[nc1] = ary[gp[i]];
        nc1++;
      }
    delete[] arx;
    delete[] ary;
    delete[] arz;
    delete[] gp;
    *numcomp = nc1;
    float* arx2 = new float[*numcomp];
    float* ary2 = new float[*numcomp];
    float* arz2 = new float[*numcomp];  
    int* gp2 = new int[*numcomp];  
    *arxR = arx2;
    *aryR = ary2;
    *arzR = arz2;  
    *gpR = gp2;  
    for (int i=0; i<*numcomp; i++)
    { arx2[i] = arx1[i];
      ary2[i] = ary1[i];
      arz2[i] = arz1[i];
      gp2[i] = gp1[i];
    }
    delete[] arx1;
    delete[] ary1;
    delete[] arz1;
    delete[] gp1;

}

void update_rois(float roixy[4], int roi[4], int numcomp, float* arx, float* ary, float* arz, int* gp, int numfs)
{
    // new roixy
    roixy[0] = roixy[1] = 1e10;
    roixy[2] = roixy[3] = -1e10;
    for (int i=0; i<numcomp; i++)
    { if (arx[i]<roixy[0]) roixy[0] = arx[i];
      if (ary[i]<roixy[1]) roixy[1] = ary[i];
      if (arx[i]>roixy[2]) roixy[2] = arx[i];
      if (ary[i]>roixy[3]) roixy[3] = ary[i];
    }
    // new roi
    roi[0] = roi[1] = 1000000;
    roi[2] = roi[3] = 0;
    for (int i=0; i<numcomp; i++)
    { int cfs = gp[i] % numfs;
      int css = gp[i] / numfs;
      if (cfs<roi[0]) roi[0] = cfs;
      if (css<roi[1]) roi[1] = css;
      if (cfs>roi[2]) roi[2] = cfs;
      if (css>roi[3]) roi[3] = css;
    }
}

size_t read_file_list(char* fName, char** nameslist, avec** rotangles, bool* angles, bool* eulerangles, bool* scalepatterns, float** scalePats, int** frameNums)
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

  // just testing
//  stList = fopen("out.txt","wt");
//  for (int i=0; i<numFi; i++)
//    if (i%10 == 0) 
//    { strncpy(str1,&filesList[i*MaxStrLen],MaxStrLen);
//      fprintf(stList,"%s\t//%d\t%0.4f\t%0.5f\n",str1, frameNum[i], rotanglel[i][0]*180./M_PI, scalePat[i]);
//    }  
//  fclose(stList);

  return numFi;
}


int file_type(char* fName1)
{
  int ftype = -1;
  char str1[MaxStrLen];
  if (strrchr(fName1,'.')!=NULL)
  { strcpy(str1,strrchr(fName1,'.')+1);
    strLower1(str1);
    if (strncmp(str1,"cbf",3)==0) ftype = 1;
    else if (strncmp(str1,"tif",3)==0 || strncmp(str1,"mccd",4)==0) ftype = 2;
    else if (strncmp(str1,"img",3)==0) ftype = 3;
    else if (strncmp(str1,"edf",3)==0) ftype = 4;
    else if (strncmp(str1,"h5",2)==0 || strncmp(str1,"cxi",3)==0) ftype = 5;
    else ftype = 0;
//      if (strncmp(str1,"raw",3)==0 || strncmp(str1,"bin",3)==0) ftype = 0;
//    else
//    { printf("File extension %s is not recognized - skipping\n",str1);
//      ftype = -1;
//    }
  } else ftype = 0;  //!!! not good

  return ftype;
}

int read_a_file(char* fName1, float** apattern, int* dimsL, int numfs, int numss, bool transpose, char* datafield, int frameaxis, int framenum, float* angle)
{
//    int ftype = 0; //0 - raw, 1 - cbf, 2 - tif, 3 - img, 4 - edf, 5 - h5???
  int ftype = file_type(fName1);
  char comment[MaxStrLen];
  bool readOK = false;

  float angleInc = 0;
  float pixel[2] = {0,0};
  float expo = 0;
  float waveLen = 0;
  float dist = 0;
  float beamxy[2] = {0,0};
  float flux = 0;
//  float angle = 0;

  if (ftype==1) readOK = ReadCBFfile(fName1, apattern, dimsL, pixel, &expo, &waveLen, &dist, beamxy, &flux, NULL,  angle, &angleInc,PanelGapValue);
  else if (ftype==2) readOK = TIFReader(fName1, apattern, &dimsL[0], &dimsL[1], comment);
  else if (ftype==3) readOK = ReadIMGfile(fName1, apattern, dimsL, pixel, &expo, &waveLen, &dist, beamxy, &flux, NULL, angle, &angleInc);
  else if (ftype==4) readOK = ReadEDFfile(fName1, apattern, dimsL, pixel, &expo, &waveLen, &dist, beamxy, &flux, NULL, angle, &angleInc);
  else if (ftype==5) readOK = ReadHDF5frame2D(fName1, datafield, apattern, frameaxis, framenum, dimsL);
  else if (ftype==0) 
  {
    if (!(numfs>0 && numss>0))
    { printf("For raw files you must set dimentions in geometry (fs,ss)!\n");
      readOK = false;
    }
    int fulldims = numfs*numss;
    size_t numEl = (size_t)fulldims;
    dimsL[0] = numfs;
    dimsL[1] = numss;
    readOK = RAWFloatReader(fName1, apattern, &numEl);
    if (numEl != fulldims) printf("Strange, read %d elements instead of expected %d\n",numEl,fulldims);
  }

  if (readOK)
  { //??? shell be done in geometry, I think
    if (transpose) TransposeArray(&dimsL[0], &dimsL[1], *apattern);  
    return dimsL[0]*dimsL[1];
  } else
    return 0;

}

int read_mask(char* maskName, int* gp, int numfs, int numss, bool transpose, char* maskfield)
{   float* amask = NULL;
    int maskfs=numfs;
    int maskss=numfs;
    int numcomp = numfs*numss;
    int dimsL[3] = {0,0,1};
    int numBads = 0;
    float angle;
    if (read_a_file(maskName, &amask, dimsL, maskfs, maskss, transpose, maskfield, 0, 0, &angle) == numcomp)
    { 
        for (int i=0; i<numcomp; i++)
          if (amask[i]<1) 
          { gp[i] = -1;
            numBads++;
          };
//      printf("Mask %s was read, number of bad pixels: %d\n", maskName, numBads);
    } else 
      printf("Mask can't be read, number of pixels in the file %d*%d != %d (number in geometry). Skipping the mask!\n", dimsL[0], dimsL[1], numcomp);
    if (amask != NULL)
      delete[] amask;
  return numBads;
}

void read_dark_gain(char* fileName, int numfs, int numss, bool transpose, char* datafield, const char* fitype, float* anarray)
{   float anangle;
    int dimsL[3] = {0,0,1};
    size_t numcomp = (size_t)numfs*(size_t)numss;
    if (read_a_file(fileName, &anarray, dimsL, numfs, numss, transpose, datafield, 0, 0, &anangle) != numcomp)
    { 
      printf("%s can't be read or number of pixels in the file %d*%d != %d (number in geometry). Not using %s!\n", fitype, dimsL[0], dimsL[1], numcomp, fitype);
      if (anarray != NULL) delete [] anarray;
      anarray = NULL;
    } else 
      printf("The %s file %s will be subtracted from each image.\n", fitype, fileName);
}


void save_first_pat(char* fName1, int dimsL[2], int numcomp, int* gp, float* apattern)
{
        int numcompFull = dimsL[0]*dimsL[1];
        float* roipattern = new float[numcompFull];
        for (int i=0; i<numcompFull; i++)
          roipattern[i] = 0.;
        for (int i=0; i<numcomp; i++)
          roipattern[gp[i]] = apattern[i];
        sprintf(fName1,"%s_ROI.tif",ExtractJustFileNameNoExtC(fName1));
        printf("Saving an example file %s with applied ROI, npix=%d\n", fName1, numcomp);
        TIFFloatWriter(fName1, roipattern, dimsL[1], dimsL[0], "roi");
        delete[] roipattern;

}

void save_2D(bool convertTo, char* fName, float* ar, int dim0, int dim1)
{
  char mName[MaxStrLen] = {0};
  size_t num = (size_t)dim0*(size_t)dim1;
  strcpy(mName, fName);
  if (convertTo == 0)
  { strcat(mName,".bin"); //!!! need to add dimentions in the file name
    SaveFloat3DtoFile(mName, ar, num);
  } else if (convertTo == 1)
  { strcat(mName,".tif");
    TIFFloatWriter(mName, ar, dim0, dim1, "SimpleMerge by O.Yefanov");
  } else if (convertTo == 2)
  { strcat(mName,".mrc");
    MRCWriterF(mName, ar, dim0, dim1, 1, NULL, NULL, "SimpleMerge by O.Yefanov");
  } else if (convertTo == 3)
  { strcat(mName,".h5");
    int dims[2] = {dim0, dim1};
    WriteHDF5(mName, "/data/data", ar, dims, 2);
  }
}  

void save_3D(bool convertTo, char* fName, float* ar3D, size_t num3D, int dims[3])
{
  char mName[MaxStrLen] = {0};
  strcpy(mName, fName);
  if (convertTo == 0)
  { strcat(mName,".bin"); //!!! need to add dimentions in the file name
    SaveFloat3DtoFile(mName, ar3D, num3D);
  } else if (convertTo == 1)
  { strcat(mName,".tif");
    TIFWriter3D(mName, ar3D, NULL, dims[0], dims[1], dims[2], "SimpleMerge by O.Yefanov");
  } else if (convertTo == 2)
  { strcat(mName,".mrc");
    MRCWriterF(mName, ar3D, dims[0], dims[1], dims[2], NULL, NULL, "SimpleMerge by O.Yefanov");
  } else if (convertTo == 3)
  { strcat(mName,".h5");
    WriteHDF5(mName, "/data/data", ar3D, dims, 3);
  }
}  

void save_3D_and_Slices(bool convertTo, char* fName, float* ar3D, size_t num3D, int dims[3])
{
  char fOutName[MaxStrLen];
  size_t dim[3] = {dims[0],dims[1],dims[2]};
  float* slice3Dxy = new float[dim[0]*dim[1]];
  for (size_t yi=0; yi<dim[1]; yi++)
    for (size_t xi=0; xi<dim[0]; xi++)
    { size_t zi = roundMy((dims[2]-1)*0.5);
      size_t _ind = xi+dim[0]*(yi+dim[1]*zi);
		  slice3Dxy[yi*dims[0]+xi] = ar3D[_ind];
    }  
  sprintf(fOutName,"%s_slice_xy",fName);//,numPo);
  save_2D(convertTo, fOutName, slice3Dxy, dim[0], dim[1]);
  delete[] slice3Dxy;

  float* slice3Dyz = new float[dim[1]*dim[2]];
  for (size_t zi=0; zi<dim[2]; zi++)
    for (size_t yi=0; yi<dim[1]; yi++)
    { size_t xi = roundMy((dims[0]-1)*0.5);
      size_t _ind = xi+dim[0]*(yi+dim[1]*zi);
		  slice3Dxy[zi*dims[1]+yi] = ar3D[_ind];
    }  
  sprintf(fOutName,"%s_slice_yz",fName);//,numPo);
  save_2D(convertTo, fOutName, slice3Dxy, dim[1], dim[2]);
  delete[] slice3Dyz;

  float* slice3Dzx = new float[dim[2]*dim[0]];
  for (size_t xi=0; xi<dim[0]; xi++)
    for (size_t zi=0; zi<dim[2]; zi++)
    { size_t yi = roundMy((dims[1]-1)*0.5);
      size_t _ind = xi+dim[0]*(yi+dim[1]*zi);
		  slice3Dxy[xi*dims[2]+zi] = ar3D[_ind];
    }  
  sprintf(fOutName,"%s_slice_zx",fName);//,numPo);
  save_2D(convertTo, fOutName, slice3Dzx, dim[2], dim[0]);
  delete[] slice3Dzx;

  save_3D(convertTo, fName, ar3D, num3D, dims);

}

void apply_lb(int numfs, int numss, int* roi, int localBG, float badVal1, float* apattern)
    { float badVal = -1e10; // change if the mask is applied before        //????????????????????
      int mar[4] = {0,0,numfs,numss};
      if (true) //??? to CHECK!
      { mar[0] = (roi[0]-localBG>0?roi[0]-localBG:0);
        mar[1] = (roi[1]-localBG>0?roi[1]-localBG:0);
        mar[2] = (roi[2]+localBG<numfs?roi[2]+localBG:numfs);
        mar[3] = (roi[3]+localBG<numss?roi[3]+localBG:numss);
      }
      if (localBG > 0)
        SubLocalBG(apattern, mar[0], mar[2], mar[1], mar[3], numfs, localBG, localBG, badVal, NULL);
//        SubLocalBG(apattern, 0, numfs, 0, numss, numfs, localBG, localBG, badVal, NULL);
      else
      { badVal = -10000;   //????????????????
        float* smAr = new float[numfs*numss];
        SubLocalBG(apattern, mar[0], mar[2], mar[1], mar[3], numfs, localBG, localBG, badVal, smAr);
//        SubLocalBG(apattern, 0, numfs, 0, numss, numfs, localBG, localBG, badVal, smAr);
        delete[] apattern;
        apattern = smAr;
      }
    }

void only_gp(int* fulldims, int numcomp, float* apattern, int* gp)
{ 
    if (*fulldims > numcomp) // need truncation using gp array
    { float* pat1 = new float[numcomp];
      for (int i=0; i<numcomp; i++)
        pat1[i] = apattern[gp[i]];
      delete[] apattern;
      apattern = pat1;  
    }
}    

void mask_some_values(float* apattern, int numcomp, float lessthan, float morethan, float equalto, float badVal)
{
  for (int i=0; i<numcomp; i++)
  { if (apattern[i] > morethan) apattern[i] = badVal;
    if (apattern[i] < lessthan) apattern[i] = badVal;
    if (apattern[i] == equalto) apattern[i] = badVal;
  }  
}

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

void SaveExampleGeometry()
{ 
    printf("Geometry couldn't be read!\n"); 
    printf("Please make a geometry file.\n");
    printf("An example geometry for a simple detector with 400x400 pixels, 75um each at 0.1m distance is saved as example.geom \n"); 
    FILE* ageom = fopen("example.geom","wt");
    fprintf(ageom,"clen = 0.1 ; Detector distance in m\n");
    fprintf(ageom,"res = 13333.3 ; Inverse pixel size in 1/m, (here 1/75um) \n\n");
    fprintf(ageom,"p0/min_fs = 0\n");
    fprintf(ageom,"p0/min_ss = 0\n");
    fprintf(ageom,"p0/max_fs = 399\n");
    fprintf(ageom,"p0/max_ss = 399\n");
    fprintf(ageom,"p0/corner_x = -199.5 ; Detector corner in pixels (== -Det_Center)\n");
    fprintf(ageom,"p0/corner_y = -199.5\n");
    fprintf(ageom,"p0/fs = +1.0x +0.0y ; Direction of fast and slow scan axis\n");
    fprintf(ageom,"p0/ss = +0.0x +1.0y\n");
    fclose(ageom);
    printf("Detailed description and more examples for existing detectors can be found at:  \n");
    printf("https://www.desy.de/~twhite/crystfel\n");
}

void ShowHelp()
{
    printf("  simerge list geometry=geom.geom [rot=1,0,0] [dims=x,y,z] [shifts=x,y,z]\n");
    printf("          [scale/voxelsize=X] [polariz=X] [nocorrectsolid] [interp=X]\n");
    printf("          [localbg=+-I] [radialsnr=X] [ringsmooth] [ringdiff] [badvalin=X]\n");
    printf("          [maskroi=file.ext] [roi=fs1,ss1,fs2,ss2] [roixy=x1,y1,x2,y2]\n");
    printf("          [mask=file.ext] [maskfield=xxx] [dark=dark.ext] [gain=gain.ext]\n");
    printf("          [euler|newnorm=x,y,z] [exactborders] [savenums] [raw,tif,mrc] \n");
//postproc    printf("           [badval3d=X] [radialsnr3d=X] [localbg3d=X] \n");
    printf("Files supported: raw/bin, tif, cbf, edf, img, h5\n");
    printf("  localbg=I - local background subtraction (I - integer): +I - subtracting median calculated in I*I window, -I - replacing with median\n");
    printf("  radialSNR=X - radial background subtraction with some SNR (usually 3-8). By defauls - off\n");
    printf("  polariz=X - polarization correction, 1 is vertical, 0 - horizontal and 0-1 - anything in between, -1 - no correction\n");
    printf("  interp=X - merge with interpolation within radius of X voxels\n");
    printf("  exactborders - calculate exact bsize and shifts of the resulting 3D array. Long...\n");
    printf("  badvalin=X - the value below which the input data is masked. Default -1e10\n");
    printf("  badvalout=X - the value to fill not-existing regions in 3D. Can be 0, but not recommended. Default -10000\n");
    printf("  dark=dark.ext and gain=gain.ext - using a DARK file (subtracting) and a GAIN file (multiplying) for each of the input files\n");
    printf("  savenums - saves the 3D array with number of pixels per voxel\n");
    printf("\nEXAMPLE: simerge files.lst geometry=example.geom rot=0,1,0 scale=0.5\n");
  
}


int main(int argc, char* argv[])
{
  printf("Does simple merge\n");

  if (argc<2)
  { ShowHelp();
    return 0;
  }

  //!!! to ADD:
  // Maybe: SimpleDeleteOutliers, 
  // Wedge??? I can try to estimate to how many pixels each radius has to be split (divide by num and add to each voxel). Only nearest neighbour
  // get back the 2nd radius for the interpolation
  // if internal mask is present (hdf5 file, specified in the geometry) - apply it in addition to other masks
  // output the shift of the 3D center and the center of rotation for further analysis (in pixels x,y,z)

  // separate tool for 3d postproc
  // separate tool for visualization (VTK)

  //!!! to CHANGE:
  // maybe remove transpose - it should be in the geometry (x,y directions)
  // remove all 3D postprocessing - into a separate program (postproc3d)

  // BUGs !!!!!!!!
  
  //!!! to TEST: 
  // CorrectSolidAngle, CorrectPolarisation
  // localBG,  MaskRingsSimple 
  // POSTPROCESSING in 3D: SubLocalBG3D, radial, MaskOutliersAtRings3D
  // take angle from the header, also print at the screen in this case (or log?)
  // LB do in the box bigger than ROI by LB radius
  // Make a possibility to do radial before applying the ROI (negative SNR)
  // MaxVal and MinVal - probably in geometry. And take into account for "save_first_file" (flag_lessthan, flag_morethan, flag_equal)
  // additional rotation (both Euler angles and )

  //!!! TESTed:
  // Scale for each pattern from the file - add one more collumn
  // INTERPOLATION,  radialBG
  // new function for Qx,y,z

  //!!! Into description
  // non-perpendicular to the beam detectors are possible, but the function ReadGeometry has to be modified. Look at the corresponding function in CrystFel, better version 0.8.0. Also the fast extimation of dims and shifts might not work - use the precise one

  // input: file.lst, [geometry=XXX], [dims=X,Y,Z out size (default - sbigest size of ROI / scale)], [roi=x1,x2,y1,y2 in px but x,y (not fs,ss)], [scale=XXX in px], [rot=x,xy,z axis], [distance=XXX in pixels!] [save as raw/tif?/mrc]
  // read a file and geometry (if no geom - recalculate for the center)
  // function to recalc the Q (from CF)
  // function to merge 

  char fName[MaxStrLen] = {0};
  char maskName[MaxStrLen] = {0};                    
  char maskROI[MaxStrLen] = {0}; //a "positive" mask
  int convertTo = 0; //0 - raw, 1 - tif, 2 - mrc; 3 - hdf5
  int dims[3] = {0,0,0};
  int shifts[3] = {-1,-1,-1}; 
  int roi[4] = {0,0,0,0}; // fs_min, ss_min, fs_max, ss_max
  float roixy[4] = {0.,0.,0.,0.}; // x_min,y_min,x_max,y_max in px
  float scale = 1.;
  float voxelsize = 0.;
  float rot[3] = {0.,0.,0.}; //From command line, needs normalization 
  float newRot[3] = {0.,0.,0.};
//  float newNorm[3] = {0.,0.,0.};
  int addRot = 0; // 1 - Euler angles, 2 - new normal
  char geometry[MaxStrLen] = {0};
  char datafield[MaxStrLen] = {0};
  char maskfield[MaxStrLen] = "/data/data";
  char darkName[MaxStrLen] = {0};
  char gainName[MaxStrLen] = {0};
  int frameaxis = -1; // along which axis in HDF5 data is stacked
  float distance = 0; //in pixels!  Not sure that it's needed - I din't make a function to merge simple data
  bool transpose = false;  //Some images hade to be transposed
  int localBG = 0;
  float radialSNR = 0;
  float outliersSNR = -1;  // function???
  float ringDiff = -1;
  int ringSmooth = 0;
  float interpRad = -1;
  int localBG3D = 0;
  float radialSNR3D = 0;
  float outlierSNR3D = -1;
  float polariz = -1.;
  float badVal = -10000.;
  float badValOut = -10000.;
  bool correctsolid = true;
  bool exactborders = false;
  bool saveNums = false;
  int debug = 0;

  float ADCthresh = 0; //??? parameter?
  bool eulerangles = false;
  char* pch;

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
    else if (strncmp(argv[i],"h5",2)==0) convertTo = 3;
    else if (strncmp(argv[i],"dims",4)==0)
    { read_params_int(argv[i], ",", dims);
      printf("dims=%d,%d,%d\n",dims[0],dims[1],dims[2]);
    }
    else if (strncmp(argv[i],"shifts",6)==0)
    { read_params_int(argv[i], ",", shifts);
      printf("shifts=%d,%d,%d\n",shifts[0],shifts[1],shifts[2]);
    }
    else if (strncmp(argv[i],"roixy",5)==0)
    { read_params(argv[i], ",", roixy);
      printf("roixy=%0.2f,%0.2f,%0.2f,%0.2f\n",roixy[0],roixy[1],roixy[2],roixy[3]);
    }
    else if (strncmp(argv[i],"roi",3)==0)
    { read_params_int(argv[i], ",", roi);
      printf("roi=%d,%d,%d,%d\n",roi[0],roi[1],roi[2],roi[3]);
    }
    else if (strncmp(argv[i],"rot",3)==0)
    { read_params(argv[i], ",", rot);
      printf("rot=%0.2f,%0.2f,%0.2f\n",rot[0],rot[1],rot[2]);
    }
    else if (strncmp(argv[i],"euler",5)==0)
    { read_params(argv[i], ",", newRot);
      addRot = 1;
      printf("Euler angles=%0.2f,%0.2f,%0.2f deg\n",newRot[0],newRot[1],newRot[2]);
    }
    else if (strncmp(argv[i],"newnorm",7)==0)
    { read_params(argv[i], ",", newRot);
      addRot = 2;
      printf("Make new normal=%0.2f,%0.2f,%0.2f deg\n",newRot[0],newRot[1],newRot[2]);
    }
    else if (strncmp(argv[i],"scale",5)==0) scale = atof(strchr(argv[i],'=')+1);
    else if (strncmp(argv[i],"voxelsize",8)==0) voxelsize = atof(strchr(argv[i],'=')+1);
    else if (strncmp(argv[i],"polariz",7)==0) polariz = atof(strchr(argv[i],'=')+1);
    else if (strncmp(argv[i],"transpose",8)==0) transpose = true;
    else if (strncmp(argv[i],"nocorrectsolid",12)==0) correctsolid = false;
    else if (strncmp(argv[i],"exactborders",10)==0) exactborders = true;
    else if (strncmp(argv[i],"savenums",7)==0) saveNums = true;
    else if (strncmp(argv[i],"distance",8)==0) distance = atof(strchr(argv[i],'=')+1);
    else if (strncmp(argv[i],"localbg",7)==0) localBG = atoi(strchr(argv[i],'=')+1);
    else if (strncmp(argv[i],"radialsnr",8)==0) radialSNR = atof(strchr(argv[i],'=')+1);
    else if (strncmp(argv[i],"ringsmooth",10)==0) ringSmooth = atoi(strchr(argv[i],'=')+1);
    else if (strncmp(argv[i],"ringdiff",8)==0) ringDiff = atof(strchr(argv[i],'=')+1);
    else if (strncmp(argv[i],"interp",6)==0) interpRad = atof(strchr(argv[i],'=')+1);
    else if (strncmp(argv[i],"badvalin",8)==0) badVal = atof(strchr(argv[i],'=')+1);
    else if (strncmp(argv[i],"badval3d",8)==0) badValOut = atof(strchr(argv[i],'=')+1);
    else if (strncmp(argv[i],"geometry",8)==0) strcpy(geometry,(strchr(argv[i],'=')+1));
    else if (strncmp(argv[i],"maskroi",7)==0) strcpy(maskROI,(strchr(argv[i],'=')+1));
    else if (strncmp(argv[i],"maskfield",9)==0) strcpy(maskfield,(strchr(argv[i],'=')+1));
    else if (strncmp(argv[i],"mask",4)==0) strcpy(maskName,(strchr(argv[i],'=')+1));
    else if (strncmp(argv[i],"dark",4)==0) strcpy(darkName,(strchr(argv[i],'=')+1));
    else if (strncmp(argv[i],"gain",4)==0) strcpy(gainName,(strchr(argv[i],'=')+1));
    else if (strncmp(argv[i],"debug",5)==0) debug = atoi(strchr(argv[i],'=')+1);       // Those should go to postproc3d
    else if (strncmp(argv[i],"localbg3d",9)==0) localBG3D = atoi(strchr(argv[i],'=')+1);       // Those should go to postproc3d
    else if (strncmp(argv[i],"radialsnr3d",10)==0) radialSNR3D = atof(strchr(argv[i],'=')+1);
    else if (strncmp(argv[i],"outliersnr3d",11)==0) outlierSNR3D = atof(strchr(argv[i],'=')+1);
    else printf("Unknown parameter %s ! Skipping...\n", argv[i]);
  }
  // some default values
  if (strlen(geometry) < 1 && distance < 0.1)
  { 
    SaveExampleGeometry();
    printf("You have to provide a geometry file! Look at an example geometry \"example.geom\". Exiting...\n");
    return 0;
  }

  // reading text file with file names and rotation angles
  size_t numFi = 0;
  char* filesList = NULL;
  avec* rotangles = NULL;
  float* scalesPat = NULL;
  int* frameNums = NULL;
  bool anglesPresent = false;
  bool scalePatterns = false;
  numFi = read_file_list(fName, &filesList, &rotangles, &anglesPresent, &eulerangles, &scalePatterns, &scalesPat, &frameNums);
  if (numFi<1)
  { printf("No files found!\n");
    return 0;
  }
  if (!anglesPresent)
    printf("Rotation angles were not found, we'll try to read from the actual files (cbf, img, edf)\n");
  
  // rotation axis
  if (fabs(rot[0])<0.01 && fabs(rot[1])<0.01 && fabs(rot[2])<0.01)
  { if (!eulerangles)
      printf("Rotation axis is not set, so 1,0,0 will be used.\n");
    rot[2] = rot[1] = 0.;
    rot[0] = 1.;  
  } else 
  { if (eulerangles)
      printf("File contains Euler angles, so rotation axis won't be used.\n");
    float rotL = sqrt(SQR(rot[0])+SQR(rot[1])+SQR(rot[2]));
    for (int i=0; i<3; i++)
      rot[i] /= rotL;
  }  

  if (scalePatterns)
    printf("will rescale patterns\n");


  float *arx,*ary,*arz;
  int *gp; // good pixels that have to be considered
  float istep = 0;
  float kwp_real = 0;
  float kwp = 0;
  int numcomp = 0;
  int numfs = 0;
  int numss = 0;
  float lessthan = 1e37;
  float morethan = -1e37;
  float equalto = -1e37;

  // read the geometry
  bool geomread = false;
  char evmaskfield[MaxStrLen] = {0};
  if (strlen(geometry) > 1)
    geomread = ReadGeometry(geometry, &arx, &ary, &arz, &istep, &numfs, &numss, &numcomp, &kwp_real, datafield, evmaskfield, &frameaxis, &lessthan, &morethan, &equalto);
  if (!geomread || fabs(arz[0]) < MinVal)  
  {
    SaveExampleGeometry();
    if (!geomread)
      printf("Geometetry couldn't be read.\n"); 
    else 
      printf("Detector distance is 0. Check \"clen\" and \"coffset\" parameters.\n"); 
    printf("Look at an example geometry \"example.geom\". Exiting...\n"); 
    return 0;
  }
  int numcompFull = numcomp;   // just to have it for creation of 3D array

  //??? if the geometry is not set - simple square geometry from a file (need to load one).  NOT SURE IF NEEDED

  //good pixels
  gp = new int[numcomp];  
  for (int i=0; i<numcomp; i++)
  { arx[i] *= istep;
    ary[i] *= istep;
    arz[i] *= istep;
    gp[i] = i;
  }

  // rescale geometry
//???  kwp *= scale/istep;
  kwp = arz[0]*scale; 

  float single_voxel = 10000*kwp_real/kwp; //20.*(kwp_real*sin(0.5*atan(1./kwp))); //in nm^-1 
//  printf("Kwp = %0.3f, kwp_real=%0.8f\n", kwp, kwp_real);
  printf("Taking into account the energy E=%0.3fkeV in the geometry file, the voxel size in 3D will be: %0.6fum^-1\n", kwp_real*12.3984, single_voxel);
  // if the voxelsize is given, calculate kwp from it:
  if (voxelsize > MinVal)
    if (kwp_real > MinVal)
    { kwp = 10000*kwp_real/voxelsize;
      printf("Voxel size of %0.3fnm^-1 is used, so scale factor is %0.3f\n", voxelsize, kwp/arz[0]);
    } 
    else 
      printf("ACHTUNG!!! If you want voxelsize please specify the energy in the geometry. Now using voxel size calculated before\n");  

  // If RadialSNR <0, do radial subtraction for FULL patterns:
  // Building an array with radii for each pixel
  int* pix_r = NULL;
  int maxRad = 0;
  if (radialSNR < -MinVal)
  { pix_r = new int[numcomp];
    float* pixelsR = new float[numcomp];
    BuildRadialArray(numcomp, arx, ary, 1, pix_r, &maxRad, pixelsR);
    if (pixelsR != NULL) delete[] pixelsR;
  } 

  // Reading gain and/or dark files if present
  float* adark = NULL;
  if (strlen(darkName)>1)
    read_dark_gain(darkName, numfs, numss, transpose, datafield, "DARK", adark);  
  float* again = NULL;
  if (strlen(gainName)>1)
    read_dark_gain(gainName, numfs, numss, transpose, datafield, "GAIN", again);  

  // Reading masks, if present (I can also take it into account in the gp array - remake arXYZ and change numcomp)
  int numbads = 0;
  if (strlen(maskName)>1)
    numbads = read_mask(maskName, gp, numfs, numss, transpose, maskfield);
  if (strlen(maskROI)>1)
    numbads = read_mask(maskROI, gp, numfs, numss, transpose, maskfield);

  // Applying ROI and ROIXY
  applying_roi(roixy, roi, numfs, numss, arx, ary, gp);
 
  // Getting rid of all pixels that are not interesting and finding new roixy
  clean_good_pixels(&numcomp, &arx, &ary, &arz, &gp);

  // Updating the ROIs (roi is used for LB, roixy - for 3D size determination)
  update_rois(roixy, roi, numcomp, arx, ary, arz, gp, numfs);

  float globTrMat[9] = {1.,0.,0.,0.,1.,0.,0.,0.,1.};
  if (addRot > 0)
  { if (addRot == 1) //Euler angles
      RotEulerMatrix(newRot[0]*M_PI/180., newRot[1]*M_PI/180., newRot[2]*M_PI/180.,  globTrMat);
    else if (addRot == 2)
    { //??? need to understand how to make rotation matrix in this case
      float cqm1[3] = {1.,0.,0.};
      float nrabs = sqrt(SQR(newRot[0])+SQR(newRot[1])+SQR(newRot[2]));
      for (int i=0; i<3; i++)
        newRot[i] /= nrabs;
      MultMatrixF(cqm1, 3, 1, newRot, 1, 3, globTrMat);
    } 
    else addRot = 0;
  } 

  // here I need to estimate the shifts in 3D and the size of 3D if not set
  // 1st method: go through rotations with edges of the ROIXY and the most distant Ewald sphere point (only 1 det distance is used)
  // 2nd method - just calculate Q for each pixel - precise but slow
  if (dims[0]<1 || dims[1]<1 || dims[2]<1 || shifts[0]<0 || shifts[1]<0 || shifts[2]<0)
    if (exactborders) calculate_boders2(gp, arx, ary, arz, numcomp, numFi, rotangles, eulerangles, rot, kwp, dims, shifts, globTrMat);
    else calculate_boders(roixy, arz[0], numFi, rotangles, eulerangles, rot, kwp, dims, shifts, globTrMat);

  // Choosing the size of the resulting array
  size_t num3D = (size_t)dims[0]*(size_t)dims[1]*(size_t)dims[2];
  if (num3D < 1) // I think it's not needed already
  { dims[0] = dims[1] = dims[2] = roundMy(scale*sqrt(numcompFull));
    num3D = (size_t)dims[0]*(size_t)dims[1]*(size_t)dims[2];
    printf("Generating output array of the size %dx%dx%d pixels\n",dims[0],dims[1],dims[2]);
  } else
    printf("Size of the 3D array is %dx%dx%d pixels, shifts: %dx%dx%d\n",dims[0],dims[1],dims[2],shifts[0],shifts[1],shifts[2]);

  // creating the main array
  float *ar3D = new float[num3D];
  int *ar3Dn = NULL;
  for (size_t i=0; i<num3D; i++)
    ar3D[i] = 0.;

  // For the interpolated array
  TPoints** allPo = NULL;  
  if (interpRad > 0.1)
  { printf("Merging with distance weighted interpolation. Max. radius is %0.2f\n", interpRad);
    allPo = new TPoints*[num3D];
    for (size_t i=0; i<num3D; i++) 
      allPo[i] = NULL;
  } else // Without interp the Num in each pixel is needed
  { ar3Dn = new int[num3D];
    for (size_t i=0; i<num3D; i++)
      ar3Dn[i] = 0;
  }

  // Building an array with radii for each pixel (here only for the good pixels)
  if (radialSNR > MinVal)
  { pix_r = new int[numcomp];
    float* pixelsR = new float[numcomp];
    BuildRadialArray(numcomp, arx, ary, 1, pix_r, &maxRad, pixelsR);
    if (pixelsR != NULL) delete[] pixelsR;
  } 

  // Building a scaling array for each pixel
  float* multFact = NULL;
  if ((polariz >= 0 && polariz <= 1) || correctsolid)
  { multFact = new float[numcomp];
    for (int i=0; i<numcomp; i++)
      multFact[i] = 1.;
    if (polariz >= 0 && polariz <= 1)
      CorrectPolarisation(multFact, numcomp, arx, ary, arz, badVal, polariz); ///!!!! not arz[0]!
    if (correctsolid)
      CorrectSolidAngle(multFact, numcomp, arx, ary, arz, badVal); ///!!!! not arz[0]!
  }    

  // While using interpolation I can't use OMP, otherwise I can try to create the same pointer twice - seg.fault :(
  int numThreads = 1;
#ifdef _OPENMP
  numThreads = omp_get_max_threads();
  printf("OMP is on. Found %d threads\n",numThreads);
  if (interpRad > 0.1)
  { //omp_set_num_threads(1);
    numThreads = 1;
    printf("But interpolation is on, so for processing files using only 1 thread\n");
  }
#endif
  
  printf("Starting the merge\n");
  // Main part. Here the actual files will be loaded and added into 3D.
#pragma omp parallel for shared(numcomp, rotangles, rot, eulerangles, arx, ary, arz, dims, shifts, kwp, scale, ar3D, ar3Dn, allPo) num_threads(numThreads)
  for (int fi=0; fi<numFi; fi+=1)
  { float* apattern = NULL;
    char fName1[MaxStrLen];
    int dimsL[3] = {0,0,1};
    float angle;

    strncpy(fName1,&filesList[fi*MaxStrLen],MaxStrLen);

    // Reading a file
    int _numcomp = read_a_file(fName1, &apattern, dimsL, numfs, numss, transpose, datafield, frameaxis, frameNums[fi], &angle);
    // Checking file dimentions
    if (_numcomp != numcompFull)
    { if (_numcomp==0 || dimsL[0]<2 || dimsL[1]<1)
        printf("The file %s couldn't be read!\n", fName1);
      else   
        printf("Number of pixels in the file %s doesn't correspond the geometry %d*%d != %d \n", fName1, dimsL[0], dimsL[1], numcompFull);
      printf("Skipping file...\n");
      if (apattern != NULL)
        delete[] apattern;
      continue;
    }

    // If there was no angular info in the list - take the angle from the file
    if (!anglesPresent)
    { rotangles[fi][0] = angle*M_PI/180.;
      printf("File %s is at the angle (from the file) %0.3f\n",fName1, angle);
    }

    // Kill pixels that are lessthan, equalto and morethan
    if (lessthan < 1e-36 || morethan > -1e-36 || equalto > -1e-36)
      mask_some_values(apattern, numcompFull, lessthan, morethan, equalto, badVal);

    // Here apply corrections like local background, that need to be applyed to the FULL image. What about mask???
    // Subtract dark and multiply by gain
    if (adark != NULL)
      for (int i=0; i<numcomp; i++)
        apattern[i] -= adark[i];
    if (again != NULL)
      for (int i=0; i<numcomp; i++)
        apattern[i] *= again[i];

    //??? Maybe apply the mask before here to the full image? Than there will be difference betweem mask and maskroi
    if (localBG != 0)
      apply_lb(numfs, numss, roi, localBG, -1e10, apattern);

    // Apply radial subtraction to full pattern  
    if (radialSNR < -MinVal)
      SubtractBgLoop(apattern, numcompFull, pix_r, -radialSNR, ADCthresh, badVal, NULL);

    // now getting only values that are good
    int fulldims = 1;
    for (int i=0; i<3; i++)
      if (dimsL[i]>0) fulldims *= dimsL[i];
    if (fulldims > numcomp) // need truncation using gp array
    { float* pat1 = new float[numcomp];
      for (int i=0; i<numcomp; i++)
        pat1[i] = apattern[gp[i]];
      delete[] apattern;
      apattern = pat1;  
    }

    // Here apply corrections 
    // Polarization and solid angle corrections
    if ((polariz >= 0 && polariz <= 1) || correctsolid)
      for (int i=0; i<numcomp; i++)
        apattern[i] *= multFact[i];

    // rescale patterns, if set in the input file list
    if (scalePatterns)
      for (int i=0; i<numcomp; i++)
        apattern[i] /= scalesPat[fi];

    // Radial background correction (for good pixels only)
    if (radialSNR > MinVal)
      SubtractBgLoop(apattern, numcomp, pix_r, radialSNR, ADCthresh, badVal, NULL);

    // Mask ice rings
    if (ringDiff > MinVal)
      MaskRingsSimple(apattern, NULL, pix_r, numcomp, badVal, ringDiff, ringSmooth);

    // MaskOutliersAtRings ????   -  maybe make a wrapper around SubtractBgLoop to get the mask???

    // Saving only the 1st pattern to see what is being merged
    if (fi == 0)
      save_first_pat(fName1, dimsL, numcomp, gp, apattern);

    // now adding to 3D
    add_pattern_to_3D(apattern, numcomp, rotangles[fi], rot, eulerangles, arx, ary, arz, dims, shifts, kwp, badVal, ar3D, ar3Dn, globTrMat, allPo);
    if (debug > 0)
      if (fi > 0 && fi%10 == 0)
        printf("Processed %d files\n",fi);

    // delete the pattern
    if (apattern != NULL) delete[] apattern;
  }

  printf("Adding patterns finished. Now processing 3D\n");    
 
  // Merge together interpolated 3D volume
  if (interpRad > 0.1)
  { 
    merge_interpolated(interpRad, dims, ar3D, allPo, badValOut);
    printf("Cleaning interp\n");  
    for (size_t i=0; i<num3D; i++)
      if (allPo[i] != NULL)
        delete[] allPo[i];
    delete[] allPo;
  } else    //normalize non-interpolater 3D
    for (size_t i=0; i<num3D; i++)
    { if (ar3Dn[i] > 0)
        ar3D[i] /= (float)ar3Dn[i];
      else 
        ar3D[i] = badValOut;
    }

  strcpy(maskName, fName);
  strcpy(maskName, ExtractJustFileNameNoExtC(maskName));
 
  if (saveNums)
    if (interpRad < 0.1) 
    { float* tmp3D = new float[num3D];
      for (size_t i=0; i<num3D; i++)
        tmp3D[i] = (float)ar3Dn[i];  
      char str1[MaxStrLen];
      sprintf(str1,"%s_nums", maskName);
//      save_3D(convertTo, str1, tmp3D, num3D, dims);
      save_3D_and_Slices(convertTo, str1, tmp3D, num3D, dims);
      delete[] tmp3D;
      printf("Saved nums as %s\n", str1);
    } else
      printf("Cannot save number of pixels per voxel in interpolation mode\n");

  // Clearning
  if (pix_r != NULL) delete[] pix_r;
  if (ar3Dn != NULL) delete[] ar3Dn;
  if (multFact != NULL) delete[] multFact;

  // Postprocessing 3D   ---  REMOVE to postprocess!
  // Radial 3D
  if (radialSNR3D > 0)
  { int* pix_r3 = NULL;
    int maxRad = 0;
    size_t dimsL[3] = {dims[0], dims[1], dims[2]};
    BuildSimpleRadialArray3D(dimsL, &pix_r3);
    SubtractBgLoop(ar3D, num3D, pix_r3, radialSNR3D, ADCthresh, badValOut, NULL);
    if (pix_r3 != NULL) delete[] pix_r3;
  }  

  // Local 3D
  if (localBG3D != 0)
    if (localBG > 0)
      SubLocalBG3D(ar3D, (size_t)dims[0], (size_t)dims[1], (size_t)dims[2], localBG3D, badValOut, NULL);
    else
    { float* smAr = new float[num3D];
      SubLocalBG3D(ar3D, (size_t)dims[0], (size_t)dims[1], (size_t)dims[2], localBG3D, badValOut, smAr);
      delete[] ar3D;
      ar3D = smAr;
    }  

  // Mask outliers
//  if (outlierSNR3D > 0)
//    MaskOutliersAtRings3D(ar3D, dims, NULL, badValOut, outlierSNR3D);

  //saving 3D
//  save_3D(convertTo, maskName, ar3D, num3D, dims);
  save_3D_and_Slices(convertTo, maskName, ar3D, num3D, dims);
/*  if (convertTo == 0)
  { strcat(maskName,".bin");
    SaveFloat3DtoFile(maskName, ar3D, num3D);
  } else if (convertTo == 1)
  { strcat(maskName,".tif");
    TIFWriter3D(maskName, ar3D, NULL, dims[0], dims[1], dims[2], "SimpleMerge by O.Yefanov");
  } else if (convertTo == 2)
  { strcat(maskName,".mrc");
    MRCWriterF(maskName, ar3D, dims[0], dims[1], dims[2], NULL, NULL, "SimpleMerge by O.Yefanov");
  }
*/
//  printf("Cleaning\n");
  if (ar3D != NULL) delete[] ar3D;
  
  return 1;
}
//---------------------------------------------------------------------------
