#include "stdio.h"
#include "time.h"
#include "stdlib.h"
#include "string.h"
#include "ctype.h"
#include "math.h"
#include <stdint.h>
//#include <sys/types.h>

//#define Hdf5
#ifdef Hdf5
#include "hdf5.h"
//#include "hdf5_hl.h"
#endif


#define dtype float
const unsigned int MaxStrLen = 1020;
const unsigned int MaxHeader = 10000;
const float MinVal = 1e-10;


int roundCpp(double xd)
{ //double _xd = xd + 103079215104.5;
  //return ((int*)&_xd)[0] >> 16;
  double t = ((xd) + 6755399441055744.0);
  return *((int *)(&t));
}

int roundMy(float x)
{ double x1 = x;
  return roundCpp(x1);
  //  return (int)round(x);
};


bool TextWriter(const char* fnam, float* data, int numEl)
{
  FILE* txtF = fopen(fnam,"wt");
  if (txtF==NULL) return false;
  for (int i=0; i<numEl; i++)
    fprintf(txtF,"%0.3f\n",data[i]);
  fclose(txtF);
  return true;
}

bool ReadCBFfile(const char* fname, dtype** outArray, int* dims, float* pixel,
                 float* expo, float* waveLen, float* dist, float* beamxy,
                 float* flux, unsigned short** badmask, float* angle,
                 float* angleInc, float badVal)                   //!!! change to arrays [2] pixel, dims, beamxy, ???
{
  const unsigned int MaxHeader = 10000;
  FILE* inpCBF = fopen(fname,"rb");
  if (inpCBF == NULL) return false;
  char* headerAr = new char[MaxHeader];
  signed char aByte;
  size_t cInd = 0;
  // Searching for the data and reading the header
  bool notCBF = false;
  while (!feof(inpCBF))
  {
    fread(&aByte,1,1,inpCBF);
    if (aByte==0xC)
    { unsigned char ch3[3];
      fread(&ch3,3,1,inpCBF);
      if (ch3[0]==0x1A && ch3[1]==0x04 && ch3[2]==0xD5) break;
    }
    //here header just copied to a big array of chars here?
    if (cInd<MaxHeader)
    { headerAr[cInd] = aByte;
      cInd++;
    } else
    { printf("Big Header! Is it realy CBF?\n");
      notCBF = true;
      break;
    }

  }
  if (feof(inpCBF) || notCBF)
  { printf("Data not found in %s\n",fname);
    return false;
  }
  // Here analyse the header
  // (can be combined with the previous one)
  char str1[MaxStrLen];
  int headerLen = cInd;
  cInd = 0;
  while (cInd<headerLen)
  { int sPo = 0;
    while (headerAr[cInd]!='\n' && cInd<headerLen && sPo<MaxStrLen-1)
    { str1[sPo] = headerAr[cInd];
      cInd++;
      sPo++;
    }
    str1[sPo]=0;
    cInd++;
    sscanf(str1,"# Pixel_size %f m x %f m", pixel,pixel+1);
    sscanf(str1,"# Exposure_time %f", expo);
    sscanf(str1,"# Wavelength %f", waveLen);
    sscanf(str1,"# Detector_distance %f", dist);
    sscanf(str1,"# Beam_xy (%f, %f) pixels", beamxy, beamxy+1);
    sscanf(str1,"# Flux %f", flux);
    sscanf(str1,"# Start_angle %f", angle);
    sscanf(str1,"# Angle_increment %f", angleInc);
    sscanf(str1,"X-Binary-Size-Fastest-Dimension: %d", dims);
    sscanf(str1,"X-Binary-Size-Second-Dimension: %d", dims+1);
  }
  size_t totLen = dims[0]*dims[1];
  if (totLen<1)
  { printf("Some Dimentions are 0!\n");
    return false;
  }
  delete headerAr;

  *outArray = new dtype[totLen];
  dtype* outAr = *outArray;

//mask  *badmask = new unsigned short[totLen];
//mask  unsigned short* badMa = *badmask;

  dtype cVal = 0;                       // what about float? Make more universal?
  // The main reading of the CBS file
  signed short aWord;
  int anInt = sizeof(char);
  cInd = 0;
  while (!feof(inpCBF) && cInd<totLen)
  { fread(&aByte,1,1,inpCBF);
    if (aByte == -128)
    { fread(&aWord,2,1,inpCBF);
      if (aWord == -32768)
      { fread(&anInt,4,1,inpCBF);
        cVal += anInt;
      } else cVal += aWord;
    } else cVal += aByte;
    outAr[cInd] = (dtype)cVal;
    if (cVal<0) outAr[cInd] = badVal;
//doesn't work :(    if (cVal<0) outAr[cInd] = PanelGapValue;
//mask    if (cVal<Threshold-1e-10)
//mask      badMa[cInd] = 0;  //mask negative intensities
//mask    else
//mask      badMa[cInd] = 1;
    cInd++;
  }
  fclose(inpCBF);

  return true;
}

bool ReadCBFheader(const char* fname, int* dims, float* pixel, float* expo,
                   float* waveLen, float* dist, float* beamxy, float* flux,
                   float* angle, float* angleInc)
{
  const unsigned int MaxHeader = 10000;
  FILE* inpCBF = fopen(fname,"rb");
  if (inpCBF == NULL) return false;
  char* headerAr = new char[MaxHeader];
  signed char aByte;
  size_t cInd = 0;
  // Searching for the data and reading the header
  bool notCBF = false;
  while (!feof(inpCBF))
  {
    fread(&aByte,1,1,inpCBF);
    if (aByte==0xC)
    { unsigned char ch3[3];
      fread(&ch3,3,1,inpCBF);
      if (ch3[0]==0x1A && ch3[1]==0x04 && ch3[2]==0xD5) break;
    }
    //here header just copied to a big array of chars here?
    if (cInd<MaxHeader)
    { headerAr[cInd] = aByte;
      cInd++;
    } else
    { printf("Big Header! Is it realy CBF?\n");
      notCBF = true;
      break;
    }

  }
  if (feof(inpCBF) || notCBF)
  { printf("Data not found in %s\n",fname);
    return false;
  }
  // Here analyse the header
  // (can be combined with the previous one)
  char str1[MaxStrLen];
  int headerLen = cInd;
  cInd = 0;
  while (cInd<headerLen)
  { int sPo = 0;
    while (headerAr[cInd]!='\n' && cInd<headerLen && sPo<MaxStrLen-1)
    { str1[sPo] = headerAr[cInd];
      cInd++;
      sPo++;
    }
    str1[sPo]=0;
    cInd++;
    sscanf(str1,"# Pixel_size %f m x %f m", pixel,pixel+1);
    sscanf(str1,"# Exposure_time %f", expo);
    sscanf(str1,"# Wavelength %f", waveLen);
    sscanf(str1,"# Detector_distance %f", dist);
    sscanf(str1,"# Beam_xy (%f, %f) pixels", beamxy, beamxy+1);
    sscanf(str1,"# Flux %f", flux);
    sscanf(str1,"# Start_angle %f", angle);
    sscanf(str1,"# Angle_increment %f", angleInc);
    sscanf(str1,"X-Binary-Size-Fastest-Dimension: %d", dims);
    sscanf(str1,"X-Binary-Size-Second-Dimension: %d", dims+1);
  }
  size_t totLen = dims[0]*dims[1];
  if (totLen<1)
  { printf("Some Dimentions are 0!\n");
    return false;
  }
  delete headerAr;

  fclose(inpCBF);

  return true;
}

bool WriteCBFfile(const char* fname, dtype* outArray, int* dims, float* waveLen,
                  float* dist, float* pixelSize, float* expo, float* beamxy,
                  float* flux)
{
  const unsigned char DataSep[4] = {0x0C,0x1A,0x04,0xD5};

  if (dims[0]<1)
  { printf("(CBF) Something wrong, dimention[0]<1, not writing\n");
    return false;
  }
  FILE* outCBF = fopen(fname,"wb");
  if (outCBF == NULL)
  { printf("(CBF) File %s couldn't be created.\n",fname);
    return false;
  }
  size_t dataSize = (size_t)dims[0]*(dims[1]>0?(size_t)dims[1]:(size_t)1);

  //the length of the compressed array (from CASS)
  int nBytes = 0;
  int pixvalue = 0;
  int diff;
  int absdiff;
  for (int iadr=0; iadr<dataSize; ++iadr)
  {
    diff = roundMy(outArray[iadr]) - pixvalue;
    pixvalue = roundMy(outArray[iadr]);

    absdiff = abs(diff);
    ++nBytes;
    if (absdiff < 128)
      continue;
    nBytes += 2;
    if (absdiff > 32768)
      nBytes += 4;
  }

  char str1[MaxStrLen];
  fprintf(outCBF,"###CBF: VERSION 1.5, simple writer made by O.Yefanov\r\n\r\n");
  strcpy(str1,fname);
  strrchr(str1,'.')[0] = 0;
  fprintf(outCBF,"%s\r\n\r\n",str1);
  fprintf(outCBF,"_array_data.header_convention \"PILATUS_1.2\"\r\n");
  fprintf(outCBF,"_array_data.header_contents\r\n");
  fprintf(outCBF,";\r\n");
  // here put comments
  fprintf(outCBF,"# Detector: PILATUS 6MF-0109\r\n");
  fprintf(outCBF,"# 2019-03-27T01:13:49.952\r\n");
  fprintf(outCBF,"# Pixel_size %f m x %f m\r\n", pixelSize,pixelSize+1);
  fprintf(outCBF,"# Exposure_time %f\r\n", expo);
  fprintf(outCBF,"# Exposure_period %f\r\n", expo);
  fprintf(outCBF,"# Wavelength %f\r\n", waveLen);
  fprintf(outCBF,"# Detector_distance %f\r\n", dist);
  fprintf(outCBF,"# Beam_xy (%f, %f) pixels\r\n", beamxy, beamxy+1);
  fprintf(outCBF,"# Flux %f\r\n", flux);
//  for (int i=0; i<numAddPar; i++)
//    fprintf(outCBF,"# %s\r\n",strPar[i]);
//  fprintf(outCBF,"\r\n");

  fprintf(outCBF,";\r\n\r\n");
  fprintf(outCBF,"_array_data.data\r\n;\r\n");
  fprintf(outCBF,"--CIF-BINARY-FORMAT-SECTION--\r\n");
  fprintf(outCBF,"Content-Type: application/octet-stream;\r\n");
  fprintf(outCBF,"     conversions=\"x-CBF_BYTE_OFFSET\"\r\n");
  fprintf(outCBF,"Content-Transfer-Encoding: BINARY\r\n");
  fprintf(outCBF,"X-Binary-Size: %d\r\n", nBytes);
  fprintf(outCBF,"X-Binary-ID: 1\r\n");
  fprintf(outCBF,"X-Binary-Element-Type: \"signed 32-bit integer\"\r\n");
  fprintf(outCBF,"X-Binary-Element-Byte-Order: LITTLE_ENDIAN\r\n");
  fprintf(outCBF,"Content-MD5:\r\n");
  fprintf(outCBF,"X-Binary-Number-of-Elements: %ld\r\n", dataSize);
  fprintf(outCBF,"X-Binary-Size-Fastest-Dimension: %d\r\n",dims[0]);
  fprintf(outCBF,"X-Binary-Size-Second-Dimension: %d\r\n",dims[1]);
  fprintf(outCBF,"X-Binary-Size-Padding: 4095\r\n\r\n");
  fwrite(DataSep,1,4,outCBF);

/*
  // determine endianness
  int step, first2, last2, first4, last4;
  union
  {
    uint32_t ii;
    char cc[4];
  } bint = {0x01020304};
  if (bint.cc[0] == 1)
    step = -1; // big endian
  else
    step =  1; // little endian

  // I can write always little end.
*/

  //write compressed data
  signed char aByte;
  signed short aWord;
  signed int anInt;
  int cVal = 0;
  for (size_t ip=0; ip<dataSize; ip++)
  { int dif = roundMy(outArray[ip]) - cVal;
    int absdif = abs(dif);
    cVal = roundMy(outArray[ip]);
    aByte = -128;
    if (absdif<128)
    { aByte = dif;
      fwrite(&aByte,sizeof(aByte),1,outCBF);
      continue;
    } else
      fwrite(&aByte,sizeof(aByte),1,outCBF);
    aWord = -32768;
    if (absdif<32768)
    { aWord = dif;
      fwrite(&aWord,sizeof(aWord),1,outCBF);
      continue;
    } else
      fwrite(&aWord,sizeof(aWord),1,outCBF);
    anInt = dif;
    fwrite(&anInt,sizeof(anInt),1,outCBF);
  }

  //padding with 0 and writing trailer
  char zerO = 0;
  for (int i=0; i<4096; i++)
    fwrite(&zerO,1,1,outCBF);
  fprintf(outCBF,"\r\n--CIF-BINARY-FORMAT-SECTION----\r\n;\r\n\r\n");
  fclose(outCBF);

  return true;
}

bool ReadIMGfile(const char* fname, dtype** outArray, int* dims, float* pixel,
                 float* expo, float* waveLen, float* dist, float* beamxy,
                 float* flux, unsigned short** badmask, float* angle,
                 float* angleInc)                   //!!! change to arrays [2] pixel, dims, beamxy, ???
{
  FILE* inpIMG = fopen(fname,"rb");
//  if (inpIMG == NULL) return false;
  const unsigned int maxHeader = 100;
  int numlin = 0;
  char aline[MaxStrLen];
  int headerbytes = 1024;
  while (!feof(inpIMG) && numlin<maxHeader && strchr(aline,'}')==NULL)
  { fgets(aline, MaxStrLen, inpIMG);
    sscanf(aline,"HEADER_BYTES=%d;", &headerbytes);
//    sscanf(aline,"# Flux %f", flux);
    sscanf(aline,"PIXEL_SIZE=%f;", pixel);
    sscanf(aline,"TIME=%f;", expo);
    sscanf(aline,"DISTANCE=%f;", dist);
    sscanf(aline,"WAVELENGTH=%f;", waveLen);
    sscanf(aline,"PHI=%f;", angle);
    sscanf(aline,"OSC_RANGE=%f;", angleInc);
    sscanf(aline,"BEAM_CENTER_X=%f;", beamxy);
    sscanf(aline,"BEAM_CENTER_Y=%f;", beamxy+1);
    sscanf(aline,"SIZE1=%d;", dims);
    sscanf(aline,"SIZE2=%d;", dims+1);
    numlin++;
  }
  if (feof(inpIMG) || numlin>=maxHeader)
  { printf("The file %s doesn't look like .img\n",fname);
    return false;
  }

  pixel[0] *= 0.001;
  pixel[1] = pixel[0];
  *dist *= 0.001;

  size_t totLen = dims[0]*dims[1];
  if (totLen<1)
  { printf("Some dimentions are 0!\n");
    return false;
  }

  *outArray = new dtype[totLen];
  dtype* outAr = *outArray;

  rewind(inpIMG);
  fseek(inpIMG, headerbytes, SEEK_SET);
  for (size_t i=0; i<totLen; i++)
  { uint16_t aval = 0;
    fread(&aval, sizeof(uint16_t), 1, inpIMG);
    outAr[i] = aval;
  }
  fclose(inpIMG);

  return true;
}


bool ReadEDFfile(const char* fname, dtype** outArray, int* dims, float* pixel,
                 float* expo, float* waveLen, float* dist, float* beamxy,
                 float* flux, unsigned short** badmask, float* angle,
                 float* angleInc)                   //!!! change to arrays [2] pixel, dims, beamxy, ???
{
  FILE* inpEDF = fopen(fname,"rb");
  if (inpEDF == NULL) 
  { printf("File %s not oppened!\n",fname);
    return false;
  }
  int maxHeader = 100;
  int numlin = 0;
  char aline[MaxStrLen];
  strcpy(aline,"");
  int headerbytes = 1024;
//DataType = UnsignedShort ;
  while (!feof(inpEDF) && numlin<maxHeader && strchr(aline,'}')==NULL)
  { fgets(aline, MaxStrLen, inpEDF);
//    sscanf(aline,"# Flux %f", flux);
    sscanf(aline,"Dim_1 = %d;", dims);
    sscanf(aline,"Dim_2 = %d;", dims+1);
    sscanf(aline,"acq_expo_time = %f;", expo);
//    sscanf(aline,"PIXEL_SIZE=%f;", pixel);
//    sscanf(aline,"DISTANCE=%f;", dist);
//    sscanf(aline,"WAVELENGTH=%f;", waveLen);
//    sscanf(aline,"PHI=%f;", angle);
//    sscanf(aline,"OSC_RANGE=%f;", angleInc);
//    sscanf(aline,"BEAM_CENTER_X=%f;", beamxy);
//    sscanf(aline,"BEAM_CENTER_Y=%f;", beamxy+1);
    numlin++;
  }
  if (feof(inpEDF) || numlin>=maxHeader)
  { printf("The file %s doesn't look like an .edf\n",fname);
    fclose(inpEDF);
    return false;
  }

  pixel[0] *= 0.001;
  pixel[1] = pixel[0];
  *dist *= 0.001;

  size_t totLen = dims[0]*dims[1];
  if (totLen<1)
  { printf("Some dimentions are 0!\n");
    fclose(inpEDF);
    return false;
  }

  *outArray = new dtype[totLen];
  dtype* outAr = *outArray;

  rewind(inpEDF);
  fseek(inpEDF, 1024, SEEK_SET);
  for (size_t i=0; i<totLen; i++)
  { uint16_t aval = 0;
    fread(&aval, sizeof(uint16_t), 1, inpEDF);
    outAr[i] = aval;
  }
  fclose(inpEDF);

  return true;
}

#undef dtype

bool TIFFloatWriter(const char* fnam, float* data, int dx, int dy, const char* comment)
{ char aChar;
  short aWord;
  int anInt;
  int commentLen = strlen(comment)+1;
  FILE* tifF = fopen(fnam,"wb");
  if (tifF == NULL) return false;
  //Header
  fputc(0x49,tifF);                // I can check endians - still from CASS
  fputc(0x49,tifF);
  aWord = 42;
  fwrite(&aWord,sizeof(aWord),1,tifF);
  anInt = 8; //offset
  fwrite(&anInt,sizeof(anInt),1,tifF);

  int numFields = 9;
  aWord = numFields; // Num of fields in IFD
  fwrite(&aWord,sizeof(aWord),1,tifF);

  aWord = 254; //NewSubfileType ???
  fwrite(&aWord,sizeof(aWord),1,tifF);
  aWord = 4; // type int
  fwrite(&aWord,sizeof(aWord),1,tifF);
  anInt = 1; // num elements
  fwrite(&anInt,sizeof(anInt),1,tifF);
  anInt = 0;
  fwrite(&anInt,sizeof(anInt),1,tifF);

  aWord = 256; // Width
  fwrite(&aWord,sizeof(aWord),1,tifF);
  aWord = 4; // type int
  fwrite(&aWord,sizeof(aWord),1,tifF);
  anInt = 1; // num elements
  fwrite(&anInt,sizeof(anInt),1,tifF);
  anInt = dx; // the width
  fwrite(&anInt,sizeof(anInt),1,tifF);

  aWord = 257; // Height
  fwrite(&aWord,sizeof(aWord),1,tifF);
  aWord = 4; // type int
  fwrite(&aWord,sizeof(aWord),1,tifF);
  anInt = 1; // num elements
  fwrite(&anInt,sizeof(anInt),1,tifF);
  anInt = dy; // the height
  fwrite(&anInt,sizeof(anInt),1,tifF);

  aWord = 258; // bits per sample
  fwrite(&aWord,sizeof(aWord),1,tifF);
  aWord = 3; // type short
  fwrite(&aWord,sizeof(aWord),1,tifF);
  anInt = 1; // num elements
  fwrite(&anInt,sizeof(anInt),1,tifF);
  aWord = 32; //
  fwrite(&aWord,sizeof(aWord),1,tifF);
  aWord = 0; // 12 byte
  fwrite(&aWord,sizeof(aWord),1,tifF);

  aWord = 262; // Photometric Interpretation ???
  fwrite(&aWord,sizeof(aWord),1,tifF);
  aWord = 3; // type short
  fwrite(&aWord,sizeof(aWord),1,tifF);
  anInt = 1; // num elements
  fwrite(&anInt,sizeof(anInt),1,tifF);
  aWord = 1; // ?
  fwrite(&aWord,sizeof(aWord),1,tifF);
  aWord = 0; // 12 byte
  fwrite(&aWord,sizeof(aWord),1,tifF);

  aWord = 273; // offset to the data
  fwrite(&aWord,sizeof(aWord),1,tifF);
  aWord = 4; // type int
  fwrite(&aWord,sizeof(aWord),1,tifF);
  anInt = 1; // num elements
  fwrite(&anInt,sizeof(anInt),1,tifF);
  anInt = 14+numFields*12+commentLen; // calculate!
  fwrite(&anInt,sizeof(anInt),1,tifF);

  aWord = 279; // data size
  fwrite(&aWord,sizeof(aWord),1,tifF);
  aWord = 4; // type int
  fwrite(&aWord,sizeof(aWord),1,tifF);
  anInt = 1; // num elements
  fwrite(&anInt,sizeof(anInt),1,tifF);
  anInt = dx*dy*sizeof(float); // calculate!
  fwrite(&anInt,sizeof(anInt),1,tifF);

  aWord = 339; // Sample Format
  fwrite(&aWord,sizeof(aWord),1,tifF);
  aWord = 3; // type short
  fwrite(&aWord,sizeof(aWord),1,tifF);
  anInt = 1; // num elements
  fwrite(&anInt,sizeof(anInt),1,tifF);
  aWord = 3; // IEEE floar
  fwrite(&aWord,sizeof(aWord),1,tifF);
  aWord = 0; // 12 byte
  fwrite(&aWord,sizeof(aWord),1,tifF);

  aWord = 305; // Software
  fwrite(&aWord,sizeof(aWord),1,tifF);
  aWord = 2; // type string
  fwrite(&aWord,sizeof(aWord),1,tifF);
  anInt = commentLen; // the length
  fwrite(&anInt,sizeof(anInt),1,tifF);
  anInt = 14+numFields*12; // offset - easy to calc
  fwrite(&anInt,sizeof(anInt),1,tifF);

  anInt = 0; // next IFD - not present
  fwrite(&anInt,sizeof(anInt),1,tifF);

  //write comment
  for (int i=0; i<commentLen; i++)
    fputc(comment[i],tifF);
//  fprintf(tifF,"%s",comment);

  // Just saving the data
  size_t numel = (size_t)dx * (size_t)dy;
  for (size_t i=0; i<numel; i++)
    fwrite(&data[i],sizeof(float),1,tifF);

  fclose(tifF);
  return true;
}

bool TIFFloatReader(const char* fnam, float** data, int* dx, int* dy, char* comment)
{ char aChar;
  short aWord;
  int anInt;
  FILE* tifF = fopen(fnam,"rb");
  if (tifF==NULL) return false;
  fread(&anInt,4,1,tifF); //Header - could be analysed for little/big endians
  fread(&anInt,4,1,tifF); // offset of IFD
  fseek(tifF,anInt,SEEK_SET);

  fread(&aWord,2,1,tifF); // number of fields
  int numFields = aWord;
  if (numFields<3) return false;

  struct tFields
  { short tag;
    short type;
    int numEl;
    short valS1;
    short valS2;
    int valI;
    char valC[4];
    float valF;
  };

  tFields*  aField = new tFields[numFields];
  for (int i=0; i<numFields; i++)  //Reading fields in IFD
  {
    fread(&aField[i].tag,2,1,tifF);
    fread(&aField[i].type,2,1,tifF);
    fread(&aField[i].numEl,4,1,tifF);
    if (aField[i].type == 3)
    { fread(&aField[i].valS1,2,1,tifF);
      fread(&aField[i].valS2,2,1,tifF);
    } else if (aField[i].type == 4 || aField[i].type == 2)
      fread(&aField[i].valI,4,1,tifF);
      else if (aField[i].type == 1)
      fread(aField[i].valC,1,4,tifF);
      else if (aField[i].type == 11)
      fread(&aField[i].valF,4,1,tifF);
  }

  int imOffset = 0;
  int dataSize = 0;
  float dataType = 0;
  bool rightType = true;                                        //!!!Corrected 131023
  for (int i=0; i<numFields; i++)  //Analysing fields in IFD
  { if (aField[i].tag == 256) *dx = aField[i].valI;
    else if (aField[i].tag == 257) *dy = aField[i].valI;
    else if (aField[i].tag == 273) imOffset = aField[i].valI;
    else if (aField[i].tag == 258)
      rightType = (aField[i].valS1==32);
    else if (aField[i].tag == 339)
      dataType = aField[i].valS1;
    else if (aField[i].tag == 279)
      dataSize = aField[i].valI;
    else if (aField[i].tag == 305) // This is string reader
    { fseek(tifF,aField[i].valI,SEEK_SET);
      if (aField[i].numEl>MaxStrLen) aField[i].numEl=MaxStrLen;
      for (size_t k=0; k<aField[i].numEl; k++)
        fread(&comment[k],1,1,tifF);
//      aField[i].valI;
    }
  }
  if (!rightType || (dataType!=2 && dataType!=3))
  { printf("Supported data types are: float or int (both 32bit)!\n");
    return false;
  }

  // creating array for the image
  if (imOffset<8 || *dx<1 || *dy<1) return false; //???
  size_t numEl = *dx * *dy;
  if (numEl<1) return false;
#ifdef linux
  *data = (float*)valloc(numEl*sizeof(float));
#else   //090615
  *data = new float[numEl];
#endif
  float* _Ar = *data;

  // Reading the actual image
  fseek(tifF,imOffset,SEEK_SET);
  if (dataType==3)
    fread(_Ar,sizeof(float),numEl,tifF);
  else if (dataType==2)
  { int32_t aval;
    for (int i=0; i<numEl; i++)
    { fread(&aval,4,1,tifF);
      _Ar[i] = (float)aval;
    }
  }

  delete aField;
  fclose(tifF);
  return true;
}

bool TIFWriter(const char* fnam, float* dataF, int* dataI, int dx, int dy, const char* comment, char tifdtype)  // tifdtype=2 - int, tifdtype=3 - float
{ char aChar;
  short aWord;
  int anInt;
  int commentLen = strlen(comment)+1;
  FILE* tifF = fopen(fnam,"wb");
  if (tifF == NULL) return false;
  //Header
  fputc(0x49,tifF);                // I can check endians - still from CASS
  fputc(0x49,tifF);
  aWord = 42;
  fwrite(&aWord,sizeof(aWord),1,tifF);
  anInt = 8; //offset
  fwrite(&anInt,sizeof(anInt),1,tifF);

  int numFields = 9;
  aWord = numFields; // Num of fields in IFD
  fwrite(&aWord,sizeof(aWord),1,tifF);

  aWord = 254; //NewSubfileType ???
  fwrite(&aWord,sizeof(aWord),1,tifF);
  aWord = 4; // type int
  fwrite(&aWord,sizeof(aWord),1,tifF);
  anInt = 1; // num elements
  fwrite(&anInt,sizeof(anInt),1,tifF);
  anInt = 0;
  fwrite(&anInt,sizeof(anInt),1,tifF);

  aWord = 256; // Width
  fwrite(&aWord,sizeof(aWord),1,tifF);
  aWord = 4; // type int
  fwrite(&aWord,sizeof(aWord),1,tifF);
  anInt = 1; // num elements
  fwrite(&anInt,sizeof(anInt),1,tifF);
  anInt = dx; // the width
  fwrite(&anInt,sizeof(anInt),1,tifF);

  aWord = 257; // Height
  fwrite(&aWord,sizeof(aWord),1,tifF);
  aWord = 4; // type int
  fwrite(&aWord,sizeof(aWord),1,tifF);
  anInt = 1; // num elements
  fwrite(&anInt,sizeof(anInt),1,tifF);
  anInt = dy; // the height
  fwrite(&anInt,sizeof(anInt),1,tifF);

  aWord = 258; // bits per sample
  fwrite(&aWord,sizeof(aWord),1,tifF);
  aWord = 3; // type short
  fwrite(&aWord,sizeof(aWord),1,tifF);
  anInt = 1; // num elements
  fwrite(&anInt,sizeof(anInt),1,tifF);
  aWord = 32; //
  fwrite(&aWord,sizeof(aWord),1,tifF);
  aWord = 0; // 12 byte
  fwrite(&aWord,sizeof(aWord),1,tifF);

  aWord = 262; // Photometric Interpretation ???
  fwrite(&aWord,sizeof(aWord),1,tifF);
  aWord = 3; // type short
  fwrite(&aWord,sizeof(aWord),1,tifF);
  anInt = 1; // num elements
  fwrite(&anInt,sizeof(anInt),1,tifF);
  aWord = 1; // ?
  fwrite(&aWord,sizeof(aWord),1,tifF);
  aWord = 0; // 12 byte
  fwrite(&aWord,sizeof(aWord),1,tifF);

  aWord = 273; // offset to the data
  fwrite(&aWord,sizeof(aWord),1,tifF);
  aWord = 4; // type int
  fwrite(&aWord,sizeof(aWord),1,tifF);
  anInt = 1; // num elements
  fwrite(&anInt,sizeof(anInt),1,tifF);
  anInt = 14+numFields*12+commentLen; // calculate!
  fwrite(&anInt,sizeof(anInt),1,tifF);

  aWord = 279; // data size
  fwrite(&aWord,sizeof(aWord),1,tifF);
  aWord = 4; // type int
  fwrite(&aWord,sizeof(aWord),1,tifF);
  anInt = 1; // num elements
  fwrite(&anInt,sizeof(anInt),1,tifF);
  anInt = dx*dy*4;//sizeof(float); // calculate!
  fwrite(&anInt,sizeof(anInt),1,tifF);

  aWord = 339; // Sample Format
  fwrite(&aWord,sizeof(aWord),1,tifF);
  aWord = 3; // type short
  fwrite(&aWord,sizeof(aWord),1,tifF);
  anInt = 1; // num elements
  fwrite(&anInt,sizeof(anInt),1,tifF);
  aWord = 3; // IEEE float
  if (tifdtype==2) aWord = 2;  // INTEGER                  // DATA TYPE
  fwrite(&aWord,sizeof(aWord),1,tifF);
  aWord = 0; // 12 byte
  fwrite(&aWord,sizeof(aWord),1,tifF);

  aWord = 305; // Software
  fwrite(&aWord,sizeof(aWord),1,tifF);
  aWord = 2; // type string
  fwrite(&aWord,sizeof(aWord),1,tifF);
  anInt = commentLen; // the length
  fwrite(&anInt,sizeof(anInt),1,tifF);
  anInt = 14+numFields*12; // offset - easy to calc
  fwrite(&anInt,sizeof(anInt),1,tifF);

  anInt = 0; // next IFD - not present
  fwrite(&anInt,sizeof(anInt),1,tifF);

  //write comment
  for (int i=0; i<commentLen; i++)
    fputc(comment[i],tifF);
//  fprintf(tifF,"%s",comment);

  // Just saving the data
  size_t numel = (size_t)dx * (size_t)dy;
  if (tifdtype==2)                                // DATA TYPE - int 32 bit
    for (size_t i=0; i<numel; i++)
    { int32_t aval = dataI[i];
      fwrite(&aval,4,1,tifF);
    } else
//    for (size_t i=0; i<numel; i++)
      fwrite(dataF,sizeof(float),numel,tifF);

  fclose(tifF);
  return true;
}

bool TIFReader(const char* fnam, float** data, int* dx, int* dy, char* comment)
{ char aChar;
  short aWord;
  int anInt;
  FILE* tifF = fopen(fnam,"rb");
  if (tifF==NULL) return false;
//  fread(&anInt,4,1,tifF); //Header - could be analysed for little/big endians
  fread(&aChar,1,1,tifF); // little/big endian?
  if (aChar != 0x49)
  { printf("Reading ONLY Intel byte order (little big endian) tifs!\n");
    return false;
  }
  fread(&aChar,1,1,tifF); // repeat
  fread(&aWord,2,1,tifF); // something

  fread(&anInt,4,1,tifF); // offset of IFD
  fseek(tifF,anInt,SEEK_SET);

  fread(&aWord,2,1,tifF); // number of fields
  int numFields = aWord;
  if (numFields<3) return false;

  struct tFields
  { short tag;
    short type;
    size_t numEl;
    short valS1;
    short valS2;
    int valI;
    char valC[4];
    float valF;
  };

  tFields*  aField = new tFields[numFields];
  for (int i=0; i<numFields; i++)  //Reading fields in IFD
  {
    fread(&aField[i].tag,2,1,tifF);
    fread(&aField[i].type,2,1,tifF);
    fread(&aField[i].numEl,4,1,tifF);
    if (aField[i].type == 3)
    { fread(&aField[i].valS1,2,1,tifF);
      fread(&aField[i].valS2,2,1,tifF);
    } else if (aField[i].type == 4 || aField[i].type == 2)
      fread(&aField[i].valI,4,1,tifF);
      else if (aField[i].type == 1)
      fread(aField[i].valC,1,4,tifF);
      else if (aField[i].type == 11)
      fread(&aField[i].valF,4,1,tifF);
  }

  int imOffset = 0;
  int dataSize = 0;
  int dataType = 0;
  int bitsType = true;
  for (int i=0; i<numFields; i++)  //Analysing fields in IFD
  { if (aField[i].tag == 256) *dx = aField[i].valI;
    else if (aField[i].tag == 257) *dy = aField[i].valI;
    else if (aField[i].tag == 273) imOffset = aField[i].valI;
    else if (aField[i].tag == 258)
      bitsType = aField[i].valS1;//==32);
    else if (aField[i].tag == 339)
      dataType = aField[i].valS1;
    else if (aField[i].tag == 279)
      dataSize = aField[i].valI;   //??? do I need it?
    else if (aField[i].tag == 305) // This is string reader
    { fseek(tifF,aField[i].valI,SEEK_SET);
      if (aField[i].numEl>MaxStrLen) aField[i].numEl=MaxStrLen;
      for (size_t k=0; k<aField[i].numEl; k++)
        fread(&comment[k],1,1,tifF);
    }
  }
  if (!(bitsType==16 || (bitsType==32 && (dataType==2 || dataType==3))))
  { printf("Supported data types are: 32bit float or 16/32 int!\n");
    return false;
  }

  // creating array for the image
  if (imOffset<8 || *dx<1 || *dy<1) return false; //???
  size_t numEl = *dx * *dy;
  if (numEl<1) return false;
  *data = new float[numEl];
  float* _Ar = *data;

  // Reading the actual image
  fseek(tifF,imOffset,SEEK_SET);
  if (dataType==3 && bitsType==32)
    fread(_Ar,sizeof(float),numEl,tifF);
  else if (dataType==2 && bitsType==32)
  { int32_t aval;
    for (size_t i=0; i<numEl; i++)
    { fread(&aval,4,1,tifF);
      _Ar[i] = (float)aval;
    }
  } else if (bitsType==16)
  { uint16_t aval;
    for (size_t i=0; i<numEl; i++)
    { fread(&aval,2,1,tifF);
      _Ar[i] = (float)aval;
    }
  }

  delete aField;
  fclose(tifF);
  return true;
}

void TIFWriteHeaderIFD(FILE* tifF, int dx, int dy, int dz, int iz, int commentLen, int tifType) //returning some address
{
  short aWord;
  int anInt;
  int numFields = 9;

  aWord = numFields; // Num of fields in IFD
  fwrite(&aWord,sizeof(aWord),1,tifF);

  aWord = 254; //NewSubfileType ???
  fwrite(&aWord,sizeof(aWord),1,tifF);
  aWord = 4; // type int
  fwrite(&aWord,sizeof(aWord),1,tifF);
  anInt = 1; // num elements
  fwrite(&anInt,sizeof(anInt),1,tifF);
  anInt = 0;
  fwrite(&anInt,sizeof(anInt),1,tifF);

  aWord = 256; // Width
  fwrite(&aWord,sizeof(aWord),1,tifF);
  aWord = 4; // type int
  fwrite(&aWord,sizeof(aWord),1,tifF);
  anInt = 1; // num elements
  fwrite(&anInt,sizeof(anInt),1,tifF);
  anInt = dx; // the width
  fwrite(&anInt,sizeof(anInt),1,tifF);

  aWord = 257; // Height
  fwrite(&aWord,sizeof(aWord),1,tifF);
  aWord = 4; // type int
  fwrite(&aWord,sizeof(aWord),1,tifF);
  anInt = 1; // num elements
  fwrite(&anInt,sizeof(anInt),1,tifF);
  anInt = dy; // the height
  fwrite(&anInt,sizeof(anInt),1,tifF);

  aWord = 258; // bits per sample
  fwrite(&aWord,sizeof(aWord),1,tifF);
  aWord = 3; // type short
  fwrite(&aWord,sizeof(aWord),1,tifF);
  anInt = 1; // num elements
  fwrite(&anInt,sizeof(anInt),1,tifF);
  aWord = 32; //
  fwrite(&aWord,sizeof(aWord),1,tifF);
  aWord = 0; // 12 byte
  fwrite(&aWord,sizeof(aWord),1,tifF);

  aWord = 262; // Photometric Interpretation ???
  fwrite(&aWord,sizeof(aWord),1,tifF);
  aWord = 3; // type short
  fwrite(&aWord,sizeof(aWord),1,tifF);
  anInt = 1; // num elements
  fwrite(&anInt,sizeof(anInt),1,tifF);
  aWord = 1; // ?
  fwrite(&aWord,sizeof(aWord),1,tifF);
  aWord = 0; // 12 byte
  fwrite(&aWord,sizeof(aWord),1,tifF);

  aWord = 273; // offset to the data                      ///!!!!!!!!!!! change
  fwrite(&aWord,sizeof(aWord),1,tifF);
  aWord = 4; // type int
  fwrite(&aWord,sizeof(aWord),1,tifF);
  anInt = 1; // num elements
  fwrite(&anInt,sizeof(anInt),1,tifF);
  anInt = 8+commentLen + (numFields*12+6)*dz + iz*dx*dy*4; // calculate!
  fwrite(&anInt,sizeof(anInt),1,tifF);

  aWord = 279; // data size
  fwrite(&aWord,sizeof(aWord),1,tifF);
  aWord = 4; // type int
  fwrite(&aWord,sizeof(aWord),1,tifF);
  anInt = 1; // num elements
  fwrite(&anInt,sizeof(anInt),1,tifF);
  anInt = dx*dy*4;//sizeof(float); // calculate!
  fwrite(&anInt,sizeof(anInt),1,tifF);

  aWord = 339; // Sample Format
  fwrite(&aWord,sizeof(aWord),1,tifF);
  aWord = 3; // type short
  fwrite(&aWord,sizeof(aWord),1,tifF);
  anInt = 1; // num elements
  fwrite(&anInt,sizeof(anInt),1,tifF);
  aWord = 3; // IEEE float
  if (tifType==2) aWord = 2;  // INTEGER                  // DATA TYPE
  fwrite(&aWord,sizeof(aWord),1,tifF);
  aWord = 0; // 12 byte
  fwrite(&aWord,sizeof(aWord),1,tifF);

  aWord = 305; // Software (comment string)        //The same for all slices
  fwrite(&aWord,sizeof(aWord),1,tifF);
  aWord = 2; // type string
  fwrite(&aWord,sizeof(aWord),1,tifF);
  anInt = commentLen; // the length
  fwrite(&anInt,sizeof(anInt),1,tifF);
  anInt = 8+(numFields*12+6)*dz; // offset - easy to calc
  fwrite(&anInt,sizeof(anInt),1,tifF);
}

bool TIFWriter3D(const char* fnam, float* dataF, int* dataI, int dx, int dy, int dz, const char* comment)  // only 1 array should be not NULL
{ char aChar;
  short aWord;
  int anInt;
  int tifType;//tifdtype=2 - int, tifdtype=3 - float
  if (dataF != NULL) tifType = 3;
  else if (dataI != NULL) tifType = 2;
  else return false;
  int commentLen = strlen(comment)+1;
  int numFields = 9;
  if (dz < 1) dz = 1;

  FILE* tifF = fopen(fnam,"wb");
  if (tifF == NULL) return false;

  //Header
  fputc(0x49,tifF);                // I can check endians - still from CASS
  fputc(0x49,tifF);
  aWord = 42;
  fwrite(&aWord,sizeof(aWord),1,tifF);
//  anInt = 8; //offset
//  fwrite(&anInt,sizeof(anInt),1,tifF);

  //multiple IFDs
  for (int iz=0; iz<dz; iz++)
  { anInt = 8 + (numFields*12+6)*iz;
    fwrite(&anInt,sizeof(anInt),1,tifF);
    TIFWriteHeaderIFD(tifF, dx, dy, dz, iz, commentLen, tifType);
  }
  anInt = 0; // next IFD - not present
  fwrite(&anInt,sizeof(anInt),1,tifF);

  //write comment
  for (int i=0; i<commentLen; i++)
    fputc(comment[i],tifF);
//  fprintf(tifF,"%s",comment);

  // Just saving the data
  size_t numel = (size_t)dx * (size_t)dy * (size_t)dz;
  if (tifType==2)                                // DATA TYPE - int 32 bit
    for (size_t i=0; i<numel; i++)
    { int32_t aval = dataI[i];
      fwrite(&aval,4,1,tifF);
    } else
    fwrite(dataF,sizeof(float),numel,tifF);

  fclose(tifF);
  return true;
}

bool TIFReader3D(const char* fnam, float** data, int* dx, int* dy, int* dz, char* comment)
{ char aChar;
  short aWord;
  int anInt;
  FILE* tifF = fopen(fnam,"rb");
  if (tifF==NULL) return false;
  fread(&anInt,4,1,tifF); //Header - could be analysed for little/big endians

  fread(&anInt,4,1,tifF); // offset of IFD
  fseek(tifF,anInt,SEEK_SET);

  fread(&aWord,2,1,tifF); // number of fields
  int numFields = aWord;
  if (numFields<3) return false;

  struct tFields
  { short tag;
    short type;
    size_t numEl;
    short valS1;
    short valS2;
    int valI;
    char valC[4];
    float valF;
  };

  tFields*  aField = new tFields[numFields];
  for (int i=0; i<numFields; i++)  //Reading fields in IFD
  {
    fread(&aField[i].tag,2,1,tifF);
    fread(&aField[i].type,2,1,tifF);
    fread(&aField[i].numEl,4,1,tifF);
    if (aField[i].type == 3)
    { fread(&aField[i].valS1,2,1,tifF);
      fread(&aField[i].valS2,2,1,tifF);
    } else if (aField[i].type == 4 || aField[i].type == 2)
      fread(&aField[i].valI,4,1,tifF);
      else if (aField[i].type == 1)
      fread(aField[i].valC,1,4,tifF);
      else if (aField[i].type == 11)
      fread(&aField[i].valF,4,1,tifF);
  }

  // Now support for 3D TIFFs
  int iz = 0;
  anInt = 1;
  while (anInt > 0)
  { fread(&anInt,4,1,tifF); // offset of IFD
    if (anInt > 0)
      fseek(tifF,anInt+(numFields*12+6)-4,SEEK_SET);
    iz++;
  }
  printf("found %d frames\n",iz);
  *dz = iz;

  int imOffset = 0;
  int dataSize = 0;
  int dataType = 0;
  int bitsType = true;
  for (int i=0; i<numFields; i++)  //Analysing fields in IFD
  { if (aField[i].tag == 256) *dx = aField[i].valI;
    else if (aField[i].tag == 257) *dy = aField[i].valI;
    else if (aField[i].tag == 273) imOffset = aField[i].valI;
    else if (aField[i].tag == 258)
      bitsType = aField[i].valS1;//==32);
    else if (aField[i].tag == 339)
      dataType = aField[i].valS1;
    else if (aField[i].tag == 279)
      dataSize = aField[i].valI;   //??? do I need it?
    else if (aField[i].tag == 305) // This is string reader
    { fseek(tifF,aField[i].valI,SEEK_SET);
      if (aField[i].numEl>MaxStrLen) aField[i].numEl=MaxStrLen;
      for (size_t k=0; k<aField[i].numEl; k++)
        fread(&comment[k],1,1,tifF);
    }
  }
  if (!(bitsType==16 || (bitsType==32 && (dataType==2 || dataType==3))))
  { printf("Supported data types are: 32bit float or 16/32 int!\n");
    return false;
  }

  // creating array for the image
  if (imOffset<8 || *dx<1 || *dy<1) return false; //???
  size_t numEl = (size_t)*dx * (size_t)*dy * (size_t)*dz;
  if (numEl<1) return false;
  *data = new float[numEl];
  float* _Ar = *data;

  // Reading the actual image
  fseek(tifF,imOffset,SEEK_SET);
  if (dataType==3 && bitsType==32)
    fread(_Ar,sizeof(float),numEl,tifF);
  else if (dataType==2 && bitsType==32)
  { int32_t aval;
    for (size_t i=0; i<numEl; i++)
    { fread(&aval,4,1,tifF);
      _Ar[i] = (float)aval;
    }
  } else if (bitsType==16)
  { uint16_t aval;
    for (size_t i=0; i<numEl; i++)
    { fread(&aval,2,1,tifF);
      _Ar[i] = (float)aval;
    }
  }

  delete aField;
  fclose(tifF);
  return true;
}

bool RAWFloatReader(const char* fnam, float** data, size_t* numEl)
{
  FILE* rawF = fopen(fnam,"rb");
  if (rawF==NULL) return false;
#ifdef linux
  *data = (float*)valloc(*numEl*sizeof(float));
#else   //090615
  *data = new float[*numEl];
#endif
  float* _Ar = *data;

  size_t i = 0;
  while (!feof(rawF))
  { fread(&_Ar[i],sizeof(float),1,rawF);
    i++;
    if (i>=*numEl) break;
  }
  *numEl = i;
  fclose(rawF);
  return true;
}

bool MRCWriterF(const char* fnam, float* dataF, int dx, int dy, int dz, float cell[6], char axis[3], const char* comment)
{ //unsigned __int32 aWord41;
  //unsigned int32_t aWord4;
  int32_t aWord4;
  float aFloat;
  char ch;

  FILE* mrcF = fopen(fnam,"wb");
  if (mrcF == NULL) return false;
  //Header
  aWord4 = dx;
  fwrite(&aWord4,sizeof(aWord4),1,mrcF);
  aWord4 = dy;
  fwrite(&aWord4,sizeof(aWord4),1,mrcF);
  aWord4 = dz;
  fwrite(&aWord4,sizeof(aWord4),1,mrcF);
  //Mode: 0 8-bit signed integer (range -128 to 127)
  //1 16-bit signed integer
  //2 32-bit signed real
  //3 transform : complex 16-bit integers
  //4 transform : complex 32-bit reals
  //6 16-bit unsigned integer
  aWord4 = 2; // for now only Float
  fwrite(&aWord4,sizeof(aWord4),1,mrcF);
  aWord4 = (int)(0.5*(dx)-0.01);//0? //NXSTART - location of first column in unit cell //0.5*(Dim[0]-1.)+0.01
  fwrite(&aWord4,sizeof(aWord4),1,mrcF);
  aWord4 = (int)(0.5*(dy)-0.01); //NYSTART
  fwrite(&aWord4,sizeof(aWord4),1,mrcF);
  aWord4 = (int)(0.5*(dz)-0.01); //NZSTART
  fwrite(&aWord4,sizeof(aWord4),1,mrcF);
  aWord4 = dx; //MX - sampling along X axis of unit cell
  fwrite(&aWord4,sizeof(aWord4),1,mrcF);
  aWord4 = dy; //MY
  fwrite(&aWord4,sizeof(aWord4),1,mrcF);
  aWord4 = dz; //MZ
  fwrite(&aWord4,sizeof(aWord4),1,mrcF);
  aFloat = cell[0]; //Cell dimentions in A
  fwrite(&aFloat,sizeof(aWord4),1,mrcF);
  aFloat = cell[1]; //
  fwrite(&aFloat,sizeof(aWord4),1,mrcF);
  aFloat = cell[2]; //
  fwrite(&aFloat,sizeof(aWord4),1,mrcF);
  aFloat = cell[3]; //Cell angles in deg
  fwrite(&aFloat,sizeof(aWord4),1,mrcF);
  aFloat = cell[4]; //
  fwrite(&aFloat,sizeof(aWord4),1,mrcF);
  aFloat = cell[5]; //
  fwrite(&aFloat,sizeof(aWord4),1,mrcF);
  aWord4 = axis[0]; //axis corresp to cols (1,2,3 for X,Y,Z)
  fwrite(&aWord4,sizeof(aWord4),1,mrcF);
  aWord4 = axis[1]; //axis corresp to rows (1,2,3 for X,Y,Z)
  fwrite(&aWord4,sizeof(aWord4),1,mrcF);
  aWord4 = axis[2]; //axis corresp to sections (1,2,3 for X,Y,Z)
  fwrite(&aWord4,sizeof(aWord4),1,mrcF);

  aWord4 = 0; //minimum density value
  fwrite(&aWord4,sizeof(aWord4),1,mrcF);
  aWord4 = 0; //maximum density value
  fwrite(&aWord4,sizeof(aWord4),1,mrcF);
  aWord4 = 0; //mean density value
  fwrite(&aWord4,sizeof(aWord4),1,mrcF);
  aWord4 = 1; //space group number
  fwrite(&aWord4,sizeof(aWord4),1,mrcF);
  aWord4 = 0; //size of extended header
  fwrite(&aWord4,sizeof(aWord4),1,mrcF);

  aWord4 = 0; //extra space used for anything
  for (int i=0; i<23; i++)
    fwrite(&aWord4,sizeof(aWord4),1,mrcF);

  aWord4 = 0; //code for the type of extended header
  fwrite(&aWord4,sizeof(aWord4),1,mrcF);
  aWord4 = 20140; //version of the MRC format
  fwrite(&aWord4,sizeof(aWord4),1,mrcF);
  aWord4 = (int)(0.5*(dx)-0.01)-dx+1; //phase origin (pixels) or origin of subvolume (A)
  fwrite(&aWord4,sizeof(aWord4),1,mrcF);
  aWord4 = (int)(0.5*(dy)-0.01)-dy+1; //
  fwrite(&aWord4,sizeof(aWord4),1,mrcF);
  aWord4 = (int)(0.5*(dz)-0.01)-dz+1; //
  fwrite(&aWord4,sizeof(aWord4),1,mrcF);
  printf("origin is %d\n",aWord4);

//  aWord4 = 0; //'MAP '
//  fwrite(&aWord4,sizeof(aWord4),1,mrcF);
/*
  ch='M';
  fwrite(&ch,sizeof(ch),1,mrcF);
  ch='A';
  fwrite(&ch,sizeof(ch),1,mrcF);
  ch='P';
  fwrite(&ch,sizeof(ch),1,mrcF);
  ch=' ';
  fwrite(&ch,sizeof(ch),1,mrcF);
*/
//  char DataInd[4] = {'M','A','P',' '};
  char DataInd[4] = {'H','K','L',' '};
  fwrite(DataInd,1,4,mrcF);


  ch=0x44; //little endian (0x11 for big)
  fwrite(&ch,sizeof(ch),1,mrcF);
  ch=0x44; //little endian (0x11 for big)
  fwrite(&ch,sizeof(ch),1,mrcF);
  ch=0;
  fwrite(&ch,sizeof(ch),1,mrcF);
  ch=0;
  fwrite(&ch,sizeof(ch),1,mrcF);

  aWord4 = 0; //rms deviation of map
  fwrite(&aWord4,sizeof(aWord4),1,mrcF);
  aWord4 = 1; //number of labels being used
  fwrite(&aWord4,sizeof(aWord4),1,mrcF);

  //write comment
  for (int i=0; i<strlen(comment); i++)
    fputc(comment[i],mrcF);

  rewind(mrcF);
  fseek(mrcF, 1024, SEEK_SET);

  // Just saving the data
  size_t numel = (size_t)dx * (size_t)dy * (size_t)dz;
    fwrite(dataF,sizeof(float),numel,mrcF);

  fclose(mrcF);
  return true;


}

bool MRCReader(const char* fnam, float** dataF, int* dx, int* dy, int* dz, float cell[6], char* comment)
{
  uint32_t aWord4, dType;
  float aFloat;
//  char ch;
  char axis[3];

  FILE* mrcF = fopen(fnam,"rb");
  if (mrcF == NULL) return false;
  //Header
  fread(&aWord4,4,1,mrcF);
  *dx = aWord4;
  fread(&aWord4,4,1,mrcF);
  *dy = aWord4;
  fread(&aWord4,4,1,mrcF);
  *dz = aWord4;

  //Mode: 0 8-bit signed integer (range -128 to 127)
  //1 16-bit signed integer
  //2 32-bit signed real
  //3 transform : complex 16-bit integers
  //4 transform : complex 32-bit reals
  //6 16-bit unsigned integer
  fread(&aWord4,4,1,mrcF);
  dType = aWord4;

  fread(&aWord4,4,1,mrcF); //OrX
  fread(&aWord4,4,1,mrcF); //OrY
  fread(&aWord4,4,1,mrcF); //OrZ
  fread(&aWord4,4,1,mrcF); //MX
  fread(&aWord4,4,1,mrcF); //MY
  fread(&aWord4,4,1,mrcF); //MZ

  fread(&aFloat,4,1,mrcF); //cell
  cell[0] = aFloat;
  fread(&aFloat,4,1,mrcF); //cell
  cell[1] = aFloat;
  fread(&aFloat,4,1,mrcF); //cell
  cell[2] = aFloat;
  fread(&aFloat,4,1,mrcF); //cell
  cell[3] = aFloat;
  fread(&aFloat,4,1,mrcF); //cell
  cell[4] = aFloat;
  fread(&aFloat,4,1,mrcF); //cell
  cell[5] = aFloat;

  fread(&aWord4,4,1,mrcF); 
  axis[0] = aWord4;
  fread(&aWord4,4,1,mrcF); 
  axis[1] = aWord4;
  fread(&aWord4,4,1,mrcF); 
  axis[2] = aWord4;

  fread(&aWord4,4,1,mrcF); //density
  fread(&aWord4,4,1,mrcF); 
  fread(&aWord4,4,1,mrcF);

  fread(&aWord4,4,1,mrcF); //space group number
  fread(&aWord4,4,1,mrcF);  //size of extended header

  for (int i=0; i<23; i++)
    fread(&aWord4,4,1,mrcF); //extra space 

  fread(&aWord4,4,1,mrcF); //code for the type of extended header
  fread(&aWord4,4,1,mrcF); //20140; //version of the MRC format
  fread(&aWord4,4,1,mrcF); //phase origin (pixels) or origin of subvolume (A)
  fread(&aWord4,4,1,mrcF); 
  fread(&aWord4,4,1,mrcF); 

  char DataInd[4] = {'H','K','L',' '};
  fread(&DataInd[0],4,1,mrcF); 

  fread(&aWord4,4,1,mrcF);  //Endian (little/big)

  fread(&aWord4,4,1,mrcF); //rms deviation of map
  fread(&aWord4,4,1,mrcF); //number of labels being used

  if (aWord4 > 0)
//    for (int i=0; i<aWord4; i++)
  //read comment
  for (int i=0; i<80; i++)
    comment[i] = fgetc(mrcF);

  rewind(mrcF);

  // Reading data
  fseek(mrcF, 1024, SEEK_SET);

  size_t numEl = (*dx>0?(size_t)*dx:1)*(*dy>0?(size_t)*dy:1)*(*dz>0?(size_t)*dz:1);
  *dataF = new float[numEl];
  float* _Ar = *dataF;


  if (dType==2)
    fread(_Ar,sizeof(float),numEl,mrcF);
/*
  else if (dataType==2 && bitsType==32)
  { int32_t aval;
    for (size_t i=0; i<numEl; i++)
    { fread(&aval,4,1,tifF);
      _Ar[i] = (float)aval;
    }
  } else if (bitsType==16)
  { uint16_t aval;
    for (size_t i=0; i<numEl; i++)
    { fread(&aval,2,1,tifF);
      _Ar[i] = (float)aval;
    }
  }
*/
  return true;
}

//------------------------------------------------------------------

bool SaveFloat3DtoFile(const char* _fname, float* inpAr, long numEl)
{ FILE *stream;
  stream = fopen(_fname, "wb");
  if (!stream) return false;
  for (long i=0; i<numEl; i++)
	  fwrite(&inpAr[i],sizeof(float),1,stream);
  fclose(stream);
  return true;
}

bool ReadFloat3DfromFileU(const char* _fname, float* inpAr, size_t* numEl)
{ FILE *stream;
  stream = fopen(_fname, "rb");
  if (!stream)
  { printf("file %s not found\n",_fname);
    return false;
  }
  long i = 0;
  while (!feof(stream))
  { fread(&inpAr[i],sizeof(float),1,stream);
    i++;
    if (i>=*numEl) break;
  }
  *numEl = i;
  fclose(stream);
  return true;

}

size_t ReadFloatFromFile(const char* _fname, float* inpAr, size_t numEl)
{ FILE *stream;
  stream = fopen(_fname, "rb");
  if (!stream)
  { printf("file %s not found\n",_fname);
    return 0;
  }
  size_t i = 0;
  while (!feof(stream))
  { fread(&inpAr[i],sizeof(float),1,stream);
    i++;
    if (i>=numEl) break;
  }
  fclose(stream);
  return i;
}

bool SaveInt1DtoTextFile(const char* _fname, int* inpAr, long numEl)
{ FILE *stream;
  stream = fopen(_fname, "wt");
  if (!stream) return false;
  for (long i=0; i<numEl; i++)
	fprintf(stream,"%d\n",inpAr[i]);
  fclose(stream);
  return true;
}

bool SaveFloat1DtoTextFile(const char* _fname, float* inpAr, long numEl)
{ FILE *stream;
  stream = fopen(_fname, "wt");
  if (!stream) return false;
  for (long i=0; i<numEl; i++)
	fprintf(stream,"%0.3f\n",inpAr[i]);
  fclose(stream);
  return true;
}

bool Save2Float1DtoTextFile(const char* _fname, float* inpAr1, float* inpAr2, long numEl)
{ FILE *stream;
  stream = fopen(_fname, "wt");
  if (!stream) return false;
  for (long i=0; i<numEl; i++)
	fprintf(stream,"%0.3f\t%0.3f\n",inpAr1[i],inpAr2[i]);
  fclose(stream);
  return true;
}


// doesn't work
bool TIFPaletteWriter(const char* fnam, float* data, int dx, int dy, char* comment,
                      float minval, float maxval, char* lut)
{ char aChar;
  short aWord;
  int anInt;
  int commentLen = strlen(comment)+1;
  int colormapLen = 256; //?????
  FILE* tifF = fopen(fnam,"wb");
  if (tifF == NULL) return false;
  //Header
  fputc(0x49,tifF);                // I can check endians - steal from CASS
  fputc(0x49,tifF);
  aWord = 42;
  fwrite(&aWord,sizeof(aWord),1,tifF);
  anInt = 8; //offset
  fwrite(&anInt,sizeof(anInt),1,tifF);

  int numFields = 9;
  aWord = numFields; // Num of fields in IFD
  fwrite(&aWord,sizeof(aWord),1,tifF);

/*
  aWord = 254; //NewSubfileType ???
  fwrite(&aWord,sizeof(aWord),1,tifF);
  aWord = 4; // type int
  fwrite(&aWord,sizeof(aWord),1,tifF);
  anInt = 1; // num elements
  fwrite(&anInt,sizeof(anInt),1,tifF);
  anInt = 0;
  fwrite(&anInt,sizeof(anInt),1,tifF);
*/
  aWord = 256; // Width
  fwrite(&aWord,sizeof(aWord),1,tifF);
  aWord = 4; // type int
  fwrite(&aWord,sizeof(aWord),1,tifF);
  anInt = 1; // num elements
  fwrite(&anInt,sizeof(anInt),1,tifF);
  anInt = dx; // the width
  fwrite(&anInt,sizeof(anInt),1,tifF);

  aWord = 257; // Height
  fwrite(&aWord,sizeof(aWord),1,tifF);
  aWord = 4; // type int
  fwrite(&aWord,sizeof(aWord),1,tifF);
  anInt = 1; // num elements
  fwrite(&anInt,sizeof(anInt),1,tifF);
  anInt = dy; // the height
  fwrite(&anInt,sizeof(anInt),1,tifF);

  aWord = 258; // bits per sample
  fwrite(&aWord,sizeof(aWord),1,tifF);
  aWord = 3; // type short
  fwrite(&aWord,sizeof(aWord),1,tifF);
  anInt = 1; // num elements
  fwrite(&anInt,sizeof(anInt),1,tifF);
  aWord = 8; // 8bit?
  fwrite(&aWord,sizeof(aWord),1,tifF);
  aWord = 0; // 12 byte
  fwrite(&aWord,sizeof(aWord),1,tifF);

  aWord = 262; // Photometric Interpretation ???
  fwrite(&aWord,sizeof(aWord),1,tifF);
  aWord = 3; // type short
  fwrite(&aWord,sizeof(aWord),1,tifF);
  anInt = 1; // num elements
  fwrite(&anInt,sizeof(anInt),1,tifF);
  aWord = 3; // ?
  fwrite(&aWord,sizeof(aWord),1,tifF);
  aWord = 0; // 12 byte
  fwrite(&aWord,sizeof(aWord),1,tifF);

  aWord = 273; // offset to the data
  fwrite(&aWord,sizeof(aWord),1,tifF);
  aWord = 4; // type int
  fwrite(&aWord,sizeof(aWord),1,tifF);
  anInt = 1; // num elements
  fwrite(&anInt,sizeof(anInt),1,tifF);
  anInt = 14+numFields*12+commentLen+6*colormapLen;
  fwrite(&anInt,sizeof(anInt),1,tifF);

  aWord = 278; // rows per strip
  fwrite(&aWord,sizeof(aWord),1,tifF);
  aWord = 4; // type int
  fwrite(&aWord,sizeof(aWord),1,tifF);
  anInt = 1; // num elements
  fwrite(&anInt,sizeof(anInt),1,tifF);
  anInt = dx;//dx*dy*sizeof(char);
  fwrite(&anInt,sizeof(anInt),1,tifF);

  aWord = 279; // strip byte count
  fwrite(&aWord,sizeof(aWord),1,tifF);
  aWord = 4; // type int
  fwrite(&aWord,sizeof(aWord),1,tifF);
  anInt = 1; // num elements
  fwrite(&anInt,sizeof(anInt),1,tifF);
  anInt = dy;//dx*dy*sizeof(char);
  fwrite(&anInt,sizeof(anInt),1,tifF);

  aWord = 305; // Software
  fwrite(&aWord,sizeof(aWord),1,tifF);
  aWord = 2; // type string
  fwrite(&aWord,sizeof(aWord),1,tifF);
  anInt = commentLen; // the length
  fwrite(&anInt,sizeof(anInt),1,tifF);
  anInt = 14+numFields*12; // offset
  fwrite(&anInt,sizeof(anInt),1,tifF);

  aWord = 320; // Color Map
  fwrite(&aWord,sizeof(aWord),1,tifF);
  aWord = 3; // type short
  fwrite(&aWord,sizeof(aWord),1,tifF);
  anInt = 3*colormapLen; // num elements
  fwrite(&anInt,sizeof(anInt),1,tifF);
  anInt = 14+numFields*12+commentLen; // offset
  fwrite(&anInt,sizeof(anInt),1,tifF);

  anInt = 0; // next IFD - not present
  fwrite(&anInt,sizeof(anInt),1,tifF);

  //write comment
  for (int i=0; i<commentLen; i++)
    fputc(comment[i],tifF);
//  fprintf(tifF,"%s",comment);
  // Different LUTs

  //write LUT
  for (int i=0; i<colormapLen; i++)
  { aWord = i; //red
    fwrite(&aWord,sizeof(aWord),1,tifF);
    aWord = i; //green
    fwrite(&aWord,sizeof(aWord),1,tifF);
    aWord = i; //blue
    fwrite(&aWord,sizeof(aWord),1,tifF);
  }

  // Converting the data to the LUT

  //??? Should I make a mode for auto min-max? For example whem max<min

  // Just saving the data
  size_t numel = (size_t)dx * (size_t)dy;
  for (size_t i=0; i<numel; i++)
  { char aval = (char)data[i];
    fwrite(&aval,sizeof(char),1,tifF);
  }
  fclose(tifF);
  return true;
}

void BuildLUT(int numcomp, float* valL, int lutLen, char* aLUT, float gamma) // input from 0 to 1
{ int lim1 = 0;
  int lim2 = 0;//valL[0];
  for (int in=1; in<numcomp; in++)
  { lim1 = lim2;
    lim2 = roundMy(exp(log(fabs(valL[in*4]))/gamma)*(double)lutLen);
    if (in==numcomp-1) lim2++;
    for (int i=lim1; i<lim2; i++)
    { aLUT[3*i+0] = (char)roundMy(((i-lim1)*(valL[in*4+1]-valL[(in-1)*4+1])/(lim2-lim1)+valL[(in-1)*4+1]));
      aLUT[3*i+1] = (char)roundMy(((i-lim1)*(valL[in*4+2]-valL[(in-1)*4+2])/(lim2-lim1)+valL[(in-1)*4+2]));
      aLUT[3*i+2] = (char)roundMy(((i-lim1)*(valL[in*4+3]-valL[(in-1)*4+3])/(lim2-lim1)+valL[(in-1)*4+3]));
    }
  }

};

void FireLUT(int lutLen, char* aLUT, float gamma) // input from 0 to 1
{ int ncomp = 8;
  float* valL = new float[ncomp*4];
  valL[0*4+0] = 0; valL[0*4+1] = 0; valL[0*4+2] = 0; valL[0*4+3] = 0;
  valL[1*4+0] = 15/255.; valL[1*4+1] = 0; valL[1*4+2] = 0; valL[1*4+3] = 91;
  valL[2*4+0] = 56/255.; valL[2*4+1] = 122; valL[2*4+2] = 0; valL[2*4+3] = 227;
  valL[3*4+0] = 96/255.; valL[3*4+1] = 195; valL[3*4+2] = 0; valL[3*4+3] = 93;
  valL[4*4+0] = 127/255.; valL[4*4+1] = 238; valL[4*4+2] = 76; valL[4*4+3] = 0;
  valL[5*4+0] = 144/255.; valL[5*4+1] = 255; valL[5*4+2] = 117; valL[5*4+3] = 0;
  valL[6*4+0] = 208/255.; valL[6*4+1] = 255; valL[6*4+2] = 234; valL[6*4+3] = 0;
  valL[7*4+0] = 255/255.; valL[7*4+1] = 255; valL[7*4+2] = 255; valL[7*4+3] = 255;
  BuildLUT(ncomp, valL, lutLen, aLUT, gamma);
}

void BwLUT(int lutLen, char* aLUT, float gamma) // input from 0 to 1
{
  int ncomp = 2;
  float* valL = new float[ncomp*4];
  valL[0*4+0] = 0; valL[0*4+1] = 0; valL[0*4+2] = 0; valL[0*4+3] = 0;
  valL[1*4+0] = 1.; valL[1*4+1] = 255; valL[1*4+2] = 255; valL[1*4+3] = 255;
  BuildLUT(ncomp, valL, lutLen, aLUT, gamma);
}

void RainbowLUT(int lutLen, char* aLUT, float gamma) // input from 0 to 1
{
  int ncomp = 5;
  float* valL = new float[ncomp*4];
//  valL[0*4+0] = -1; valL[0*4+1] = 0; valL[0*4+2] = 0; valL[0*4+3] = 0; //???
  valL[0*4+0] = 0; valL[0*4+1] = 0; valL[0*4+2] = 0; valL[0*4+3] = 255;
  valL[1*4+0] = 63/255.; valL[1*4+1] = 0; valL[1*4+2] = 255; valL[1*4+3] = 255;
  valL[2*4+0] = 127/255.; valL[2*4+1] = 0; valL[2*4+2] = 255; valL[2*4+3] = 0;
  valL[3*4+0] = 191/255.; valL[3*4+1] = 255; valL[3*4+2] = 255; valL[3*4+3] = 0;
  valL[4*4+0] = 255/255.; valL[4*4+1] = 255; valL[4*4+2] = 0; valL[4*4+3] = 0;
  BuildLUT(ncomp, valL, lutLen, aLUT, gamma);
}

bool GIFPaletteWriter(const char* fnam, float* data, int dx, int dy, float minval, float maxval,
                      float gamma, int lutN, bool scalebar, bool includeHL)
{ char aChar;
  short aWord;

//  if (includeHL) fprintf(fnam,"%s
  FILE* gifF = fopen(fnam,"wb");
  if (gifF == NULL) return false;
  //Header      
  fputc('G',gifF); fputc('I',gifF); fputc('F',gifF);  //"GIF89a"
  fputc('8',gifF); fputc('9',gifF); fputc('a',gifF);
  aWord = dx;
  fwrite(&aWord,sizeof(aWord),1,gifF);
  aWord = dy;
  fwrite(&aWord,sizeof(aWord),1,gifF);
//  fputc(0xF7,gifF); // Color table for 256 colors
  fputc(0xF6,gifF); // Color table for 128 colors
  fputc(0x00,gifF); // BG color #0
  fputc(0x00,gifF); //  pixel aspect ratio
  // Now color table (LUT)
  int colormapLen = 127; //?????
  if (fabs(gamma)<MinVal) gamma = 1.;
  char* aLUT = new char[colormapLen*3];
  if (lutN==1) FireLUT(colormapLen, aLUT, gamma);
  else if (lutN=2) RainbowLUT(colormapLen, aLUT, gamma);
  else BwLUT(colormapLen, aLUT, gamma);

  for (int i=0; i<=colormapLen; i++)
    for (int j=0; j<3; j++)
    { aChar = aLUT[3*i+j]; //red
      fwrite(&aChar,sizeof(aChar),1,gifF);
    }
  delete aLUT;

// a scale bar
  short sdx = roundMy(0.02*dx);
  short sdy = roundMy(0.1*dy);
  short corx = dx-2*sdx;
  short cory = (includeHL?0.5*sdx:0);

  short* scaleAr = new short[sdx*(sdy+2*sdx)];
  for (int i=0; i<sdx*(sdy+2*sdx); i++) scaleAr[i] = -1;
  for (int iy=sdx; iy<sdy+sdx; iy++)
    for (int ix=0; ix<sdx; ix++)
      scaleAr[iy*sdx+ix] = (short)roundMy((sdy+sdx-iy)*colormapLen/(double)sdy);

  if (includeHL)   // integration of leters H and L in highest level color
  { int bord = roundMy(0.1*sdx);
    int wCh = roundMy(0.2*sdx);
    for (int iy=0; iy<sdx; iy++)
      for (int ix=0; ix<sdx; ix++)
        if ((((ix>bord && ix<bord+wCh) || (ix>sdx-bord-wCh-1 && ix<sdx-bord-1)) && (iy>=0 && iy<sdx-1)) ||
              ((iy>0.5*sdx-wCh/2-1 && iy<0.5*sdx+wCh/2-1) && (ix>bord && ix<sdx-bord-1)))
          scaleAr[iy*sdx+ix] = colormapLen;
        else scaleAr[iy*sdx+ix] = -1;

    for (int iy=sdy+sdx; iy<sdy+2*sdx; iy++)
      for (int ix=0; ix<sdx; ix++)
        if (((ix>bord && ix<bord+wCh) && (iy>sdy+sdx && iy<sdy+2*sdx)) ||
              ((iy>sdy+2*sdx-wCh && iy<sdy+2*sdx) && (ix>bord && ix<sdx-bord-1)))
          scaleAr[iy*sdx+ix] = colormapLen;
        else scaleAr[iy*sdx+ix] = -1;
  }

// max value
  int totnum = dx*dy;
  if (maxval<minval+MinVal)
    for (int i=0; i<totnum; i++)
      if (data[i]>maxval) maxval = data[i];
//  maxval *= 0.7;
  if (maxval<minval+MinVal) return false;
// image ----------------------------------
  // Array of the image data
  char* aPic = new char[totnum];
  for (int i=0; i<totnum; i++)
  { double aval = (float)colormapLen*(data[i]-minval)/maxval;
    if (aval<0) aPic[i] = 0;
    else if (aval>colormapLen) aPic[i] = colormapLen;
    else aPic[i] = (char)roundMy(aval);
  }

  // integration of color bar
  int minimsize = 50;
  if (scalebar && dx>minimsize && dy>minimsize)
    for (int iy=0; iy<sdy+2*sdx; iy++)
      for (int ix=0; ix<sdx; ix++)
        if (scaleAr[iy*sdx+ix]>=0) aPic[(iy+cory)*dx+ix+corx] = (char)scaleAr[iy*sdx+ix];
  delete scaleAr;

  fputc(0x2C,gifF); //  Image block
  aWord = 0; // corner X
  fwrite(&aWord,sizeof(aWord),1,gifF);
  aWord = 0; // corner Y
  fwrite(&aWord,sizeof(aWord),1,gifF);
  aWord = dx; // size X
  fwrite(&aWord,sizeof(aWord),1,gifF);
  aWord = dy; // size Y
  fwrite(&aWord,sizeof(aWord),1,gifF);
  fputc(0x00,gifF); // no local color table
//  fputc(0x08,gifF); // I think this is 8bit per pixel
  fputc(0x07,gifF); // I think this is 7bit per pixel

  // now image
  short sizebl = 126;//254
  for (int i=0; i<totnum; i+=sizebl)
  { short cursize = sizebl;
    if (totnum-i < sizebl) cursize = (char)(totnum-i);
    aChar = cursize+1;
    fwrite(&aChar,sizeof(aChar),1,gifF);
    fputc(0x80,gifF); // clear
    for (int j=0; j<cursize; j++)
    { if (i+j>=totnum) break;
      fwrite(&aPic[i+j],sizeof(aChar),1,gifF);
    }
  }
  fputc(0x01,gifF); // end of image
  fputc(0x81,gifF); // end of image
  fputc(0x00,gifF); // end of image
  delete aPic;

  fputc(0x3B,gifF); // end of file
  fclose(gifF);

  return true;
}
// addition to GIF
/*
  int minimsize = 50;
  if (dx>minimsize && dy>minimsize)
  {
  // graphic control extension
    fputc(0x21,gifF); // extension
    fputc(0xF9,gifF); // label
    fputc(0x04,gifF); // block size
    fputc(0x00,gifF); // method of disposal
    fputc(0xff,gifF); // time delay (word)
    fputc(0x00,gifF); //
    fputc(0x00,gifF); // transp. color
    fputc(0x00,gifF); // end

  // works as animated gif
    fputc(0x2C,gifF); //  Image block
    fwrite(&corx,sizeof(aWord),1,gifF); // corner X
    fwrite(&cory,sizeof(aWord),1,gifF); // corner Y
    fwrite(&sdx,sizeof(aWord),1,gifF);  // size X
    fwrite(&sdy,sizeof(aWord),1,gifF);  // size Y
    fputc(0x00,gifF); // no local color table
//  fputc(0x08,gifF); // I think this is 8bit per pixel
    fputc(0x07,gifF); // I think this is 7bit per pixel

    for (int i=0; i<totnum; i+=sizebl)
    { short cursize = sizebl;
      if (totnum-i < sizebl) cursize = (char)(totnum-i);
      aChar = cursize+1;
      fwrite(&aChar,sizeof(aChar),1,gifF);
      fputc(0x80,gifF); // clear
      for (int j=0; j<cursize; j++)
      { if (i+j>=totnum) break;
        char aval = scaleAr[i+j];
        if (aval<0) aChar = 0;
        else if (aval>colormapLen) aChar = colormapLen;
        else aChar = aval;
        fwrite(&aChar,sizeof(aChar),1,gifF);
      }
    }

    fputc(0x01,gifF); // end of image
    fputc(0x81,gifF); // end of image
    fputc(0x00,gifF); // end of image
*/
/*
// text block - doesn't work
  fputc(0x21,gifF); //  Image block
  fputc(0x01,gifF); //  Image block
  fputc(0x0C,gifF); //  Image block
  corx = 10;
  cory = 10;
  sdx = 50;
  sdy = 20;
  fwrite(&corx,sizeof(aWord),1,gifF); // corner X
  fwrite(&cory,sizeof(aWord),1,gifF); // corner Y
  fwrite(&sdx,sizeof(aWord),1,gifF);  // size X
  fwrite(&sdy,sizeof(aWord),1,gifF);  // size Y
  aChar = 5;// cell width
  fputc(aChar,gifF);
  aChar = 8;// cell height
  fputc(aChar,gifF);
  fputc((char)colormapLen,gifF); //foregraund
  fputc(0,gifF); //background
  char aString[255];
  strcpy(aString,"a Label 1 2");
  aChar = (char)strlen(aString);
  fputc(aChar,gifF);
  for (int i=0; i<strlen(aString); i++)
    fputc(aString[i],gifF);
  fputc(0x00,gifF); //end of block
  }
*/

int GetHDF5dims(char* fname, char* dataset, int** dims, int* numdims)
{
#ifdef Hdf5
  hid_t file_id = H5Fopen (fname, H5F_ACC_RDONLY, H5P_DEFAULT);
  if (file_id < 0)
  { printf("ERROR: File %s cannot be opened\n",fname);
    return 0;
  }

  if (!H5Lexists(file_id, dataset, H5P_DEFAULT ))
  { printf("ERROR: Could not open field %s\n",dataset);
    H5Fclose(file_id);
    return 0;
  }

  hid_t dat_id = H5Dopen2(file_id, dataset, H5P_DEFAULT);
  hid_t dataspace = H5Dget_space(dat_id);
  *numdims = H5Sget_simple_extent_ndims(dataspace);
  hsize_t* dimsL = new hsize_t[*numdims];
  H5Sget_simple_extent_dims(dataspace, dimsL, NULL);
  for(int i=0; i<*numdims; i++)
    (*dims)[i] = (int)dimsL[i];

  delete[] dimsL;
  H5Sclose(dataspace);
  H5Dclose(dat_id);
  H5Fclose(file_id);

  return 1;
#endif
  return 0;
}


int ReadHDF5frame(char* fname, char* dataset, float** outArray, int im_dim, int offset, int numframes, int** dims, int* numdims)
// im_dim - which dimention is the image number; offset - an offset to that position
{
#ifdef Hdf5
  hid_t file_id = H5Fopen (fname, H5F_ACC_RDONLY, H5P_DEFAULT);
  if (file_id < 0)
  { printf("ERROR: File %s cannot be opened\n",fname);
    return 0;
  }

  if (!H5Lexists(file_id, dataset, H5P_DEFAULT ))
  { printf("ERROR: Could not open field %s\n",dataset);
    return 0;
  }

  hid_t dat_id = H5Dopen2(file_id, dataset, H5P_DEFAULT);
  hid_t dataspace = H5Dget_space(dat_id);
  int ndims = H5Sget_simple_extent_ndims(dataspace);
  hsize_t* dimsL = new hsize_t[ndims];
  *numdims = ndims;
  H5Sget_simple_extent_dims(dataspace, dimsL, NULL);

  if (im_dim >=0 && im_dim < ndims)
    if (offset > dimsL[im_dim])
    { delete[] dimsL;
      printf("ERROR: Dimention %d is %d and it is smaller than offset %d!\n", im_dim, dimsL[im_dim], offset);
      return 0;
    }

  hsize_t* slabstart = new hsize_t[ndims];
  hsize_t* slabcount = new hsize_t[ndims];

  size_t fulldims = 1;

  //??? this is to read the whole files (3D for example)
//  if (offset < 0 && im_dim >= 0 && im_dim < ndims) 
//    offset = dimsL[im_dim];

//??? Not yet sure it is right - to check with other .h5 files. I think it is fine for 2D and 3D files only
//??? I think I need to sent exactly which axis is fs and which ss and make an image from those only
  if (im_dim >= 0 && im_dim < ndims) (*numdims)--;
  *dims = new int[*numdims];
  int k = 0;
  for (int i=0; i<ndims; i++)
  { //printf("%ld ",dimsL[i]);
    if (i != im_dim) 
    { slabstart[i] = 0;
      slabcount[i] = dimsL[i];
      if (dimsL[i] > 0)
        fulldims *= (size_t)dimsL[i];
      (*dims)[k] = (int)dimsL[i];
      k++;
    } else
    { slabstart[i] = offset;
      slabcount[i] = numframes;
      fulldims *= slabcount[i];
    }
  }
//?  if (im_dim >= 0 && im_dim < ndims) 
//?    (*dims)[*numdims-1] = (int)dimsL[im_dim];

  hid_t memspace = H5Screate_simple(ndims,slabcount,NULL);

  if (H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, slabstart, NULL, slabcount, NULL) < 0)
  { printf("ERROR: selecting hyperslab failed\n");
    return 0;
  }

  *outArray = new float[fulldims];
  float* outAr = *outArray;

  if (H5Dread(dat_id, H5T_NATIVE_FLOAT, memspace, dataspace, H5P_DEFAULT, outAr) < 0)
  { printf("ERROR: Cannot read data\n");
    return 0;
  }

  H5Sclose(memspace);
  H5Sclose(dataspace);
  H5Dclose(dat_id);
  H5Fclose(file_id);

  delete[] dimsL;
  delete[] slabstart;
  delete[] slabcount;

  return 1;
#endif
  return 0;
}

// run with im_dim=-1 to read the whole hdf5!
int ReadHDF5frame2D(char* fname, char* dataset, float** outArray, int im_dim, int offset, int* dims)
{
#ifndef Hdf5
  printf("Compiled without hdf5 library\n");
  return 0;
#endif

  int* dimsL;
  int numdims;
  if (ReadHDF5frame(fname, dataset, outArray, im_dim, offset, 1, &dimsL, &numdims) == 1)
  { if (numdims != 2)
    { printf("File %s is not a 2D or 3D file!\n",fname);
      delete[] dimsL;
      return 0;
    }
//    printf("Dims: %d, %d\n", dimsL[0], dimsL[1]);
    dims[0] = dimsL[0];
    dims[1] = dimsL[1];
    delete[] dimsL;
    return 1;
  }  
  return 0;
}

int WriteHDF5(char* fname, char* dataset, float* outArray, int* dims, int numDims)
{
#ifdef Hdf5
  hid_t file_id = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

  hsize_t* dimsL = new hsize_t[numDims];  
  for (int i=0; i<numDims; i++)
    dimsL[i] = dims[i];
  hid_t dataspace_id = H5Screate_simple(numDims, dimsL, NULL);

  hid_t datatype_id = H5Tcopy(H5T_NATIVE_FLOAT);
  hid_t status = H5Tset_order(datatype_id, H5T_ORDER_LE);

  hid_t dataset_id = H5Dcreate2(file_id, dataset, datatype_id, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, outArray);

  delete[] dimsL;
  H5Sclose(dataspace_id);
  H5Tclose(datatype_id);
  H5Dclose(dataset_id);
  H5Fclose(file_id);

  return 1;
#endif
  return 0;
}

//int ReadHDF5frame(char* fname, char* dataset, float** outArray, int im_dim, int offset, int numframes, int** dims, int* numdims)
/*
int ReadHDF5(char* fname, char* dataset, float** outArray, int** dims, int* numdims)
{
#ifdef Hdf5
  hid_t file_id = H5Fopen (fname, H5F_ACC_RDONLY, H5P_DEFAULT);
  if (file_id < 0)
  { printf("ERROR: File %s cannot be opened\n",fname);
    return 0;
  }

  if (!H5Lexists(file_id, dataset, H5P_DEFAULT ))
  { printf("ERROR: Could not open field %s\n",dataset);
    return 0;
  }

  hid_t dat_id = H5Dopen2(file_id, dataset, H5P_DEFAULT);
  hid_t dataspace = H5Dget_space(dat_id);
  int ndims = H5Sget_simple_extent_ndims(dataspace);
  hsize_t* dimsL = new hsize_t[ndims];
  *numdims = ndims;
  H5Sget_simple_extent_dims(dataspace, dimsL, NULL);
  hsize_t fuldims = 1;
  for (int i=0; i<*numdims; i++)
    fulldims *= dimsL[i];
  if (fulldims < 1)
  { delete[] dimsL;
    H5Sclose(dataspace);
    H5Dclose(dat_id);
    H5Fclose(file_id);
    printf("ERROR: fulldims == 0!\n");
    return 0;
  }
  *dims = new int[*numdims];
  for (int i=0; i<*numdims; i++)
    (*dims)[i] = (int)dimsL[i];

  *outArray = new float[fulldims];
  float* outAr = *outArray;



  if (H5Dread(dat_id, H5T_NATIVE_FLOAT, memspace, dataspace, H5P_DEFAULT, outAr) < 0)
  { printf("ERROR: Cannot read data\n");
    return 0;
  }

  H5Sclose(memspace);
  H5Sclose(dataspace);
  H5Dclose(dat_id);
  H5Fclose(file_id);

  delete[] dimsL;
  delete[] slabstart;
  delete[] slabcount;

  return 1;
#endif
  return 0;
}
*/