//#define _USE_MATH_DEFINES
//#include <cmath>
#include <math.h>
#include <ctype.h>
#include <stdarg.h>
#include <stdint.h>
#include <algorithm>
//#include "peakfinder8.cpp"
#include "fileformats.h" 

//_control87(MCW_EM, MCW_EM);

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

//const int MaxStrLen = 1020;
const float BstpReg = -10001.;
const int MASK_BAD = -3;
const int MASK_GOOD = 1;
const char MASK_BRAGG = 2;
//const float MinVal = 1e-10;

#define MIN(A,B) (A<B?A:B)
#define MAX(A,B) (A>B?A:B)
#define SQR(A) ((A)*(A))
#define AMP(A,B) sqrt(A*A+B*B)
#define INTENS(A,B) (A*A+B*B)
#define INTENSCPX(A) (A[0]*A[0]+A[1]*A[1])

//#define alglib          // some STL library, rather useful
#ifdef alglib
#include "specialfunctions.cpp"
#endif

#define double2int(i, d) \
    {double t = ((d) + 6755399441055744.0); i = *((int *)(&t));}

#define double2intA(i,n)  __asm {__asm fld n   __asm fistp i}

//const double MinVal = 1e-10;
//const double MaxVal = 1e+30;
//const int MaxStrLen = 255;

//const char MASK_BAD = 0;
//const char MASK_GOOD = 1;               
//const char MASK_BRAGG = 2;


#define SWAP(A,B) {A=A^B;B=A^B;A=A^B;}
void SwapFloat(float* A, float* B)
{ float T = *A;
  *A = *B;
  *B = T;
}

void BubbleSortFloat(float* A, int* B, int n)
{ for (int i=n-1; i>0; i--)
	for (int j=0; j<i; j++)
	  if (A[j] > A[j+1])
	  { SwapFloat(&A[j],&A[j+1]);
		SWAP(B[j],B[j+1]);
	  }
};

void InvBubbleSortFloat(float* A, int* B, int n)
{ for (int i=n-1; i>0; i--)
	for (int j=0; j<i; j++)
	  if (A[j] < A[j+1])
	  { SwapFloat(&A[j],&A[j+1]);
		SWAP(B[j],B[j+1]);
	  }
};

#undef SWAP

bool PresentInArrayInt(int digit, int* ar, int numel)
{ for (int i=0; i<numel; i++)
	if (ar[i] == digit) return true;
  return false;
}


double sqr(double inp)
{ return inp*inp;
};

void rotate2D(float* x1, float* y1, float ang)
{ float xin = *x1;
  *x1 = xin*cos(ang)-*y1*sin(ang);
  *y1 = xin*sin(ang)+*y1*cos(ang);
}

void rotate2Dd(double* x1, double* y1, double ang)
{ double xin = *x1;
  *x1 = xin*cos(ang)-*y1*sin(ang);
  *y1 = xin*sin(ang)+*y1*cos(ang);
}

int round1(double inp)
{ double dif = inp - (int)inp;
  if (fabs(dif)>=0.5) return (int)inp + (dif>0.?1:-1);//fabs(dif);
  else return (int)inp;
};


double powerMy(double x, double y)
{ return 1.;
};

int Sign(double inp)
{ return (inp>MinVal?1:(inp<-MinVal?-1:0));
};

bool Near(double inp1, double inp2)
{ return (inp1+1e-10>inp2&&inp1-1e-10<inp2);
};

bool NearEps(double inp1, double inp2, double eps)
{ return (inp1+eps>inp2&&inp1-eps<inp2);
};

float BilenearInterp(float cx, float cy, float fx0y0, float fx0y1, float fx1y0, float fx1y1)
{ return fx0y0*(1-cx)*(1-cy)+fx1y0*cx*(1-cy)+fx0y1*(1-cx)*cy+fx1y1*cx*cy;
};

float BilenearInterp1(float cx, float cy, float fx0y0, float fx0y1, float fx1y0, float fx1y1)
{
 int cxup = roundMy(cx+0.5-MinVal);
 int cxdn = roundMy(cx-0.5+MinVal);
 int cyup = roundMy(cy+0.5-MinVal);
 int cydn = roundMy(cy-0.5+MinVal);
 return fx0y0*(cxup-cx)*(cyup-cy)+fx1y0*(cx-cxdn)*(cyup-cy)+fx0y1*(cxup-cx)*(cy-cydn)+fx1y1*(cx-cxdn)*(cy-cydn);
};

double MasMax(int Numb, double* inp)
{ double tempMax = -1e300;
  for (int i=0; i<Numb; i++)
	if (inp[i]>tempMax) tempMax = inp[i];
  return tempMax;
};

float MasMaxFl(int Numb, float* inp)
{ float tempMax = -1e35;
  for (int i=0; i<Numb; i++)
	if (inp[i]>tempMax) tempMax = inp[i];
  return tempMax;
};

int NOD(int x, int y)
// greatest common divisor = NOD in russian
{ if (x!=0) return NOD(y % x,x);
  else return y;
};

#define SWAP(A,B) {A=A^B;B=A^B;A=A^B;}
void BubbleSort(int A[], int n)
{ for (int i=n-1; i>0; i--)
	for (int j=0; j<i; j++)
	  if (A[j] > A[j+1]) SWAP(A[j],A[j+1]);
};
#undef SWAP

bool root2(double a, double b, double c, double* x1, double* x2)
{
  double D = b*b - 4.*a*c;
  double MinVal = 1e-10;
  if (fabs(D)<MinVal) D=0;
  if (D >=0)
	if (fabs(a)>MinVal)
	{ *x1 = 0.5*(-b+sqrt(D))/a;
	  *x2 = 0.5*(-b-sqrt(D))/a;
	  return true;
	} else if (fabs(b)>MinVal)
	  { *x1 = -c/b;
		*x2 = -c/b;
		return true;
	  } else return false;
  else if (fabs(b)<MinVal && fabs(a)>MinVal)
  { if (-c/a>-MinVal)
	{ if (fabs(c/a)<MinVal) c = 0;
	  *x1 = sqrt(-c/a);
	  *x2 = -*x1;
		  return true;
	} else return false;
  } else return false;
  return true;
};

// Equations from wolfram website f(ab) = kx + b          /* float badVal,*/
float LeastSquaresFit(float* xAr, float* yAr, size_t numEl, float* k, float* b)
{ if (numEl<1) return 0.;
  float xAv = 0;
  float yAv = 0;
  float x2Su = 0;
  float y2Su = 0;
  float xySu = 0;

  for (size_t i=0; i<numEl; i++)
//  if (xAr[i]>badVal+1 && yAr[i]>badVal+1)
    { xAv += xAr[i];
      yAv += yAr[i];
      x2Su += SQR(xAr[i]);
      y2Su += SQR(yAr[i]);
      xySu += xAr[i]*yAr[i];
      // counter++; // for badval
    }
  xAv /= (float)numEl;
  yAv /= (float)numEl;

  float ssxx = x2Su - numEl*SQR(xAv);
  float ssyy = y2Su - numEl*SQR(yAv);
  float ssxy = xySu - numEl*xAv*yAv;

  if (fabs(ssxx)<MinVal || fabs(ssyy)<MinVal) return 0.;

  *k = ssxy/ssxx;
  *b = yAv - *k * xAv;

  return SQR(ssxy)/ssxx/ssyy;
};



#define elem_type float
elem_type quick_select(elem_type arr[], size_t n)
{   if (n<1) return 0;
    size_t low, high, median, middle, ll, hh;
    elem_type A;

    low = 0 ; high = n-1 ; median = (low + high) / 2;
  for (;;) {
        if (high <= low) /* One element only */
            return arr[median] ;

        if (high == low + 1) {  /* Two elements only */
            if (arr[low] > arr[high])
                {A=arr[low]; arr[low]=arr[high]; arr[high]=A;}
            return arr[median] ;
        }

    /* Find median of low, middle and high items; swap into position low */
    middle = (low + high) / 2;
    if (arr[middle] > arr[high]) {A=arr[middle]; arr[middle]=arr[high]; arr[high]=A;}
    if (arr[low] > arr[high])    {A=arr[low]; arr[low]=arr[high]; arr[high]=A;}
    if (arr[middle] > arr[low])  {A=arr[middle]; arr[middle]=arr[low]; arr[low]=A;}

    /* Swap low item (now in position middle) into position (low+1) */
    {A=arr[middle]; arr[middle]=arr[low+1]; arr[low+1]=A;}

    /* Nibble from each end towards middle, swapping items when stuck */
    ll = low + 1;
    hh = high;
    for (;;) {
        do ll++; while (arr[low] > arr[ll]) ;
        do hh--; while (arr[hh]  > arr[low]) ;

        if (hh < ll)
        break;

        {A=arr[ll]; arr[ll]=arr[hh]; arr[hh]=A;}
    }

    /* Swap middle item (in position low) back into correct position */
    {A=arr[low]; arr[low]=arr[hh]; arr[hh]=A;}

    /* Re-set active partition */
    if (hh <= median)
        low = ll;
    if (hh >= median)
        high = hh - 1;
  }
}
#undef elem_type

float medianCutoff(float arr[], size_t n, float cutoff)  // cutoff from 0 to 1
{ if (n<1) return 0;
  const float maxval = 1e30;
  int upperlim = roundMy(cutoff*((double)n+0.01));
//  if (upperlim<1) upperlim=1;
  float _minim=arr[0];
  for (int jj=0; jj<upperlim; jj++)
  { _minim = maxval;
    int cmin = 0;
    for (int j=0; j<n; j++)
      if (arr[j]<_minim)
      { _minim = arr[j];
        cmin = j;
      }
    arr[cmin] = maxval;
  }
  return _minim;
}

float medianCutoffNew(float arr[], size_t n, float cutoff)  // cutoff from 0 to 1
{ if (n<1) return 0;
  const float maxval = 1e30;
  int upperlim = (int)(cutoff*n);
  if (upperlim<1) upperlim=1;
  float _minim=0;
  for (int jj=0; jj<upperlim; jj++)
  { _minim = arr[0];
    int cmin = 0;
    for (int j=1; j<n; j++)
      if (arr[j]<_minim)
      { _minim = arr[j];
        cmin = j;
      }
    arr[cmin] = maxval;
  }
  return _minim;
}

float medianCutoffStd(float arr[], size_t n, float cutoff)  // cutof
{
    if (n < 1)
        return 0;

    int valueToSelect = roundMy(cutoff * ((float) n + 0.01)) - 1;
    std::nth_element(arr, arr + valueToSelect, arr + n);
    return arr[valueToSelect];
}

float stdDev(float* arr, size_t N, float badVal)
{ float aver = 0;
  float aver2 = 0;
  size_t numAdd = 0;
  for (size_t i=0; i<N; i++)
    if (arr[i]>badVal)
    { aver += arr[i];
      aver2 += arr[i]*arr[i];
      numAdd++;
    }
  if (numAdd>0)
  { aver /= (float)numAdd;
    aver2 /= (float)numAdd;
  }
  return sqrt(aver2 - aver*aver);
};

float average(float* arr, size_t N, float badVal)
{ float aver = 0;
  size_t numAdd = 0;
  for (size_t i=0; i<N; i++)
    if (arr[i]>badVal)
    { aver += arr[i];
      numAdd++;
    }
  if (numAdd>0)
    aver /= (float)numAdd;
  return aver;
};

void moment(float* arr, size_t N, float badVal, float result[6])
{ double aver = 0;
  double aver2 = 0;
  double aver3 = 0;
  size_t numAdd = 0;
  for (size_t i=0; i<N; i++)
    if (arr[i]>badVal)
    { aver += arr[i];
      aver2 += arr[i]*arr[i];
      aver3 += aver2*arr[i];
      numAdd++;
    }
  if (numAdd>0)
  { aver /= (double)numAdd;
    aver2 /= (double)numAdd;
    aver3 /= (double)numAdd;
  }
  result[0] = aver;               //mean
  result[1] = aver2 - aver*aver;  //variance
  result[5] = sqrt(result[1]);    //stddev
  if (fabs(result[1]*result[5])>MinVal)
    result[2] = (aver3 - 3*aver*result[1] - aver*aver*aver)/(result[1]*result[5]); //skewness
  else
    result[2] = 0;
  result[3] = 0; //kurtosis
  result[4] = 0; //mean absolute
};

void moment1(float* arr, size_t N, float badVal, float result[6])
{ double aver = 0;
  double aver2 = 0;
  double aver3 = 0;
  size_t numAdd = 0;
  for (size_t i=0; i<N; i++)
    if (arr[i]>badVal)
    { aver += arr[i];
      numAdd++;
    }
  if (numAdd>0)
  { aver /= (double)numAdd;
  }

  // now variance and skew:
  for (size_t i=0; i<N; i++)
    if (arr[i]>badVal)
    { double sq = SQR(arr[i] - aver);
      aver2 += sq;
      aver3 += sq*(arr[i] - aver);
    }
  if (numAdd>1)
  {
    aver2 /= (double)(numAdd-1.);
    aver3 /= (double)(numAdd);
  }

  result[0] = aver;               //mean
  result[1] = aver2;              //variance
  result[5] = sqrt(result[1]);    //stddev
  if (fabs(result[1]*result[5])>MinVal)
    result[2] = (aver3)/(result[1]*result[5]); //skewness
  else
    result[2] = 0;
  result[3] = 0; //kurtosis
  result[4] = 0; //mean absolute
};


// ----------------------------- STRINGS -----------------------------------

void ErrorTerm(char* errorMessage)
{ printf("%s\n",errorMessage);
  exit(0);
};

void TrimSp(char* str1)
{ char str2[255];
  int po = 0;
  for (int i=0; i<strlen(str1); i++)
	if (!isspace(str1[i]))
	{ str2[po] = str1[i];
	  po++;
	}
  str2[po] = 0;
  strcpy(str1,str2);
}

void TrimNoCh(char* str1)
{ //char str2[MaxStrLen];
  int cstlen = strlen(str1);
  for (int i=cstlen-1; i>0; i--)
//	if (!isalnum(str1[i]) && str1[i]!='-')
	if (!isalnum(str1[i]) && str1[i]!='-' && str1[i] != '/' && str1[i] != '\\' && str1[i] != '.' && str1[i] != '+' && str1[i] != '\%')
	  str1[i] = 0;
	else break;
  while ((!isalnum(str1[0]) && str1[0] != '/' && str1[0] != '\\' && str1[0] != '.' && str1[0] != '-' && str1[0] != '+' && str1[0] != '\%') && strlen(str1)>0)
	for (int i=0; i<strlen(str1); i++)
	  str1[i] = str1[i+1];
}

void strLower(char *out, char* inp)
{ for (int i=0; i<strlen(inp); i++)
	out[i] = tolower(inp[i]);
}

void strLower1(char *inp)
{ for (int i=0; i<strlen(inp); i++)
	inp[i] = tolower(inp[i]);
}

char** strlow1(char* inp)  //bad - returning pointer to local variable!
{ char* out = new char[255];
  for (int i=0; i<=strlen(inp); i++)
	out[i] = tolower(inp[i]);
  return &out;
}

// output to a file and to the screen
int myprintf(FILE* afile, const char * format, ... )
{// Here you can redfine your input before continuing to compy with standard inputs
    char astr[MaxStrLen];
    va_list args;
    va_start(args, format);
    vsprintf(astr,format, args);// This still uses standaes formating
    va_end(args);
    printf(astr);
    if (afile) fprintf(afile,astr);
    return 0;// Before return you can redefine it back if you want...
}

// ----------------------------- FILES -----------------------------------

char* ExtractFileNameC(char* fullName)
{ return (strrchr(fullName,'\\')!=NULL?strrchr(fullName,'\\')+1:(strrchr(fullName,'/')!=NULL?strrchr(fullName,'/')+1:fullName));
//{ char* sss = new char[255];
//  strcpy(sss, fullName);
//  if (strrchr(fullName,'\\')!=NULL) sss = strrchr(fullName,'\\')+1;
//  if (strrchr(sss,'/')!=NULL) sss = strrchr(fullName,'/')+1;
//  if (strrchr(fullName,'\\')!=NULL) sss = strrchr(fullName,'\\')+1;
//  return sss;
}

char* ExtractFileNameNoExtC(char* fullName)
{ char* out = new char[255];
  strcpy(out,fullName);
  if (strrchr(ExtractFileNameC(out),'.') != NULL) *strrchr(ExtractFileNameC(out),'.') = 0;
  return out;
}

char* ExtractJustFileNameNoExtC(char* fullName)
{ char* out = new char[255];
  strcpy(out,ExtractFileNameC(fullName));
  if (strrchr(out,'.') != NULL) *strrchr(out,'.') = 0;
  return out;
}

void ExtractFileNameOnly(char* fullName)   //BAD!
{ //fullName = ExtractFileNameC(fullName); 
  char* sss;
  if (strrchr(fullName,'\\') != NULL) sss = strrchr(fullName,'\\')+1;
  else sss = fullName;
  if (strrchr(sss,'/') != NULL) fullName = strrchr(sss,'/')+1;
  else fullName = sss;
  printf("fullName = %s\n",fullName);
  if (strrchr(fullName,'.') != NULL) *strrchr(fullName,'.') = 0;  
  printf("fullName = %s\n",fullName);
}

char* ExtractFilePathC(char* fullName)
{ char* out = new char[255];
  strcpy(out,fullName);
  if (strrchr(out,'\\') != NULL) *strrchr(out,'\\') = 0;
  else if (strrchr(out,'/') != NULL) *strrchr(out,'/') = 0;
  return out;
}

/*
char** ExtractFileNameNoExtC1(char* fullName)  //bad - returning pointer to local variable!
{ char* out = new char[255];
  strcpy(out,fullName);
  if (strchr(ExtractFileNameC(out),'.') != NULL) *strchr(ExtractFileNameC(out),'.') = 0;
  return &out;
}
*/
const char* ExtractFileExtC(char* fullName)
{ return (strrchr(ExtractFileNameC(fullName),'.')==NULL ? "" : *strlow1(strrchr(ExtractFileNameC(fullName),'.')));
}

void ExtractFilePathC(char* fullName, char* fullPath)
{ int _n = strlen(fullName)-strlen(ExtractFileNameC(fullName));
  strcpy(fullPath, fullName);
  fullPath[_n] = 0;
}

bool fileExists(char* fName)
{ FILE *_fileEx = fopen(fName, "rb");
  if (_fileEx != NULL)
  { fclose(_fileEx);
	return true;
  }
  return false;
}

int CheckInpFileDimsCubeFloat(char* fnam, int dims)
{ int _dim = 0;
  FILE *stream;
  stream = fopen(fnam, "rb");
  if (!stream) return 0;
  fseek (stream, 0, 2 );
  size_t flen = ftell(stream)/4;
  fclose(stream);
  if (dims=3)
  { int trydim = roundMy(pow(flen,1./3.));
    if (flen % trydim == 0)
	{ _dim = trydim;
	} else
	  printf("The input file is not 3D cube!\n");
  } else // check for square
  { int trydim = roundMy(sqrt(flen));
    if (flen % trydim == 0)
	{ _dim = trydim;
	} else
	  printf("The input file is not 2D square!\n");
  }
  return _dim;
};

#include <stdio.h>
long FileSize(char* filnam)
{
   FILE *stream;
//   stream = fopen("250.txt", "w+");
   stream = fopen(filnam, "r+");
   long curpos, length;
   curpos = ftell(stream);
   fseek(stream, 0L, SEEK_END);
   length = ftell(stream);
   fseek(stream, curpos, SEEK_SET);
   fclose(stream);
   return length;
}

#include <fstream>
int FileSize1(const char* sFileName)
{
  std::ifstream f;
  f.open(sFileName, std::ios_base::binary | std::ios_base::in);
  if (!f.good() || f.eof() || !f.is_open()) { return 0; }
  f.seekg(0, std::ios_base::beg);
  std::ifstream::pos_type begin_pos = f.tellg();
  f.seekg(0, std::ios_base::end);
  return static_cast<int>(f.tellg() - begin_pos);
}


//------------------------ MATRIXES 2x2 ----------------------------------------

void MultMatrix2(double InpMatrix1[2][2], double InpMatrix2[2][2], double OutMatrix[2][2])
{ for (int i=0; i<=1; i++)
	for (int k=0; k<=1; k++)
	  OutMatrix[i][k] = 0.0;
  for (int i=0; i<=1; i++)
	for (int k=0; k<=1; k++)
	  for (int j=0; j<=1; j++)
		OutMatrix[i][k] = OutMatrix[i][k] + InpMatrix1[i][j]*InpMatrix2[j][k];
};

void AddMatrix2(double InpMatrix1[2][2], double InpMatrix2[2][2], double OutMatrix[2][2], int _sub)
{ double sub;
  if (_sub == -1) sub = -1.;
  else sub = 1.;
  for (int i=0; i<=1; i++)
	for (int j=0; j<=1; j++)
		OutMatrix[i][j] = InpMatrix1[i][j] + sub*InpMatrix2[i][j];
};

void CopyMatrix2(double InpMatrix[2][2], double OutMatrix[2][2])
{ for (int i=0; i<=1; i++)
    for (int k=0; k<=1; k++)
        OutMatrix[i][k] = InpMatrix[i][k];
};

//double** InvertMatrix(double InpMatrix[2][2])
void InvertMatrix2(double InpMatrix[2][2], double OutMatrix[2][2])
{ double InpDet = InpMatrix[1][1]*InpMatrix[0][0]-InpMatrix[1][0]*InpMatrix[0][1];
  OutMatrix[0][0] = InpMatrix[1][1]/InpDet;
  OutMatrix[0][1] = -InpMatrix[0][1]/InpDet;
  OutMatrix[1][0] = -InpMatrix[1][0]/InpDet;
  OutMatrix[1][1] = InpMatrix[0][0]/InpDet;
};


// ------------------ Matrix (other) ---------------

//| a11 a12 a13 |-1             |   a33a22-a32a23  -(a33a12-a32a13)   a23a12-a22a13  |
//| a21 a22 a23 |    =  1/DET * | -(a33a21-a31a23)   a33a11-a31a13  -(a23a11-a21a13) |
//| a31 a32 a33 |               |   a32a21-a31a22  -(a32a11-a31a12)   a22a11-a21a12  |
//with DET  =  a11(a33a22-a32a23)-a21(a33a12-a32a13)+a31(a23a12-a22a13)
//| inp[0] inp[1] inp[2] |-1             |   inp[8]inp[4]-inp[7]inp[5]  -(inp[8]inp[1]-inp[7]inp[2])   inp[5]inp[1]-inp[4]inp[2]  |
//| inp[3] inp[4] inp[5] |    =  1/DET * | -(inp[8]inp[3]-inp[6]inp[5])   inp[8]inp[0]-inp[6]inp[2]  -(inp[5]inp[0]-inp[3]inp[2]) |
//| inp[6] inp[7] inp[8] |               |   inp[7]inp[3]-inp[6]inp[4]  -(inp[7]inp[0]-inp[6]inp[1])   inp[4]inp[0]-inp[3]inp[1]  |
//with DET  =  inp[0](inp[8]inp[4]-inp[7]inp[5])-inp[3](inp[8]inp[1]-inp[7]inp[2])+inp[6](inp[5]inp[1]-inp[4]inp[2])
bool InvertMatrix3(float* inp, float* out)
{ double det = inp[0]*(inp[8]*inp[4]-inp[7]*inp[5])-inp[3]*(inp[8]*inp[1]-inp[7]*inp[2])+inp[6]*(inp[5]*inp[1]-inp[4]*inp[2]);
  if (fabs(det)<1e-10) return false;
  det = 1./det;
  out[0] =  (inp[8]*inp[4]-inp[7]*inp[5])*det;
  out[1] = -(inp[8]*inp[1]-inp[7]*inp[2])*det;
  out[2] =  (inp[5]*inp[1]-inp[4]*inp[2])*det;
  out[3] = -(inp[8]*inp[3]-inp[6]*inp[5])*det;
  out[4] =  (inp[8]*inp[0]-inp[6]*inp[2])*det;
  out[5] = -(inp[5]*inp[0]-inp[3]*inp[2])*det;
  out[6] =  (inp[7]*inp[3]-inp[6]*inp[4])*det;
  out[7] = -(inp[7]*inp[0]-inp[6]*inp[1])*det;
  out[8] =  (inp[4]*inp[0]-inp[3]*inp[1])*det;
  return true;
};

float DetMatrix3(float* inp)
{ return inp[0]*(inp[8]*inp[4]-inp[7]*inp[5])-inp[3]*(inp[8]*inp[1]-inp[7]*inp[2])+inp[6]*(inp[5]*inp[1]-inp[4]*inp[2]);
};

int CopyMatrix(float* InpMatrix1, float* OutMatrix, int num)
{ for (int i=0; i<num; i++)
	OutMatrix[i] = InpMatrix1[i];
  return 1;
};

void TransposeArray(int* d0, int* d1, float* im)
{ float *tmpAr = new float[*d1**d0];
  for(int i = 0; i<*d0**d1; i++)
	  tmpAr[i] = im[i];

  for(int j = 0; j<*d1; j++)
	  for(int k = 0; k<*d0; k++)
	    im[j+*d1*k] = tmpAr[k+*d0*j];
  int tmp = *d0;
  *d0 = *d1;
  *d1 = tmp;

  delete[] tmpAr;
};


bool MultMatrixF(float* inp1, int n, int m, float* inp2, int m1, int l, float* out)
{ if (m1 != m) return false;
  int i,j,k;
	for(i=0; i<n; i++)
	  for(j=0; j<l; j++)
	  { out[j+l*i] = 0.;
		for(k=0; k<m; k++)
		  out[j+l*i] = out[j+l*i]+inp1[k+m*i]*inp2[j+l*k];
	  }
  return true;
}

bool MultMatrixF1(float* inp1, int n, int m, float* inp2, int m1, int l, float* out)
{ if (m1 != m) return false;
  double* outD = new double[l*n];
  int i,j,k;
	for(i=0; i<n; i++)
	  for(j=0; j<l; j++)
	  { outD[j+l*i] = 0.;
		for(k=0; k<m; k++)
		  outD[j+l*i] += inp1[k+m*i]*inp2[j+l*k];
	  }
  for (int i=0; i<l*n; i++)
	out[i] = outD[i];
  delete outD;
  return true;
}

bool MultMatrix(double* inp1, int n, int m, double* inp2, int m1, int l, double* out)
{ if (m1 != m) return false;
  int i,j,k;
	for(i=0; i<n; i++)
	  for(j=0; j<l; j++)
	  { out[j+l*i] = 0.;
		for(k=0; k<m; k++)
		  out[j+l*i] = out[j+l*i]+inp1[k+m*i]*inp2[j+l*k];
	  }
  return true;
}

int MultMatrixBad(double* InpMatrix1, int m, int n, double* InpMatrix2, int n1, int r,
			   double* OutMatrix)   // Fignja nerabochaja
{ if (n1 != n) return -1;
  for (int i=1; i<=m; i++)
	for (int k=1; k<=r; k++)
	  OutMatrix[i,k] = 0;
  for (int i=1; i<=m; i++)
	for (int k=1; k<=r; k++)
	  for (int j=1; j<=n; j++)
		OutMatrix[i,k] = OutMatrix[i,k] + InpMatrix1[i,j]*InpMatrix2[j,k];
  return 1;
};

int AddMatrix(double* InpMatrix1, int m, int n, double* InpMatrix2,
			  double* OutMatrix, int _sub)    // Fignja nerabochaja
{ double sub;
  if (_sub == -1) sub = -1.;
  else sub = 1.;
  for (int i=0; i<m; i++)
	for (int j=0; j<n; j++)
	  OutMatrix[i*m+j] = InpMatrix1[i*m+j] + sub*InpMatrix2[i*m+j];
  return 1;
};

// --------------------------- VECTORS -----------------------------------

int Rotating(double angle, double AinX,   double AinY,   double AinZ,
						   double BX,     double BY,     double BZ,
						   double* AoutX, double* AoutY, double* AoutZ)
{
  double AbsB = sqrt(BX*BX+BY*BY+BZ*BZ);
  if (AbsB<MinVal) return -1;
  double Skal = (1.-cos(angle))*(BX*AinX+BY*AinY+BZ*AinZ)/AbsB;
  *AoutX = AinX*cos(angle)+(BX*Skal+(BY*AinZ-BZ*AinY)*sin(angle))/AbsB;
  *AoutY = AinY*cos(angle)+(BY*Skal+(BZ*AinX-BX*AinZ)*sin(angle))/AbsB;
  *AoutZ = AinZ*cos(angle)+(BZ*Skal+(BX*AinY-BY*AinX)*sin(angle))/AbsB;
  return 1.;
};

//  / cos(phi)  sin(phi) 0 \ / 1      0         0    \ /  cos(psi)  sin(psi) 0 \
//  |-sin(phi)  cos(phi) 0 | | 0  cos(the)  sin(the) | | -sin(psi)  cos(psi) 0 |
//  \    0         0     1 / \ 0 -sin(the)  cos(the) / \     0         0     1 /
// Rotation clockwize
int RotateEuler1(double phi,    double the,    double psi,
				 double AX,     double AY,     double AZ,
				 double* AoutX, double* AoutY, double* AoutZ)
{ double rotM1[9], rotM2[9], rotM3[9], resM[9], inpC[3], outC[3];
  for (int i=0; i<9; i++)
  { rotM1[i] = 0.;
	rotM2[i] = 0.;
	rotM3[i] = 0.;
  }
  rotM1[8] =  rotM2[0] = rotM3[8] = 1.;
  rotM1[0] =  rotM1[4] = cos(phi);
  rotM1[1] = sin(phi); rotM1[3] = -rotM1[1];
  rotM2[4] =  rotM2[8] = cos(the);
  rotM2[5] = sin(the); rotM2[7] = -rotM2[5];
  rotM3[0] =  rotM3[4] = cos(psi);
  rotM3[1] = sin(psi); rotM3[3] = -rotM3[1];
  MultMatrix(rotM1, 3, 3, rotM2, 3, 3, resM);
  MultMatrix(resM, 3, 3, rotM3, 3, 3, rotM1);
  inpC[0] = AX; inpC[1] = AY; inpC[2] = AZ;
  MultMatrix(rotM1, 3, 3, inpC, 3, 1, outC);
  *AoutX = outC[0]; *AoutY = outC[1]; *AoutZ = outC[2];
  return 1;
};

int RotateEuler(double phi,    double the,    double psi,
				double AX,     double AY,     double AZ,
				double* AoutX, double* AoutY, double* AoutZ)
{
//(72�, 	0�, 	0�) 	= 	(40�, 	0�, 	32�) 	- singular alignment
//(45�, 	60�, 	-30�) 	= 	(-135�, 	-60�, 	150�) = (225, 300, ?330) 	- bistable flip
  double resM[9], inpC[3], outC[3];
  double c1,c2,c3,s1,s2,s3;
  c1=cos(the);
  c2=cos(psi);
  c3=cos(phi);
  s1=sin(the);
  s2=sin(psi);
  s3=sin(phi);

  resM[0]=c2*c3-c1*s2*s3;
  resM[1]=-c2*s3-c1*s2*c3;
  resM[2]=s1*s2;
  resM[3]=s2*c3+c1*c2*s3;
  resM[4]=-s2*s3+c1*c2*c3;
  resM[5]=-s1*c2;
  resM[6]=s1*s3;
  resM[7]=s1*c3;
  resM[8]=c1;

  inpC[0] = AX; inpC[1] = AY; inpC[2] = AZ;
  MultMatrix(inpC, 1, 3, resM, 3, 3, outC);
//  MultMatrix(resM, 3, 3, inpC, 3, 1, outC);
  *AoutX = outC[0]; *AoutY = outC[1]; *AoutZ = outC[2];
  return 1;
};

int RotEulerMatrix(float phi, float the, float psi, float* resM)
{ double c1,c2,c3,s1,s2,s3;
  c1=cos(the);
  c2=cos(psi);
  c3=cos(phi);
  s1=sin(the);
  s2=sin(psi);
  s3=sin(phi);
  resM[0]=c2*c3-c1*s2*s3;
  resM[1]=-c2*s3-c1*s2*c3;
  resM[2]=s1*s2;
  resM[3]=s2*c3+c1*c2*s3;
  resM[4]=-s2*s3+c1*c2*c3;
  resM[5]=-s1*c2;
  resM[6]=s1*s3;
  resM[7]=s1*c3;
  resM[8]=c1;
  return 1;
};

int RotEulerMatrixD(double phi, double the, double psi, double* resM)
{ double c1,c2,c3,s1,s2,s3;
  c1=cos(the);
  c2=cos(psi);
  c3=cos(phi);
  s1=sin(the);
  s2=sin(psi);
  s3=sin(phi);
  resM[0]=c2*c3-c1*s2*s3;
  resM[1]=-c2*s3-c1*s2*c3;
  resM[2]=s1*s2;
  resM[3]=s2*c3+c1*c2*s3;
  resM[4]=-s2*s3+c1*c2*c3;
  resM[5]=-s1*c2;
  resM[6]=s1*s3;
  resM[7]=s1*c3;
  resM[8]=c1;
  return 1;
};

int SmallAnglesRotMatrix(float phi, float the, float psi, float* resM)   //
{
  resM[0]=1.;
  resM[1]=psi+phi*the;
  resM[2]=phi*psi-the;
  resM[3]=-psi;
  resM[4]=1.-phi*the*psi;
  resM[5]=phi+the*psi;
  resM[6]=the;
  resM[7]=-phi;
  resM[8]=1.;
  return 1;
};

void FindEuler_new(double *resMO, double *phiO, double *theO, double *psiO)  //doesn't work for some reason
{
  if (fabs(resMO[8])>0.99999)
  { if (resMO[4]> 1.) resMO[4] = 1.;
	if (resMO[4]<-1.) resMO[4] = -1.;
	*psiO=acos(resMO[4]);
	*theO=*phiO=0.;
  } else
  {
	double s1,s2,s3,c1,c2,c3;
	s1 = 1./(1.-resMO[8]*resMO[8]);
	s2 = -2.*resMO[2]*resMO[5]*s1;
	c2 = (SQR(resMO[5])-SQR(resMO[2]))*s1;
	s3 = 2.*resMO[6]*resMO[7]*s1;
	c3 = (SQR(resMO[7])-SQR(resMO[6]))*s1;
	*psiO = 0.5*(c2>0?asin(s2):M_PI-asin(s2)) + (resMO[2]<0.?M_PI:0.);
	*phiO = 0.5*(c3>0?asin(s3):M_PI-asin(s3)) + (resMO[6]<0.?M_PI:0.);
	*theO = acos(resMO[8]);
	if (*psiO<0.) *psiO += 2.*M_PI;
	if (*phiO<0.) *phiO += 2.*M_PI;
  }
}

void FindEuler(double *resMO, double *phiO, double *theO, double *psiO)
{ double s1,s2,s3,c1,c2,c3;
  if (resMO[8]> 1.) resMO[8] = 1.;
  if (resMO[8]<-1.) resMO[8] = -1.;
  *theO=acos(resMO[8]);
  s1=sin(*theO);
  if (fabs(s1)<1e-3)
  { if (resMO[4]> 1.) resMO[4] = 1.;
	if (resMO[4]<-1.) resMO[4] = -1.;
	*psiO=acos(resMO[4]);
	*theO=*phiO=0.;
  } else
  { c2 = -resMO[5]/s1;
	s2 =  resMO[2]/s1;
	c3 =  resMO[7]/s1;
	s3 =  resMO[6]/s1;
	if (c2> 1.) c2 = 1.;
	if (c2<-1.) c2 = -1.;
	if (c3> 1.) c3 = 1.;
	if (c3<-1.) c3 = -1.;
	if (s2>0.) *psiO = acos(c2);
	else *psiO = 2.*M_PI - acos(c2);
	if (s3>0.) *phiO = acos(c3);
	else *phiO = 2.*M_PI - acos(c3);
  }
}


// / 1     0         0    \ /  cos(the) 0 sin(the) \ / cos(psi) -sin(psi) 0 \
// | 0 cos(phi) -sin(phi) | |     0     1     0    | | sin(psi)  cos(psi) 0 |
// \ 0 sin(phi)  cos(phi) / \ -sin(the) 0 cos(the) / \    0         0     1 /
// Rotation clockwize
int Rotate3ang(double phi,    double the,    double psi,
			   double AX,     double AY,     double AZ,
			   double* AoutX, double* AoutY, double* AoutZ)
{ double rotM1[9], rotM2[9], rotM3[9], resM[9], inpC[3], outC[3];
  for (int i=0; i<9; i++)
  { rotM1[i] = 0.;
	rotM2[i] = 0.;
	rotM3[i] = 0.;
  }
  rotM1[0] =  rotM2[4] = rotM3[8] = 1.;
  rotM1[4] =  rotM1[8] = cos(phi);
  rotM1[7] = sin(phi); rotM1[5] = -rotM1[7];
  rotM2[0] =  rotM2[8] = cos(the);
  rotM2[2] = sin(the); rotM2[6] = -rotM2[2];
  rotM3[0] =  rotM3[4] = cos(psi);
  rotM3[3] = sin(psi); rotM3[1] = -rotM3[3];
  MultMatrix(rotM1, 3, 3, rotM2, 3, 3, resM);
  MultMatrix(resM, 3, 3, rotM3, 3, 3, rotM1);
  inpC[0] = AX; inpC[1] = AY; inpC[2] = AZ;
  MultMatrix(rotM1, 3, 3, inpC, 3, 1, outC);
  *AoutX = outC[0]; *AoutY = outC[1]; *AoutZ = outC[2];
  return 1;
};

int Rotate3angMatrix(float aX, float aY, float aZ, float* resM)
{
  float rotM1[9], rotM2[9], rotM3[9], resMt[9];
  for (int i=0; i<9; i++)
  { rotM1[i] = 0.;
	rotM2[i] = 0.;
	rotM3[i] = 0.;
  }
  rotM1[0] =  rotM2[4] = rotM3[8] = 1.;
  rotM1[4] =  rotM1[8] = cos(aX);
  rotM1[7] = sin(aX); rotM1[5] = -rotM1[7];
  rotM2[0] =  rotM2[8] = cos(aY);
  rotM2[2] = sin(aY); rotM2[6] = -rotM2[2];
  rotM3[0] =  rotM3[4] = cos(aZ);
  rotM3[3] = sin(aZ); rotM3[1] = -rotM3[3];
  MultMatrixF(rotM1, 3, 3, rotM2, 3, 3, resMt);
  MultMatrixF(resMt, 3, 3, rotM3, 3, 3, resM);
  return 1;
};

int MirrorMatrix(bool mx, bool my, bool mz, float* resM)
{ resM[1] = resM[2] = resM[3] = resM[5] = resM[6] = resM[7] = 0;
  resM[0] = resM[4] = resM[8] = 1;
  if (mx) resM[0] = -1;
  if (my) resM[4] = -1;
  if (mz) resM[8] = -1;
  return 1;
}


//mat(1,0,0,0,cos(a),sin(a),0,-sin(a),cos(a),3,3)*mat(cos(b),0,-sin(b),0,1,0,sin(b),0,cos(b),3,3)*mat(cos(g),sin(g),0,-sin(g),cos(g),0,0,0,1,3,3)

//mat(cos(b)*cos(g), cos(b)*sin(g), -sin(b),
//   sin(a)*sin(b)*cos(g)-cos(a)*sin(g), sin(a)*sin(b)*sin(g)+cos(a)*cos(g), sin(a)*cos(b),
//   cos(a)*sin(b)*cos(g)+sin(a)*sin(g), cos(a)*sin(b)*sin(g)-sin(a)*cos(g), cos(a)*cos(b))

int RotCartMatrix(float alf, float bet, float gam, float* resM)
{ double ca,cb,cg,sa,sb,sg;
  ca=cos(alf);
  cb=cos(bet);
  cg=cos(gam);
  sa=sin(alf);
  sb=sin(bet);
  sg=sin(gam);
  resM[0]=cb*cg;
  resM[1]=cb*sg;
  resM[2]=-sb;
  resM[3]=sa*sb*cg-ca*sg;
  resM[4]=sa*sb*sg+ca*cg;
  resM[5]=sa*cb;
  resM[6]=ca*sb*cg+sa*sg;
  resM[7]=ca*sb*sg-sa*cg;
  resM[8]=ca*cb;
  return 1;
};

int RotAroundAxisMatrix(float the, float x, float y, float z, float* resM)
{ double ct = cos(the);
  double st = sin(the);

  resM[0]=ct+(1-ct)*x*x;
  resM[1]=(1-ct)*x*y-st*z;
  resM[2]=(1-ct)*x*z+st*y;
  resM[3]=(1-ct)*y*x+st*z;
  resM[4]=ct+(1-ct)*y*y;
  resM[5]=(1-ct)*y*z-st*x;
  resM[6]=(1-ct)*z*x-st*y;
  resM[7]=(1-ct)*z*y+st*x;
  resM[8]=ct+(1-ct)*z*z;
  return 1;
};



int RotateEuler3DVolume(float* vol1, float* vol2, int* dim, float fillWith,
                        float phi, float tet, float psi)
{     int size2 = dim[0]*dim[1]*dim[2];
      float cen2[3] = {0.5*(dim[0]-1),0.5*(dim[1]-1),0.5*(dim[2]-1)};
      for (size_t i=0; i<size2; i++) vol2[i] = fillWith;

      for (int iz=0; iz<dim[2]; iz++)
       for (int iy=0; iy<dim[1]; iy++)
        for (int ix=0; ix<dim[0]; ix++)
        { double cof[3];
          int coi[3] = {0,0,0};
          RotateEuler(phi, tet, psi, ix-cen2[0], iy-cen2[1], iz-cen2[2],
                      &cof[0], &cof[1], &cof[2]);
          bool outOf = false;
          for (int i=0; i<3; i++)
          { coi[i] = roundMy(cof[i]+cen2[i]);
            if (coi[i]<0 || coi[i]>=dim[i]) outOf = true;
          }
          if (outOf) continue;
          if (vol1[ix+dim[0]*(iy+dim[1]*iz)]<=fillWith) continue;
          vol2[coi[0]+dim[0]*(coi[1]+dim[1]*coi[2])] = vol1[ix+dim[0]*(iy+dim[1]*iz)];
        }
      return 1;
}


int RotateEuler3DVolumeInterp(float* vol1, float* vol2, int* dim, float fillWith,
                   float phi, float tet, float psi, float radiusF, float radiusE)
{ int size2 = dim[0]*dim[1]*dim[2];
  float cen2[3] = {0.5*(dim[0]-1),0.5*(dim[1]-1),0.5*(dim[2]-1)};
  for (size_t i=0; i<size2; i++) vol2[i] = fillWith;

  struct TPoints
  { float x;
    float y;
    float z;
    float I;
    TPoints* ne;
  };
  TPoints** AllPo = new TPoints*[size2];
  for (int i=0; i<size2; i++) AllPo[i] = NULL;

  TPoints** cPo;

  for (int iz=0; iz<dim[2]; iz++)
   for (int iy=0; iy<dim[1]; iy++)
    for (int ix=0; ix<dim[0]; ix++)
    { double cof[3];
      int coi[3] = {0,0,0};
      RotateEuler(phi, tet, psi, ix-cen2[0], iy-cen2[1], iz-cen2[2],
                      &cof[0], &cof[1], &cof[2]);
      bool outOf = false;
      for (int i=0; i<3; i++)
      { coi[i] = roundMy(cof[i]+cen2[i]);
        if (coi[i]<0 || coi[i]>=dim[i]) outOf = true;
      }
      if (outOf) continue;

//      int _ind = coi[2]+dim[2]*(coi[1]+dim[1]*coi[0]);
      int _ind = coi[0]+dim[0]*(coi[1]+dim[1]*coi[2]);
//-      if (vol1[_ind]<=fillWith) continue;
      float valI = vol1[ix+dim[0]*(iy+dim[1]*iz)];
//      float valI = vol1[iz+dim[2]*(iy+dim[1]*ix)];
      if (valI<=fillWith) continue;

      cPo = &AllPo[_ind];
      while (*cPo != NULL)
        cPo = &(*cPo)->ne;

      *cPo = new TPoints;
      (*cPo)->ne = NULL;
      (*cPo)->z = cof[2]+cen2[2];
      (*cPo)->y = cof[1]+cen2[1];
      (*cPo)->x = cof[0]+cen2[0];
      (*cPo)->I = valI;
//?      (*cPo)->I = vol1[ix+dim[0]*(iy+dim[1]*iz)];
//          vol2[coi[0]+dim[0]*(coi[1]+dim[1]*coi[2])] = vol1[ix+dim[0]*(iy+dim[1]*iz)];
    }

// actual interpolation

  /// radiusE - different radius for empty points!
  float radiusF2 = radiusF*radiusF;
  float radiusE2 = radiusE*radiusE;
  float _timeS = time(0);
//  int* vol3Dnum = new int[size2];
//  for (int i=0; i<size2; i++) vol3Dnum[i]=0;
#pragma omp parallel for
//private(_dist, wi, sumWi, _cI, _cIi, _cIr, cPo)
  for (int zi=0; zi<dim[2]; zi++)
  {
    double _cI, sumWi, wi, _dist;

#ifdef _OPENMP
	  printf("Thread %d executes zi=%d of %d at %d\n", omp_get_thread_num(),zi,dim[2],(int)(time(0)-_timeS));
#endif
    for (int yi=0; yi<dim[1]; yi++)
	   for (int xi=0; xi<dim[0]; xi++)
	  { 
      sumWi = 0.;
	    _cI = 0.;

	    bool bigRad = false;

      float radius;
      float radius2;
	    if (AllPo[xi+dim[0]*(yi+dim[1]*zi)]==NULL)
      { radius = radiusE;
        radius2 = radiusE2;
	    } else
      { radius = radiusF;
        radius2 = radiusF2;
      }

	    for (int xii=xi-radius; xii<=xi+radius; xii++)
	     for (int yii=yi-radius; yii<=yi+radius; yii++)
		    for (int zii=zi-radius; zii<=zi+radius; zii++)
		    { if (xii<0 || yii<0 || zii<0 || xii>=dim[0] || yii>=dim[1] || zii>=dim[2]) continue;
		      cPo = &AllPo[xii+dim[0]*(yii+dim[1]*zii)];
		      while (*cPo != NULL)
		      { _dist = (SQR((*cPo)->x-xi)+SQR((*cPo)->y-yi)+SQR((*cPo)->z-zi));
			      if (_dist > radius2)
			      { cPo = &(*cPo)->ne;
			        continue;
			      }
			      if (_dist<MinVal) // a point at this location found!
			      { sumWi = 1.;
			        _cI = (*cPo)->I;
			        xii=dim[0]; yii=dim[1]; zii=dim[2];
			        break; // must be up to sumWi=0
			      }
			      _dist = sqrt(_dist);

			      wi = SQR((radius-_dist)/(radius*_dist));
			      sumWi += wi;
			      _cI += wi * (*cPo)->I;
        		cPo = &(*cPo)->ne;
		      }
		    }
	    if (sumWi>MinVal) //???
        vol2[xi+dim[0]*(yi+dim[1]*zi)] = _cI/sumWi;
	    else
        vol2[xi+dim[0]*(yi+dim[1]*zi)] = fillWith;
	  }
  }

//  delete vol3Dnum;
  delete AllPo;
  return 1;
}

double AbsVec(double AX, double AY, double AZ)
{
  return sqrt(AX*AX+AY*AY+AZ*AZ);
};

double AbsVec2(double AX, double AY, double AZ)
{
  return AX*AX+AY*AY+AZ*AZ;
};

double MultScal(double AX, double AY, double AZ, double BX, double BY, double BZ)
{
  return AX*BX+AY*BY+AZ*BZ;
};

void AddVec(double AX, double AY, double AZ, double BX, double BY, double BZ,
             double* CX, double* CY, double* CZ)
{
  *CX = AX+BX;
  *CY = AY+BY;
  *CZ = AZ+BZ;
};

double AnglVec(double AX, double AY, double AZ, double BX, double BY, double BZ)
{
  double mmm = MultScal(AX,AY,AZ,BX,BY,BZ)/(AbsVec(AX,AY,AZ)*AbsVec(BX,BY,BZ));
  if (mmm >= 1) return 0.;
  else return acos(mmm);
};

void MultVect(double AX, double AY, double AZ, double BX, double BY, double BZ,
			 double* CX, double* CY, double* CZ)
{
  *CX = AY*BZ-AZ*BY;
  *CY = AZ*BX-AX*BZ;
  *CZ = AX*BY-AY*BX;
};

void ProectParalVect(double AX, double AY, double AZ, double BX, double BY, double BZ,
             double* CX, double* CY, double* CZ)
{       // �������� ������� ����������� �������
  double absB = AbsVec(BX,BY,BZ);
  double scal;
  if (absB < 1e-10) scal = 0;
  else scal = MultScal(AX,AY,AZ,BX,BY,BZ)/(absB*absB);
  *CX = scal*BX;
  *CY = scal*BY;
  *CZ = scal*BZ;
};

void ProectPerpVect(double AX, double AY, double AZ, double BX, double BY, double BZ,
             double* CX, double* CY, double* CZ)
{       // �������� ������� �� ��������� ���������������� �������
  double tCX, tCY, tCZ;
  double absB = AbsVec(BX,BY,BZ);
  if (absB < 1e-10)
  { *CX = 0.; *CY = 0.; *CZ = 0.;}
  else
  { MultVect(AX,AY,AZ,BX,BY,BZ,&tCX,&tCY,&tCZ);
	tCX /= (absB*absB);
    tCY /= (absB*absB);
    tCZ /= (absB*absB);
    MultVect(BX,BY,BZ,tCX,tCY,tCZ,CX,CY,CZ);
  };
};

void MakeOrt(double* AX, double* AY, double* AZ)
{
  double absB = AbsVec(*AX,*AY,*AZ);
  *AX /= (absB);
  *AY /= (absB);
  *AZ /= (absB);
}

// ------------------------- END OF VECTORS -------------------------------

// ----------------------- Background and noise filters-------------------------

bool MedianFilter1D(float* inpAr, int numX, int radX, float badVal)
{
  float* arr1 = new float[(2*radX+1)];
  float* arrOld = new float[numX];
  for (int i=0; i<numX; i++)
    arrOld[i] = inpAr[i];

  for (int xi=0; xi<numX; xi++)
    { int nuel = 0;
      for (int bxi=-radX; bxi<=radX; bxi++)
          if (xi+bxi>=0 && xi+bxi<numX)
            if (arrOld[xi+bxi]>badVal+1)
            { arr1[nuel] = arrOld[xi+bxi];
              nuel++;
            }
      if (nuel>0)
        inpAr[xi] = medianCutoff(arr1, nuel, 0.5);
      else inpAr[xi] = badVal;
    }
  delete[] arr1;
  delete[] arrOld;
  return true;
}
                                                                                    // to delete after: RadNumComp, RadialArray, shellElem2

bool AverageFilter1D(float* inpAr, int numX, int radX, float badVal)
{
  float* arrOld = new float[numX];
  for (int i=0; i<numX; i++)
    arrOld[i] = inpAr[i];

  for (int xi=0; xi<numX; xi++)
    { int nuel = 0;
      float sum = 0.;
      for (int bxi=-radX; bxi<=radX; bxi++)
          if (xi+bxi>=0 && xi+bxi<numX)
            if (arrOld[xi+bxi]>badVal+1)
            { sum += arrOld[xi+bxi];
              nuel++;
            }
      if (nuel>0)
        inpAr[xi] = sum/(float)nuel;
      else inpAr[xi] = badVal;
    }
  delete[] arrOld;
  return true;
}
                                                                                    // to delete after: RadNumComp, RadialArray, shellElem2
bool MergeSeveralShells(size_t numEl, int maxShellsMerged, int minNumInShell, int* pix_r,
                        int* maxRad, int* numPerShell1, int*** shellElem1)
{
  if (maxShellsMerged <= 1) return false;
    int** shellElLoc = *shellElem1;

    int* numPerShell2 = new int[*maxRad];      // new one - now use this
    for (int i=0; i<*maxRad; i++)
      numPerShell2[i] = 0;
    int** shellElem2 = new int*[*maxRad];
    int numCombined = 1;
    int maxRadNew = 0;
    for (int i=0; i<*maxRad; i+=numCombined)
    { //numCombined = 1;  //?
      int curInShell = 0;//numPerShell1[i];
      int j;
      for (j=0; j<maxShellsMerged; j++)
      { if (i+j>=*maxRad)
          break;
        if (curInShell<minNumInShell)
          curInShell += numPerShell1[i+j];
        else
          break;
      }
      numCombined = j; // doesn't look right
      if (curInShell<1) continue;

      numPerShell2[maxRadNew] = curInShell;
      shellElem2[maxRadNew] = new int[numPerShell2[maxRadNew]];
      int kN=0;
      for (j=0; j<numCombined; j++)
        for (int k=0; k<numPerShell1[i+j]; k++)
        { shellElem2[maxRadNew][kN] = shellElLoc[i+j][k];
//        { shellElem2[maxRadNew][kN] = *shellElem1[i+j][k];
          kN++;
          if (kN>numPerShell2[maxRadNew])
          { printf("ups\n");
            break;
          }
        }
    // here kN must be = numPerShell[maxRadNew] = curInShell
      maxRadNew++;
    }
    for (int i=0; i<*maxRad; i++)
      if (shellElLoc[i] != NULL)
        delete shellElLoc[i];
    delete shellElLoc;

    //output
    for (int i=0; i<maxRadNew; i++) numPerShell1[i] = numPerShell2[i];
    delete numPerShell2;
    *maxRad = maxRadNew;
    *shellElem1 = shellElem2;

    // change pix_r to new bins
    if (pix_r != NULL)
    { for (int i=0; i<numEl; i++) pix_r[i] = -1;
      for (int i=0; i<*maxRad; i++)
        for (int k=0; k<numPerShell1[i]; k++)
          pix_r[shellElem2[i][k]] = i;
    }

  return true;
}

// building array with radial coordinates for each pixel
bool BuildRadialArray(size_t numEl, float *cox, float* coy, float istep,
					  int* pix_r, int* maxRad, float* pixelsR)
{
  // determine number of components in radial array
  *maxRad = 0;
  for (size_t i=0; i<numEl; i++)
  { float distR = sqrt(SQR(cox[i])+SQR(coy[i]))*istep;
	if (pixelsR != NULL) pixelsR[i] = distR;
	size_t _ind = roundMy(distR);
	if (_ind > *maxRad) *maxRad = _ind;
  }
  *maxRad += 1;

  // forming the array
  for (int i=0; i<numEl; i++)
  { int _ind = roundMy(sqrt(SQR(cox[i])+SQR(coy[i]))*istep);
  //should never happen
	if (_ind>=*maxRad)
	  {printf("In radial averaging components are too far!\n"); continue;}
	pix_r[i] = _ind;
  }
  return true;
}

// building 2D array with radial coordinates for each pixel
//maskAr, radCoord could be NULL if not needed                   // pix_r - maybe leave it just internaly
bool BuildRadial2DArray(size_t numEl, char* maskAr, int maxShellsMerged,
						float *cox, float* coy, float istep, int* maxRad,
						int** radCoord, int** numPerShell, int*** shellElem, float* pixelsR)
{ //int maxShellsMerged = 5;  // must be >=1
  int minNumInShell = 50*roundMy(0.5*numEl/(*maxRad)); // also must be >=1                     //??? probably should be a parameter

  //checking maxRad
/*  size_t numEl = (size_t)numX*(size_t)numY;
  int maxRad = 0;
  for (size_t i=0; i<numEl; i++)
	if (pix_r[i]>maxRad) maxRad = pix_r[i];
  if (maxRad<1) return 0;
  maxRad += 1;
*/
  // creating the array
  int* pix_r = new int[numEl];
  for (int i=0; i<numEl; i++)
	pix_r[i] = -1;
  BuildRadialArray(numEl, cox, coy, istep, pix_r, maxRad, pixelsR);

  if (maxShellsMerged<1 || minNumInShell<1)
  { printf("Error, some input values are wrong\n");
    delete pix_r;
    return false;
  }

  // arrays with number of elements for each shell
  int* numPerShell1 = new int[*maxRad];
  for (int i=0; i<*maxRad; i++)
    numPerShell1[i] = 0;
  for (int i=0; i<numEl; i++)
    numPerShell1[pix_r[i]]++;

  // 2d array for storing indexes for each shell
  int** shellElem1 = new int*[*maxRad];
  for (int i=0; i<*maxRad; i++)
  { shellElem1[i] = new int[numPerShell1[i]];
    for (int j=0; j<numPerShell1[i]; j++)
      shellElem1[i][j] = 0;
  }

  // again clear array
  for (int i=0; i<*maxRad; i++)
	numPerShell1[i] = 0;
  // actual forming of the 2D array
  for (size_t i=0; i<numEl; i++)
  { if (maskAr != NULL) if (maskAr[i]<MinVal) continue;  // taking into account mask if provided
    shellElem1[pix_r[i]][numPerShell1[pix_r[i]]] = i;
    numPerShell1[pix_r[i]]++;
  }

  // here I can merge several shells if there are not enough values
//already there  if (maxShellsMerged > 1)
    MergeSeveralShells(numEl, maxShellsMerged,minNumInShell, pix_r, maxRad,
					   numPerShell1, &shellElem1);

  *shellElem = shellElem1;
  *numPerShell = numPerShell1;
  if (radCoord != NULL) *radCoord = pix_r;
  else delete pix_r;

  return true;

}

// building array with radial coordinates for each pixel
int BuildRadialArray3D(size_t numEl, float *cox, float* coy, float* coz,
						float istep, int* pix_r, float* pixelsR)
{ int maxRad = 0;
  if (pixelsR==NULL && cox==NULL && coy==NULL && coz==NULL)
    return 0;
  if (pixelsR!=NULL && cox==NULL && coy==NULL && coz==NULL)
    for (size_t i=0; i<numEl; i++)
    { float distR = pixelsR[i];
      size_t _ind = roundMy(distR);
      if (_ind > maxRad) maxRad = _ind;   // determine number of components in radial array
      pix_r[i] = _ind;
    }
  else
    for (size_t i=0; i<numEl; i++)
	  { float distR = sqrt(SQR(cox[i])+SQR(coy[i])+(coz!=NULL?SQR(coz[i]):0))*istep;
	    if (pixelsR != NULL) pixelsR[i] = distR;
	    size_t _ind = roundMy(distR);
	    if (_ind > maxRad) maxRad = _ind;   // determine number of components in radial array
	    pix_r[i] = _ind;
	  }
  maxRad += 1;
  return maxRad;
}

// building 2D array with radial coordinates for each pixel
//maskAr, radCoord could be NULL if not needed                   // pix_r - maybe leave it just internaly
int BuildRadial2DArray3D(size_t dims[3], char* maskAr, int maxShellsMerged,
						 float *cox, float* coy, float* coz, float istep,
						 int** radCoord, int** numPerShell,
						 int*** shellElem, float* pixelsR)
{ //int maxShellsMerged = 5;  // must be >=1
  //int minNumInShell = 50*roundMy(0.5*numEl/(*maxRad)); // also must be >=1                     //??? probably should be a parameter
  //checking maxRad
/*  size_t numEl = (size_t)numX*(size_t)numY;
  int maxRad = 0;
  for (size_t i=0; i<numEl; i++)
    if (pix_r[i]>maxRad) maxRad = pix_r[i];
  if (maxRad<1) return 0;
  maxRad += 1;
*/
  for (int i=0; i<3; i++) if (dims[i]<1) dims[i] = 1;
  size_t numEl = dims[0]*dims[1]*dims[2];
  int maxRad = 0;

  // if simple 2D-3D array, create coordinates

  int* pix_r = new int[numEl];
  if (cox==NULL && coy==NULL && coz==NULL)
  { float* pixR = pixelsR;
	if (pixR==NULL)
    { pixR = new float[numEl];
      float cen[3] = {0,0,0};
      for (int i=0; i<3; i++)
        cen[i] = 0.5*(dims[i]-1.);
      size_t el = 0;
      for (size_t zi=0; zi<dims[2]; zi++)
       for (size_t yi=0; yi<dims[1]; yi++)
        for (size_t xi=0; xi<dims[0]; xi++)
        { pixR[el] = sqrt(SQR(cen[0]-xi)+SQR(cen[1]-yi)+SQR(cen[2]-zi));
          el++;
        }

    }
    maxRad = BuildRadialArray3D(numEl, cox, coy, coz, istep, pix_r, pixR);
    if (pixelsR==NULL) delete pixR;
  } else
  // creating the array
    maxRad = BuildRadialArray3D(numEl, cox, coy, coz, istep, pix_r, pixelsR);

  int minNumInShell = 50*roundMy(0.5*numEl/((float)maxRad)); // also must be >=1                     //??? probably should be a parameter

  if (maxShellsMerged<1 || minNumInShell<1)
  { printf("Error, some input values are wrong\n");
    delete pix_r;
	return -1;
  }

  // arrays with number of elements for each shell
  int* numPerShell1 = new int[maxRad];
  for (int i=0; i<maxRad; i++)
    numPerShell1[i] = 0;
  for (int i=0; i<numEl; i++)
	numPerShell1[pix_r[i]]++;

  // 2d array for storing indexes for each shell
  int** shellElem1 = new int*[maxRad];
  for (int i=0; i<maxRad; i++)
  { shellElem1[i] = new int[numPerShell1[i]];
    for (int j=0; j<numPerShell1[i]; j++)
      shellElem1[i][j] = 0;
  }

  // again clear array
//  int* numPerShell4 = new int[*maxRad];

  for (int i=0; i<maxRad; i++)
    numPerShell1[i] = 0;
  // actual forming of the 2D array
  for (size_t i=0; i<numEl; i++)
  { if (maskAr != NULL) if (maskAr[i]<MinVal) continue;  // taking into account mask if provided
    shellElem1[pix_r[i]][numPerShell1[pix_r[i]]] = i;
    numPerShell1[pix_r[i]]++;
  }

  // here I can merge several shells if there are not enough values
//already there  if (maxShellsMerged > 1)
  if (maxShellsMerged>1) MergeSeveralShells(numEl, maxShellsMerged,minNumInShell, pix_r, &maxRad,
                       numPerShell1, &shellElem1);

  *shellElem = shellElem1;
  *numPerShell = numPerShell1;
  if (radCoord != NULL) *radCoord = pix_r;
  else delete pix_r;

  return maxRad;

}

int BuildSimpleRadialArray3D(size_t* dims, int** radialAr)
{
  // determine number of components in radial array
//  int cen[3] = {0.5*(dims[0]-1.)+MinVal, 0.5*(dims[1]-1.)+MinVal, 0.5*(dims[2]-1.)+MinVal};
  float cen[3];
  for (int i=0; i<3; i++) 
    cen[i] = 0.5*(dims[i]-1.);
  int maxRad = roundMy(0.5*sqrt(SQR(dims[0])+SQR(dims[1])+SQR(dims[2])))+2;
  size_t numEl = dims[0]*dims[1]*dims[2];
  int* pix_r = new int[numEl];
  *radialAr = pix_r;

  // forming the array
  for (int zi=0; zi<dims[2]; zi++)
   for (int yi=0; yi<dims[1]; yi++)
    for (int xi=0; xi<dims[0]; xi++)
    { int _ind = roundCpp(sqrt(SQR(xi-cen[0])+SQR(yi-cen[1])+SQR(zi-cen[2])));
      //should never happen
      if (_ind>=maxRad)
        printf("In radial averaging components are too far!\n");
      else
        pix_r[xi+dims[0]*(yi+dims[1]*zi)] = _ind;
    }
  return maxRad;
}

float LoopMedianCutoffBG(float* data, long pix_nn, float bstpReg, float cutoff)
{
//! pix_r - array of radii of each pixel
//! hitfinderMinSNR - ?  it is =8 in cheetah.ini
//! ADCthresh -?  Some minimal threshold
  float hitfinderMinSNR = 8;
  float ADCthresh = bstpReg; //???

  const int loopsNum = 5; // number of loops to make

  // Allocate and zero arrays
  float rthreshold = 1e30;
  float lowthreshold = bstpReg+1;
  float roffset;
  float rsigma;
  long rcount;

  // Compute sigma and average of data values at each radius
  // From this, compute the ADC threshold to be applied at each radius
  // Iterate a few times to reduce the effect of positive outliers (ie: peaks)
  for(long counter=0; counter<loopsNum-1; counter++)
  {
    roffset = 0;
    rsigma = 0;
    rcount = 0;

    for (long i=0;i<pix_nn;i++)
		if(data[i] < rthreshold && data[i] > lowthreshold)
        { roffset += data[i];
          rsigma += (data[i]*data[i]);
          rcount += 1;
        }

      if (rcount == 0)
      { roffset = 0;
        rsigma = 0;
        rthreshold = 1e9;
        lowthreshold = -1e9;
      }
      else
      { float thisoffset = roffset/(float)rcount;
        roffset = thisoffset;
        float someval = rsigma/(float)rcount - thisoffset*thisoffset;
        if (someval>MinVal) rsigma = sqrt(someval);
        else rsigma = 0;
//        rsigma = sqrt(rsigma/rcount - thisoffset*thisoffset);
        rthreshold = roffset + hitfinderMinSNR*rsigma;
		lowthreshold = roffset - hitfinderMinSNR*rsigma;
        if(rthreshold < ADCthresh) rthreshold = ADCthresh;
      }
  } // main loop

  // last run with median
  // need arrays for medians (each radii):
  roffset = bstpReg;
  if (rcount>0)
  {
    float* radmed = new float[rcount];
    rsigma = 0;
    rcount = 0;

    for (long i=0;i<pix_nn;i++)
      if(data[i] < rthreshold && data[i] > lowthreshold)
	  { radmed[rcount] = data[i];
        rcount++;
      }

    if (rcount>0)
      roffset = medianCutoff(radmed, rcount, cutoff);
    delete radmed;
  }
  return roffset;
}

// it looks like this function just reject outliers
long LoopAverageBG(float* data, int pix_nn, float bstpReg, unsigned int loopsNum, float hitfinderMinSNR, int cleandata)
{
//! hitfinderMinSNR - ?  it is =8 in cheetah.ini
//  float hitfinderMinSNR = 8;
//  const int loopsNum = 5; // number of loops to make
  if (pix_nn<1) return 0;
  // Allocate and zero arrays
  float rthreshold = 1e30;
  float lowthreshold = bstpReg+1;
  float roffset;
  float rsigma;
  int rcount;

  // Compute sigma and average of data values at each radius
  // From this, compute the ADC threshold to be applied at each radius
  // Iterate a few times to reduce the effect of positive outliers (ie: peaks)
  for(unsigned int counter=0; counter<loopsNum; counter++)
  {
    roffset = 0;
    rsigma = 0;
    rcount = 0;

    for (int i=0; i<pix_nn; i++)
      if (data[i]>lowthreshold && data[i]<rthreshold)
	  { roffset += data[i];
        rsigma += (data[i]*data[i]);
        rcount += 1;
      }

    if (rcount == 0) break;

    float thisoffset = roffset/(float)rcount;
    roffset = thisoffset;
    float someval = rsigma/(float)rcount - thisoffset*thisoffset;
    if (someval>MinVal) rsigma = sqrt(someval);
    else rsigma = 0;
//    rsigma = sqrt(rsigma/rcount - thisoffset*thisoffset);
    rthreshold = roffset + hitfinderMinSNR*rsigma;
    lowthreshold = roffset - hitfinderMinSNR*rsigma;
  } // main loop
  if (rcount>0)
  { rcount = 0;
    for (int i=0; i<pix_nn; i++)
      if(data[i] < rthreshold && data[i] > lowthreshold)
	  { data[rcount] = data[i];                            // throwing away all bad values
        rcount++;
      }
//       else if (cleandata==1)                             // no idea what is it for
//        data[i] = bstpReg;
  };

  return rcount;
}

// if radialcurve==NULL, just subtract, if not - no subtraction, just export radial curve
// also if maskAr==NULL, then all data<bstpReg is considered as mask
bool SubtractRadial(float* data, size_t numEl, char* maskAr, int* numPerShell, int** shellElem,
					int maxRad, float bstpReg, float** radialcurve, int mode, int curveSmooth)
{
  float* radAver;
  if (radialcurve==NULL) radAver = new float[maxRad];
  else radAver = *radialcurve;

#pragma omp parallel for schedule(dynamic, 1) shared(radAver)
  for (int i=0; i<maxRad; i++)
  { radAver[i] = bstpReg;
    if (numPerShell[i]<1) continue;
    float* dataSh = new float[numPerShell[i]];
    int numSh = 0;
    for (int k=0; k<numPerShell[i]; k++)
    { bool toskip = false;
      if (maskAr != NULL)
      { if (maskAr[shellElem[i][k]]<=0)
          toskip = true;
      } else
        if (data[shellElem[i][k]]<=bstpReg+1)
          toskip = true;
      if (toskip) continue;
      dataSh[numSh] = data[shellElem[i][k]];
      numSh++;
	}
    if (numSh>=1)
    { radAver[i] = 0;
      if (mode==0) radAver[i] = medianCutoff(dataSh, numSh, 0.5);     //median
      else if (mode==1) radAver[i] = LoopMedianCutoffBG(dataSh, numSh, bstpReg, 0.5); //Anton's 5 loop
      else   // simple average
      { float ave = 0.;
        for (int ii=0; ii<numSh; ii++)
          ave += dataSh[ii];
        radAver[i] = ave/(float)numSh;
      }
    }
	delete dataSh;
  }

  // Here we exclude radial bins without data
  float* smoothRad = new float[maxRad];
  int* nSmoothRad = new int[maxRad];
  int radNum = 0;
  for (int i=0; i<maxRad; i++)                  //  what to do here? I can't exclude radial bins
//    if (radAver[i] > bstpReg+1)                 // Shell I fill with previous instead?
      { smoothRad[radNum] = radAver[i];
        nSmoothRad[radNum] = i;
        radNum++;
      };

//  SaveFloat1DtoTextFile("radial.txt", smoothRad, maxRad);
  // and now we smooth radial curve
  if (curveSmooth>0 && radNum>0)
    AverageFilter1D(smoothRad, radNum, curveSmooth, bstpReg);

//  SaveFloat1DtoTextFile("radialSm.txt", smoothRad, maxRad);

  // and get smoothed values back to the radAver array
//  for (int i=0; i<maxRad; i++) radAver[i] = 0;
  for (int i=0; i<radNum; i++)
	radAver[nSmoothRad[i]] = smoothRad[i];
  delete smoothRad;
  delete nSmoothRad;

  // subtracting or exporting radcurve
  if (radialcurve == NULL)
//    *radialcurve = radAver;
//  else
  { for (int i=0; i<maxRad; i++)
      for (int k=0; k<numPerShell[i]; k++)
      { bool toskip = false;
        if (maskAr != NULL)
		{ if (maskAr[shellElem[i][k]]<=0)
            toskip = true;
        } else
          if (data[shellElem[i][k]]<=bstpReg+1)
			toskip = true;
        if (toskip) continue;
        data[shellElem[i][k]] -= radAver[i];
//deb        data[shellElem[i][k]] = radAver[i];
      }
    delete radAver;
  }

  return true;
}

 //  Create a buffer for image data so we don't nuke the main image by mistake
void SubtractBgLoop(float* data, long pix_nn, int *pix_r, float hitfinderMinSNR, float ADCthresh, float bstpReg, float* radialcurve)
{
//! pix_r - array of radii of each pixel
//! hitfinderMinSNR - ?  it is =8 in cheetah.ini
//! ADCthresh -?  Some minimal threshold

  const int loopsNum = 5; // number of loops to make

  // Apply mask (multiply data by 0 to ignore regions - this makes data below threshold for peak finding)
//mask  for(long i=0;i<pix_nn;i++) temp[i] *= mask[i];

  // Determine noise and offset as a funciton of radius
  float fminr, fmaxr;
  long lminr, lmaxr;
  fminr = 1e9;
  fmaxr = -1e9;

  // Figure out radius bounds
  for(long i=0;i<pix_nn;i++){
    if (pix_r[i] > fmaxr)
	  fmaxr = pix_r[i];
    if (pix_r[i] < fminr)
      fminr = pix_r[i];
  }
  lmaxr = (int)ceil(fmaxr)+1;
  lminr = 0;

  // Allocate and zero arrays
  float* rsigma = (float*) calloc(lmaxr, sizeof(float));
  float* roffset = (float*) calloc(lmaxr, sizeof(float));
  long* rcount = (long*) calloc(lmaxr, sizeof(long));
  float* lowthreshold = (float*) calloc(lmaxr, sizeof(float));
  float* rthreshold = (float*) calloc(lmaxr, sizeof(float));

  for(long i=0; i<lmaxr; i++)
  { rthreshold[i] = 1e9;
    lowthreshold[i] = -1e9;
  }

  // Compute sigma and average of data values at each radius
  // From this, compute the ADC threshold to be applied at each radius
  // Iterate a few times to reduce the effect of positive outliers (ie: peaks)
  for(long counter=0; counter<loopsNum; counter++) {
    for(long i=0; i<lmaxr; i++) {
	  roffset[i] = 0;
      rsigma[i] = 0;
      rcount[i] = 0;
    }
    for(long i=0;i<pix_nn;i++){
//mask      if(mask[i] != 0)
      if (data[i]>bstpReg+1)
      { long thisr = pix_r[i];
        if(data[i] < rthreshold[thisr] && data[i] > lowthreshold[thisr]) {
          roffset[thisr] += data[i];
          rsigma[thisr] += (data[i]*data[i]);
          rcount[thisr] += 1;
		}
      }
    }
    for(long i=0; i<lmaxr; i++) {
      if(rcount[i] == 0) {
        roffset[i] = 0;
        rsigma[i] = 0;
        rthreshold[i] = 1e9;
        lowthreshold[i] = -1e9;
      }
      else {
        float thisoffset = roffset[i]/(float)rcount[i];        //!!!!! this is average. I need array for each i at last loop
        roffset[i] = thisoffset;
        float someval = rsigma[i]/(float)rcount[i] - thisoffset*thisoffset;
        if (someval>MinVal) rsigma[i] = sqrt(someval);
        else rsigma[i] = 0;
//        rsigma[i] = sqrt(rsigma[i]/rcount[i] - thisoffset*thisoffset);;
        rthreshold[i] = roffset[i] + hitfinderMinSNR*rsigma[i];
        lowthreshold[i] = roffset[i] - hitfinderMinSNR*rsigma[i];
        if(rthreshold[i] < ADCthresh) rthreshold[i] = ADCthresh;
         //rthreshold[i] = ADCthresh;  // For testing
      }
    }
  }

//  FILE* afile = fopen("antaver.txt","wt");
//  for(long i=0; i<lmaxr; i++)
//    fprintf(afile,"%0.2f\n",roffset[i]);
//  fclose(afile);

  if (radialcurve != NULL)
    for(long i=0;i<lmaxr;i++)
      radialcurve[i] = roffset[i];
  else
    for(long i=0;i<pix_nn;i++)
      if (data[i]>bstpReg+1)
		data[i] -= roffset[pix_r[i]];
  free (rsigma);
  free (roffset);
  free (rcount);
  free (rthreshold);
  free (lowthreshold);

  return;
}

float SubtractBgLoopMedianCutoff(float* data, long pix_nn, int *pix_r, float hitfinderMinSNR, float ADCthresh, float bstpReg, float* radialcurve, float cutoff, int curveSmooth)
{
//! pix_r - array of radii of each pixel
//! hitfinderMinSNR - ?  it is =8 in cheetah.ini
//! ADCthresh -?  Some minimal threshold

  const int loopsNum = 5; // number of loops to make
//  const int curveSmooth = 10;

  // Apply mask (multiply data by 0 to ignore regions - this makes data below threshold for peak finding)
//mask  for(long i=0;i<pix_nn;i++) temp[i] *= mask[i];

  // Determine noise and offset as a funciton of radius
  float fminr, fmaxr;
  long lminr, lmaxr;
  fminr = 1e9;
  fmaxr = -1e9;

  // Figure out radius bounds
  for(long i=0;i<pix_nn;i++){
    if (pix_r[i] > fmaxr)
      fmaxr = pix_r[i];
    if (pix_r[i] < fminr)
      fminr = pix_r[i];
  }
  lmaxr = (int)ceil(fmaxr)+1;
  lminr = 0;
//  printf("max rad = %d\n",lmaxr);

  // Allocate and zero arrays
//  float *rsigma = (float*) calloc(lmaxr, sizeof(float));
//  float *roffset = (float*) calloc(lmaxr, sizeof(float));
//  long *rcount = (long*) calloc(lmaxr, sizeof(long));
//  float *rthreshold = (float*) calloc(lmaxr, sizeof(float));
//  float *lowthreshold = (float*) calloc(lmaxr, sizeof(float));

  float *rsigma = new float[lmaxr];
  float *roffset = new float[lmaxr];
  long *rcount = new long[lmaxr];
  float *rthreshold = new float[lmaxr];
  float *lowthreshold = new float[lmaxr];

  for(long i=0; i<lmaxr; i++)
  { rthreshold[i] = 1e9;
    lowthreshold[i] = -1e9;
  }

  // Compute sigma and average of data values at each radius
  // From this, compute the ADC threshold to be applied at each radius
  // Iterate a few times to reduce the effect of positive outliers (ie: peaks)
  for(long counter=0; counter<loopsNum-1; counter++)
  {
    for(long i=0; i<lmaxr; i++)
    { roffset[i] = 0;
      rsigma[i] = 0;
      rcount[i] = 0;
    }

    for (long k=0; k<pix_nn; k++)
//mask      if(mask[i] != 0)
      if (data[k] > bstpReg+1.)
      { long thisr = pix_r[k];
		if(data[k] < rthreshold[thisr] && data[k] > lowthreshold[thisr])
        { roffset[thisr] += data[k];
          rsigma[thisr] += (data[k]*data[k]);
          rcount[thisr] += 1;
        }
      }

    for(long i=0; i<lmaxr; i++)
    { if (rcount[i] == 0)
      { roffset[i] = 0;
        rsigma[i] = 0;
        rthreshold[i] = 1e9;
        lowthreshold[i] = -1e9;
      }
      else
      { float thisoffset = roffset[i]/(float)rcount[i];        //!!!!! this is average. I need array for each i at last loop
        roffset[i] = thisoffset;
        float someval = rsigma[i]/(float)rcount[i] - thisoffset*thisoffset;
        if (someval>MinVal) rsigma[i] = sqrt(someval);
        else rsigma[i] = 0;
        rthreshold[i] = roffset[i] + hitfinderMinSNR*rsigma[i];
        lowthreshold[i] = roffset[i] - hitfinderMinSNR*rsigma[i];
        if(rthreshold[i] < ADCthresh) rthreshold[i] = ADCthresh;
         //rthreshold[i] = ADCthresh;  // For testing
	  }
    }
  } // main loop

  // last run with median
  // need arrays for medians (each radii):
  long *rcountL = (long*) calloc(lmaxr, sizeof(long));
//  float** radmed = new float*[pix_nn];
  float **radmed = (float**) calloc(lmaxr, sizeof(float*));

  for(long i=0; i<lmaxr; i++)
  { roffset[i] = 0;
	rsigma[i] = 0;
    rcountL[i] = rcount[i];
    rcount[i] = 0;
    radmed[i] = (float*) calloc(rcountL[i], sizeof(float));
  }

  for (long i=0;i<pix_nn;i++)
    if (data[i]>bstpReg+1)
    { long thisr = pix_r[i];
      if(data[i] < rthreshold[thisr] && data[i] > lowthreshold[thisr] && rcount[thisr] < rcountL[thisr])
      { radmed[thisr][rcount[thisr]] = data[i];
        rcount[thisr] += 1;
      }
    }

//  FILE *stream;
//  stream = fopen("radCount1.txt", "wt");
//  for (long i=0; i<lmaxr; i++) fprintf(stream,"%d\n",rcount[i]);
//  fclose(stream);
//  printf("saved\n");

#pragma omp parallel for schedule(dynamic, 1)
  for(long i=0; i<lmaxr; i++)
  { if (rcount[i] == 0)
	  roffset[i] = bstpReg;
    else
      roffset[i] = medianCutoff(radmed[i], rcount[i], cutoff);
  }

  for(long i=0; i<lmaxr; i++)
    free (radmed[i]);
  free (rcountL);
  free (radmed);

//  FILE* afile = fopen("antaver.txt","wt");
//  for(long i=0; i<lmaxr; i++)
//    fprintf(afile,"%0.2f\n",roffset[i]);
//  fclose(afile);

  float averBG = average(roffset+lmaxr/4,lmaxr/2, bstpReg);
//  printf("BG: %8.3f\n", averBG);
  //hear let's smooth the radial curve
//  SaveFloat1DtoTextFile("mean_before.txt", roffset, lmaxr);
  if (curveSmooth>0) MedianFilter1D(roffset, lmaxr, curveSmooth, bstpReg);
//  SaveFloat1DtoTextFile("mean_after05.txt", roffset, lmaxr);

  for(long i=0;i<pix_nn;i++)
    if (data[i]>bstpReg+1 && roffset[pix_r[i]]>bstpReg+1)
      data[i] -= roffset[pix_r[i]];
  if (radialcurve != NULL)
    for(long i=0;i<lmaxr;i++)
      radialcurve[i] = (roffset[i]>bstpReg?roffset[i]:0);
//  free (rsigma);
//  free (roffset);
//  free (rcount);
//  free (rthreshold);
//  free (lowthreshold);

  delete rsigma;
  delete roffset;
  delete rcount;
  delete rthreshold;
  delete lowthreshold;

  return averBG;
}

void MaxBg(float* data, long pix_nn, int *pix_r, float hitfinderMinSNR, float ADCthresh, float bstpReg, float* radialcurve)
{
//! pix_r - array of radii of each pixel
//! hitfinderMinSNR - ?  it is =8 in cheetah.ini
//! ADCthresh -?  Some minimal threshold
  if (radialcurve==NULL) return;
  // Determine noise and offset as a funciton of radius
  float fminr, fmaxr;
  long lminr, lmaxr;
  fminr = 1e9;
  fmaxr = -1e9;

  // Figure out radius bounds
  for(long i=0;i<pix_nn;i++){
    if (pix_r[i] > fmaxr)
      fmaxr = pix_r[i];
    if (pix_r[i] < fminr)
      fminr = pix_r[i];
  }
  lmaxr = (int)ceil(fmaxr)+1;
  lminr = 0;

  for(long i=0; i<lmaxr; i++)
  { radialcurve[i] = -1e9;
  }

    for(long i=0;i<pix_nn;i++)
      if (data[i]>bstpReg+1)
      { long thisr = pix_r[i];
		if(data[i] > radialcurve[thisr])
		  radialcurve[thisr] = data[i];
	  }

  return;
}


float Correlation(float* vol3D1, float* vol3D2, float badVal, size_t numel)
{
  double _sum = 0;
  double _sum1 = 0;
  double _sum2 = 0;

  for (size_t i = 0; i < numel; i++)
	if (vol3D1[i] > badVal && vol3D2[i] > badVal)
	{ _sum  += vol3D1[i]*vol3D2[i];
	  _sum1 += vol3D1[i]*vol3D1[i];
	  _sum2 += vol3D2[i]*vol3D2[i];
	}

  if (_sum1*_sum2>MinVal)
	_sum /= sqrt(_sum1*_sum2);
  return (float)_sum;
}

float CorrelationRadial(float* vol3D1, float* vol3D2, size_t numel, int* radAr, int radSize, float badVal, float* radCor)
{
  double _sum = 0;
  double _sum1 = 0;
  double _sum2 = 0;
  double* radSum1 = new double[radSize];
  double* radSum2 = new double[radSize];
//  int* radNum = new int[radSize];
  for (int i=0; i<radSize; i++)
  { radCor[i] = 0;
    radSum1[i] = 0;
    radSum2[i] = 0;
//    radNum[i] = 0;
  }

  for (size_t i=0; i<numel; i++)
	if (vol3D1[i] > badVal+1 && vol3D2[i] > badVal+1)
	{ double sum  = vol3D1[i]*vol3D2[i];
	  double sum1 = vol3D1[i]*vol3D1[i];
	  double sum2 = vol3D2[i]*vol3D2[i];
      _sum += sum;
      _sum1 += sum1;
      _sum2 += sum2;
      int bin = radAr[i];
      radCor[bin] += sum;
      radSum1[bin] += sum1;
      radSum2[bin] += sum2;
//      radNum[bin]++;
	}

  if (_sum1*_sum2>MinVal)
	_sum /= sqrt(_sum1*_sum2);
  else _sum = 0;
  for (int i=0; i<radSize; i++)
    if (radSum1[i]*radSum2[i]>MinVal)
      radCor[i] /= sqrt(radSum1[i]*radSum2[i]);
    else radCor[i] = 0;

//  delete[] radNum;
  delete[] radSum1;
  delete[] radSum2;
  return (float)_sum;
}



bool SubLocalBG3D(float* inpAr, size_t dimX, size_t dimY, size_t dimZ, int rad, float badVal, float* smoothedAr)
{
  size_t numComp = (size_t)dimX*(size_t)dimY*(size_t)dimZ;
  float* arrOld = new float[numComp];
//  float* arr1 = new float[(2*rad+1)*(2*rad+1)*(2*rad+1)];

  for (size_t i=0; i<numComp; i++)
    arrOld[i] = inpAr[i];

  if (smoothedAr != NULL)
    for (size_t i=0; i<numComp; i++)
      smoothedAr[i] = badVal;

#pragma omp parallel for shared(inpAr, arrOld, smoothedAr)
  for (long zi=0; zi<dimZ; zi++)
  { float* arr1 = new float[(2*rad+1)*(2*rad+1)*(2*rad+1)];
    for (long yi=0; yi<dimY; yi++)
      for (long xi=0; xi<dimX; xi++)
      { int nuel = 0;
        if (arrOld[xi+dimX*(yi+dimY*zi)]<badVal+1) continue;
        for (int bxi=-rad; bxi<=rad; bxi++)
         for (int byi=-rad; byi<=rad; byi++)
          for (int bzi=-rad; bzi<=rad; bzi++)
            if (xi+bxi>=0 && xi+bxi<dimX && yi+byi>=0 && yi+byi<dimY && zi+bzi>=0 && zi+bzi<dimZ)
              if (arrOld[xi+bxi+dimX*(yi+byi+dimY*(zi+bzi))]>badVal+1)
              { arr1[nuel] = arrOld[xi+bxi+dimX*(yi+byi+dimY*(zi+bzi))];
                nuel++;
              }
        if (nuel>0)
        { //+float theval = quick_select(arr1, nuel);
          //nth_element()
          float theval = medianCutoffStd(arr1, nuel, 0.5);
          if (smoothedAr != NULL)
            smoothedAr[xi+dimX*(yi+dimY*zi)] = theval;
          else
            inpAr[xi+dimX*(yi+dimY*zi)] = arrOld[xi+dimX*(yi+dimY*zi)] - theval;
        } //else
          //printf("was ist das?\n");
      }
    delete arr1;
  }
//  delete arr1;
  delete arrOld;
  return true;
}

bool SubLocalBG(float* inpAr, int stx, int enx, int sty, int eny, int fNumX, int radX, int radY, float badVal, float* smoothedAr)
{
  float* arr1 = new float[(2*radX+1)*(2*radY+1)];
  int numX = enx-stx;
  int numY = eny-sty;
  size_t numComp = (size_t)numX*(size_t)numY;
  float* arrOld = new float[numComp];
  for (int xi=0; xi<numX; xi++)
    for (int yi=0; yi<numY; yi++)
      arrOld[xi+numX*yi] = inpAr[stx+xi+fNumX*(sty+yi)];

  if (smoothedAr != NULL)
    for (size_t i=0; i<numComp; i++)
      smoothedAr[i] = badVal;
//  for (size_t i=0; i<numComp; i++)
//    arrOld[i] = inpAr[i];

  for (int yi=0; yi<numY; yi++)
    for (int xi=0; xi<numX; xi++)
    { int nuel = 0;
      if (arrOld[xi+numX*yi]<badVal+1) continue;
      for (int bxi=-radX; bxi<=radX; bxi++)
        for (int byi=-radY; byi<=radY; byi++)
          if (xi+bxi>=0 && xi+bxi<numX && yi+byi>=0 && yi+byi<numY)
            if (arrOld[xi+bxi+numX*(yi+byi)]>badVal+1)
            { arr1[nuel] = arrOld[xi+bxi+numX*(yi+byi)];
              nuel++;
            }
      if (nuel>0)
      { //+float theval = quick_select(arr1, nuel);
        float theval = medianCutoff(arr1, nuel, 0.5);
        if (smoothedAr != NULL)
          smoothedAr[stx+xi+fNumX*(sty+yi)] = theval;
        else
          inpAr[stx+xi+fNumX*(sty+yi)] = arrOld[xi+numX*yi] - theval;
      }
    }
  delete arr1;
  delete arrOld;
  return true;
}

bool SubLocalBGtiled(float* inpAr, int numX, int numY, int tileX, int tileY, int radX, int radY, float badVal, float* smoothedAr)
{
  size_t numComp = (size_t)numX*(size_t)numY;
  if (smoothedAr != NULL)
    for (int i=0; i<numComp; i++)
      smoothedAr[i] = badVal;
  else return false;

  int nTilesX = roundMy(numX/(float)tileX);
  int nTilesY = roundMy(numY/(float)tileY);

  for (int npy=0; npy<nTilesY; npy++)
    for (int npx=0; npx<nTilesX; npx++)
    { int stx = npx*tileX;
      int enx = (npx+1)*tileX;
      int sty = npy*tileY;
      int eny = (npy+1)*tileY;

      SubLocalBG(inpAr, stx, enx, sty, eny, numX, radX, radY, badVal, smoothedAr);
    }
  return true;
}

float CalculateAverage(float* data, char* mask, size_t numel, float badVal, char goodVal)
{ float sum = 0.;
  size_t num = 0;
  for (size_t i=0; i<numel; i++)
    if (data[i]>badVal+1 && (mask!=NULL?mask[i]==goodVal:true))
    { sum += data[i];
      num++;
    }
  if (num>0) sum /= (float)num;
  return sum;
}

size_t DilateMask(float* data, char *mask, char* inpMask, int tileX, int tileY, int nTilesX, int nTilesY, float badVal, int maskRad)
{ size_t nummasked = 0;
  size_t maxFS = nTilesX*tileX;
  for (int npy=0; npy<nTilesY; npy++)
    for (int npx=0; npx<nTilesX; npx++)
    { int stx = npx*tileX;
      int enx = (npx+1)*tileX;
      int sty = npy*tileY;
      int eny = (npy+1)*tileY;
      for (int yi=sty; yi<eny; yi++)
        for (int xi=stx; xi<enx; xi++)
          if (inpMask[xi+maxFS*yi]>0)
//          if (peakpixel[xi+maxFS*yi]>0)
            for (int iss = -maskRad; iss <= maskRad; iss++)
              for (int ifs = -maskRad; ifs <= maskRad; ifs++)
                if (ifs+xi>=stx && ifs+xi<enx && iss+yi>=sty && iss+yi<eny)
                  if (sqrt(iss*iss+ifs*ifs) < maskRad+MinVal)
                  { if (mask == NULL) data[xi+ifs+maxFS*(yi+iss)] = badVal;
                    else mask[xi+ifs+maxFS*(yi+iss)] = MASK_BRAGG;
                    nummasked++;
                  }
    }
  return nummasked;
}


int MaskRings(float* data, char *mask, int numPoRad, int** radShells, int* radialNumComp, float badVal, float ringDiff, int smF)
{ int smoothF = smF;
  if (smoothF<MinVal) smoothF = numPoRad/30; //???
  float difference = ringDiff;//0.2;

  // make radial average curve
  float* radaver = new float[numPoRad];
  int* irad = new int[numPoRad];
  for (int ring=0; ring<numPoRad; ring++)
  { radaver[ring] = 0.;
    irad[ring] = 0;
    for (int pix=0; pix<radialNumComp[ring]; pix++)
    { if (mask != NULL)
      { if (mask[radShells[ring][pix]]!=MASK_GOOD) continue;
      } else
        if (data[radShells[ring][pix]]<badVal+1) continue;
      radaver[ring] += data[radShells[ring][pix]];
      irad[ring]++;
    }
    if (irad[ring]>0) radaver[ring] /= (float)irad[ring];
    else radaver[ring] = badVal-1;
  }
  // smooth it

  float* smradaver = new float[numPoRad];
  for (int ring=0; ring<numPoRad; ring++) smradaver[ring] = radaver[ring];
  MedianFilter1D(smradaver, numPoRad, smoothF, badVal);

  // mask all rings where difference between original and smoothed is too big
  for (int ring=0; ring<numPoRad; ring++)
    if (irad[ring]>0)
      if (fabs(smradaver[ring])>MinVal)
        if (fabs((smradaver[ring]-radaver[ring])/smradaver[ring])>difference)
          for (int pix=0; pix<radialNumComp[ring]; pix++)
            if (mask == NULL) data[radShells[ring][pix]] = badVal;
            else mask[radShells[ring][pix]] = MASK_BAD;

  delete[] irad;
  delete[] radaver;
  delete[] smradaver;
}

int MaskRingsSimple(float* data, char *mask, int *pix_r, int numPo, float badVal, float ringDiff, int smF)
//int MaskRingsSimple(float* data, char *mask, int numPoRad, int** radShells, int* radialNumComp, float badVal, float ringDiff, int smF)
{
  // find max rad
  int numPoRad = 0;
  for (int pix=0; pix<numPo; pix++)
    if (pix_r[pix] > numPoRad)
      numPoRad = pix_r[pix];
  numPoRad++;

  int smoothF = smF;
  if (smoothF<MinVal) smoothF = numPoRad/30; //???
  float difference = ringDiff;//0.2;

  // make radial average curve
  float* radaver = new float[numPoRad];
  int* irad = new int[numPoRad];
  for (int ring=0; ring<numPoRad; ring++)
  { radaver[ring] = 0.;
    irad[ring] = 0;
  }

  //calculate radial averaged curve
  for (int pix=0; pix<numPo; pix++)
  { int ring = pix_r[pix];
    if (mask != NULL)
    { if (mask[pix]!=MASK_GOOD) continue;
    } //else
    if (data[pix]<badVal+1) continue;
    radaver[ring] += data[pix];
    irad[ring]++;
  }
  for (int ring=0; ring<numPoRad; ring++)
    if (irad[ring]>0) radaver[ring] /= (float)irad[ring];
    else radaver[ring] = badVal-1;

  // smooth it
  float* smradaver = new float[numPoRad];
  for (int ring=0; ring<numPoRad; ring++)
    smradaver[ring] = radaver[ring];
  MedianFilter1D(smradaver, numPoRad, smoothF, badVal);

  // mark bad rings in irad
  for (int ring=0; ring<numPoRad; ring++)
    if (irad[ring]>0)
      if (fabs(smradaver[ring])>MinVal)
        if (fabs((smradaver[ring]-radaver[ring])/smradaver[ring])>difference)
          irad[ring] = -1;

  // mask all rings where difference between original and smoothed is too big
  for (int pix=0; pix<numPo; pix++)
    if (irad[pix_r[pix]]<0)
      if (mask == NULL) data[pix] = badVal;
      else mask[pix] = MASK_BAD;

  delete[] irad;
  delete[] radaver;
  delete[] smradaver;
}


//?size_t SimpleDeleteOutliers(float* data, char* mask, int pix_nn, float badVal,
size_t SimpleDeleteOutliers(float* data, int pix_nn, float badVal,
                            unsigned int loopsNum, float hitfinderMinSNR)
{
//! hitfinderMinSNR - ?  it is =8 in cheetah.ini
//  float hitfinderMinSNR = 8;
//  const int loopsNum = 5; // number of loops to make

  // Allocate and zero arrays
  float rthreshold = 1e30;
  float lowthreshold = badVal+1;
  float roffset;
  float rsigma;
  size_t rcount;

  // Compute sigma and average of data values at each radius
  // From this, compute the ADC threshold to be applied at each radius
  // Iterate a few times to reduce the effect of positive outliers (ie: peaks)
  for(unsigned int counter=0; counter<loopsNum; counter++)
  {
    roffset = 0;
    rsigma = 0;
    rcount = 0;

    for (int i=0; i<pix_nn; i++)
      if (data[i]>lowthreshold && data[i]<rthreshold)
      { roffset += data[i];
        rsigma += (data[i]*data[i]);
        rcount += 1;
      }

    if (rcount == 0) break;

    float thisoffset = roffset/(float)rcount;
    roffset = thisoffset;
    float someval = rsigma/(float)rcount - thisoffset*thisoffset;
    if (someval>MinVal) rsigma = sqrt(someval);
    else rsigma = 0;
    rthreshold = roffset + hitfinderMinSNR*rsigma;
    lowthreshold = roffset - hitfinderMinSNR*rsigma;
  }
  if (rcount>0)
  { rcount = 0;
    for (int i=0; i<pix_nn; i++)
      if (!(data[i] < rthreshold && data[i] > lowthreshold))
      { data[i] = badVal;                            // throwing away all bad values
//      { if (mask!=NULL) data[i] = badVal;                            // throwing away all bad values
//        else mask[i] = MASK_BAD;
        rcount++;
      }
  };

  return rcount;
}

int Make2dRadArrayFrom1d(int* rad1d, int numPoRad, int*** radShells, int** radialNumComp)       // not ready!
{ int* _radNum = new int[numPoRad];
  radialNumComp = &_radNum;
  for (int ring=0; ring<numPoRad; ring++)
    ;

}

size_t PeakFinder3D8(size_t dim[3], float *data, char *mask, char* peakPixel,
                  float hitfinderMinSNR, long hitfinderMinPixCount, long hitfinderMaxPixCount,
                  int numPoRad, int** radShells, int* radialNumComp, float badVal)
//?                  float minrad, float maxrad)

{ char MASK0 = 0;    // good pixels
  char MASK1 = 1;    // outliers found but not accounted
  char MASK2 = 2;    // already accounted outliers
  char MASK3 = 3;    // braggs
  char searchR = 1;

  size_t numEl = dim[0]*dim[1]*dim[2];
  size_t* inx = new size_t[numEl]; //array with connected pixels;
  size_t* iny = new size_t[numEl]; //array with connected pixels;
  size_t* inz = new size_t[numEl]; //array with connected pixels;

  char* peakpixel;
  if (peakPixel==NULL)
    peakpixel = new char[numEl];
  else peakpixel = peakPixel;
  for (size_t i=0; i<numEl; i++)
    peakpixel[i] = MASK0;

  // find all outliers on ring-by-ring basis
  int loopsNum = 5;
  size_t totalmasked = 0;
  if (hitfinderMinSNR<MinVal) hitfinderMinSNR = 8;
  // make radial average curve
  for (int ring=0; ring<numPoRad; ring++)
  {
    float* raddata = new float[radialNumComp[ring]];
	int* radnums = new int[radialNumComp[ring]];
    int numel = 0;
    for (int pix=0; pix<radialNumComp[ring]; pix++)
    { if (mask != NULL)
      { if (mask[radShells[ring][pix]]!=MASK_GOOD) continue;
      }
      if (data[radShells[ring][pix]]<badVal+1) continue;

      raddata[numel] = data[radShells[ring][pix]];
      radnums[numel] = radShells[ring][pix];
      numel++;
    }
    size_t nummasked = 0;
    if (numel>0)
      nummasked = SimpleDeleteOutliers(raddata, numel, badVal, loopsNum, hitfinderMinSNR);
    if (nummasked>0)
      for (int i=0; i<numel; i++)
        if (raddata[i] < badVal+1)
          peakpixel[radnums[i]] = MASK1;

    totalmasked += nummasked;
    delete[] raddata;
    delete[] radnums;
  }

  size_t* conPeaks = new size_t[hitfinderMaxPixCount+1];
  size_t numFound = 0; //debug
  size_t numBraggPeaks = 0; //debug
  // search for connected pixels
  for (size_t zi=0; zi<dim[2]; zi++)
    for (size_t yi=0; yi<dim[1]; yi++)
      for (size_t xi=0; xi<dim[0]; xi++)
      { size_t ind = xi+dim[0]*(yi+dim[1]*zi);
        if (peakpixel[ind] != MASK1) continue;
        size_t numCon = 0;
        conPeaks[numCon] = ind;
        peakpixel[ind] = MASK2;
        inx[numCon] = xi;
        iny[numCon] = yi;
        inz[numCon] = zi;
        numCon++;
        numFound++;
        for (size_t cCon=0; cCon<numCon; cCon++) //counter<nat && counter <= hitfinderMaxPixCount
        { size_t cxi = inx[cCon];
          size_t cyi = iny[cCon];
          size_t czi = inz[cCon];
          for (int zii=-searchR; zii<=searchR; zii++)
            for (int yii=-searchR; yii<=searchR; yii++)
              for (int xii=-searchR; xii<=searchR; xii++)
              { //if (xii==0 && yii==0 && zii==0) continue;
                if (cxi+xii<0 || cxi+xii>=dim[0] || cyi+yii<0 || cyi+yii>=dim[1] || czi+zii<0 || czi+zii>=dim[2]) continue;
                ind = (size_t)(cxi+xii)+dim[0]*((size_t)(cyi+yii)+dim[1]*(size_t)(czi+zii));
                if (peakpixel[ind] != MASK1) continue;
                if (numCon <= hitfinderMaxPixCount) conPeaks[numCon] = ind;
                inx[numCon] = cxi+xii;
                iny[numCon] = cyi+yii;
                inz[numCon] = czi+zii;
                peakpixel[ind] = MASK2;
                numCon++;
                numFound++;
              }
        }
        if (numCon>=hitfinderMinPixCount && numCon<=hitfinderMaxPixCount)
        { for (size_t i=0; i<numCon; i++)
            peakpixel[conPeaks[i]] = MASK3;
          numBraggPeaks++;
        }
      }

  delete[] conPeaks;
  size_t actuallyMasked = 0;
  for (size_t i=0; i<numEl; i++)
    if (peakpixel[i]==MASK2)
    { data[i] = badVal-1;
      actuallyMasked++;
    }
    else if (peakpixel[i]==MASK1)
      printf("How???\n");

  printf("%d outliers, %d accounted, %d Bragg peaks, %d masked\n",totalmasked, numFound, numBraggPeaks, actuallyMasked);
  if (peakPixel==NULL)
    delete[] peakpixel;
  delete[] inx;
  delete[] iny;
  delete[] inz;

  return actuallyMasked;

}


size_t ClearOutliersD8(size_t dim[3], float *data, float* killBraggs, float badVal)
{
    int radSize = 0;
    int* radAr = NULL;
    int* numPerShell = NULL;
    int** shellElem = NULL;

    radSize = BuildRadial2DArray3D(dim, NULL, 1, NULL, NULL, NULL, 1, &radAr,
                                    &numPerShell, &shellElem, NULL);
  // here clean the data
    if (radSize>1)
    {
      size_t masked = PeakFinder3D8(dim, data, NULL, NULL, killBraggs[2], roundMy(killBraggs[3]),
                      roundMy(killBraggs[4]), radSize, shellElem, numPerShell, badVal);
      printf("Masked %d pixels\n",masked);
    }
    if (radAr != NULL)
      delete radAr;
    if (numPerShell != NULL)
      delete numPerShell;
    if (shellElem != NULL)
      for (int i=0; i<radSize; i++)
        if (shellElem[i] != NULL)
          delete shellElem[i];

}


size_t MaskOutliersAtRings(float* data, char *mask, int numPoRad, int** radShells,
                           int* radialNumComp, float badVal, float hitfinderMinSNR)
{
//  float difference = ringDiff;//0.2;
  int loopsNum = 5;
  size_t totalmasked = 0;
  if (hitfinderMinSNR<MinVal) hitfinderMinSNR = 8;
  // make radial average curve
  for (int ring=0; ring<numPoRad; ring++)
  {
    float* raddata = new float[radialNumComp[ring]];
    int* radnums = new int[radialNumComp[ring]];
    int numel = 0;
    for (int pix=0; pix<radialNumComp[ring]; pix++)
    { if (mask != NULL)
      { if (mask[radShells[ring][pix]]!=MASK_GOOD) continue;
      }
      if (data[radShells[ring][pix]]<badVal+1) continue;

      raddata[numel] = data[radShells[ring][pix]];
      radnums[numel] = radShells[ring][pix];
      numel++;
    }
    size_t nummasked = 0;
    if (numel>0)
      nummasked = SimpleDeleteOutliers(raddata, numel, badVal, loopsNum, hitfinderMinSNR);
    if (nummasked>0)
      for (int i=0; i<numel; i++)
        if (raddata[i] < badVal+1)
          if (mask == NULL) data[radnums[i]] = badVal;
          else mask[radnums[i]] = MASK_BAD;

    totalmasked += nummasked;
    delete[] raddata;
    delete[] radnums;
  }

  return totalmasked;
}

bool MaskOutliersAtRings3D(float* data, size_t dim[3], char *mask, float badVal, float hitfinderMinSNR)
{ if (hitfinderMinSNR < MinVal) return false;
  { printf("Cleaning the data with SNR threshold %0.2f\n",hitfinderMinSNR);
    int radSize = 0;
    int* radAr = NULL;
    int* numPerShell = NULL;
    int** shellElem = NULL;

    radSize = BuildRadial2DArray3D(dim, NULL, 1, NULL, NULL, NULL, 1, &radAr,
                                    &numPerShell, &shellElem, NULL);
    if (radSize>1)
  // here clean the data
    { size_t masked = MaskOutliersAtRings(data, mask, radSize, shellElem, numPerShell, badVal, hitfinderMinSNR);
      printf("Masked %d pixels\n",masked);
      //UpdateDataFromMask(pattern->Im, cmask, numComp, panelsGap, MASK_BAD);
    }
    if (radAr != NULL)
      delete radAr;
    if (numPerShell != NULL)
      delete numPerShell;
    if (shellElem != NULL)
      for (int i=0; i<radSize; i++)
        if (shellElem[i] != NULL)
          delete shellElem[i];
  }
}

int UpdateMaskFromData(float* inpAr, char* amask, size_t numel, float badVal, char maskVal)
{
  for (size_t i=0; i<numel; i++)
    if (inpAr[i]<badVal+1)
      amask[i] = maskVal;
}

int UpdateDataFromMask(float* inpAr, char* amask, size_t numel, float badVal, char maskVal)
{
  for (size_t i=0; i<numel; i++)
    if (amask[i] == maskVal)
      inpAr[i] = badVal;
}
                                                                // NOT READY
bool DeleteOutliers(float* inpAr, int stx, int enx, int sty, int eny, int fNumX, int radX, int radY, float badVal, float threshold, float* correctedAr)
{
  float* arr1 = new float[(2*radX+1)*(2*radY+1)];
  int numX = enx-stx;
  int numY = eny-sty;
  size_t numComp = (size_t)numX*(size_t)numY;
  float* arrOld = new float[numComp];
  for (int xi=0; xi<numX; xi++)
    for (int yi=0; yi<numY; yi++)
      arrOld[xi+numX*yi] = inpAr[stx+xi+fNumX*(sty+yi)];

  for (int xi=0; xi<numX; xi++)
    for (int yi=0; yi<numY; yi++)
    { int nuel = 0;
      if (arrOld[xi+numX*yi]<badVal+1) continue;
      for (int bxi=-radX; bxi<=radX; bxi++)
        for (int byi=-radY; byi<=radY; byi++)
          if (xi+bxi>=0 && xi+bxi<numX && yi+byi>=0 && yi+byi<numY)
            if (arrOld[xi+bxi+numX*(yi+byi)]>badVal+1)
            { arr1[nuel] = arrOld[xi+bxi+numX*(yi+byi)];
              nuel++;
            }
      if (nuel>0)
      {
        // HERE!!!
//        LoopAverageBG(dataShells[i], numPerShell1[i], bstpReg, 5, 3., 1);         // NOT READY
        float theval = medianCutoff(arr1, nuel, 0.5);
        if (correctedAr != NULL)
          correctedAr[stx+xi+fNumX*(sty+yi)] = theval;
        else
          inpAr[stx+xi+fNumX*(sty+yi)] = arrOld[xi+numX*yi] - theval;
      }
    }
  delete arr1;
  delete arrOld;
  return true;
}

bool DeleteOutliersTiledNotReady(float* inpAr, int numX, int numY, int tileX, int tileY, int radX, int radY, float badVal, float threshold, float* correctedAr)
{
  size_t numComp = (size_t)numX*(size_t)numY;

  if (correctedAr != NULL)
    for (int i=0; i<numComp; i++)
      correctedAr[i] = badVal;

  int nTilesX = roundMy(numX/(float)tileX);
  int nTilesY = roundMy(numY/(float)tileY);

  for (int npx=0; npx<nTilesX; npx++)
    for (int npy=0; npy<nTilesY; npy++)
    { int stx = npx*tileX;
      int enx = (npx+1)*tileX;
      int sty = npy*tileY;
      int eny = (npy+1)*tileY;

      DeleteOutliers(inpAr, stx, enx, sty, eny, numX, radX, radY, badVal, threshold, correctedAr);
    }
  return true;
}

bool DeleteOutliersArray(float* inpAr, int numX, int numY, int tileX, int tileY, int radmedian, float badVal, float tolerance, int numIter)
{
  size_t numComp = (size_t)numX*(size_t)numY;
  float* blurred = new float[numComp];
  float* difference = new float[numComp];

  for (int iter=0; iter<numIter; iter++)
  {
    for (size_t i=0; i<numComp; i++) blurred[i] = badVal;
    SubLocalBGtiled(inpAr, numX, numY, tileX, tileY, radmedian, radmedian, badVal, blurred);

    for (size_t i=0; i<numComp; i++)
      if (inpAr[i]>badVal+1 && blurred[i]>badVal+1)
        difference[i] = inpAr[i] - blurred[i];
      else
        difference[i] = badVal;
  // standard deviation
    float threshold = tolerance*stdDev(difference,numComp,badVal+1);

    for (size_t i=0; i<numComp; i++)
      if (difference[i]>threshold) inpAr[i] = badVal;
  }
  delete difference;
  delete blurred;
  return true;
}

void MinMax(float* inp, int numb, float badVal, float* minV, float* maxV);
bool histogram(float* ar, size_t numel, int* har, int nbins, float minV, float maxV);



float lininterp(float x, float x1, float y1, float x2, float y2)
{ if (fabs(x2-x1)>MinVal)
    return y1+(y2-y1)*(x-x1)/(x2-x1);
  else
    return 0;
}

// found in internet 
float quadrinterp(float x, float x1, float y1, float x2, float y2, float x3, float y3)
{ float a = ((y3-y1)*(x2-x1)-(y2-y1)*(x3-x1))/((x3*x3-x1*x1)*(x2-x1)-(x2*x2-x1*x1)*(x3-x1));
  float b = (y2-y1-a*(x2*x2-x1*x1))/(x2-x1);
  float c = y1-(a*x1*x1+b*x1);
  return a*x*x+b*x+c;
}




// ----------------------- Some other tools -------------------------------

float PolarisationFactorDet(float detx, float dety, float detz, float degree)
{ if ((fabs(detx)<MinVal && fabs(dety)<MinVal)) return 1.;

  float pdist2i = 1/(detx*detx + dety*dety + detz*detz);
  float pol = 1 - (SQR(dety)*(1.-degree) + SQR(detx)*degree)*pdist2i;
  if (fabs(pol)<MinVal) pol=1;

  return pol;
}

void CorrectPolarisation(float* im, int numcomp, float* cox, float* coy, float* coz, float badReg, float poldegree)
{ if (fabs(coz[0])>MinVal)
    for (int i=0; i<numcomp; i++)
      if (im[i]>badReg+1) im[i] /= PolarisationFactorDet(cox[i],coy[i],coz[i],poldegree);
}

bool MakePolarisationArray(float* pol, int numcomp, float* cox, float* coy, float detdist, float poldegree)
{ if (fabs(detdist)<MinVal || pol==NULL) return false;
  for (int i=0; i<numcomp; i++)
    pol[i] = PolarisationFactorDet(cox[i],coy[i],detdist,poldegree);
  return true;
}

void CorrectSolidAngle(float* im, int numcomp, float* cox, float* coy, float* coz, float badReg)
{ //if (fabs(detdist)>MinVal)
    for (int i=0; i<numcomp; i++)
      if (im[i]>badReg+1)
      { float teta = atan2(sqrt(SQR(cox[i])+SQR(coy[i])),coz[i]);
        float sol = cos(teta);
        sol *= SQR(sol);
        if (sol<MinVal) sol=0;
        else sol = 1/sol;
        im[i] *= sol;
      }
}

void MinMax(float* inp, int numb, float badVal, float* minV, float* maxV)
{   *maxV = -1e35;
    *minV = 1e35;
    for (int i=0; i<numb; i++)
      if (inp[i]>badVal+1)
      { if (inp[i]>*maxV) *maxV = inp[i];
        if (inp[i]<*minV) *minV = inp[i];
      }
}

int Minimum(float* inp, int numb)
{   int minNum = 0;
    float minV = 1e35;
    for (int i=0; i<numb; i++)
      if (inp[i]<minV)
      { minNum = i;
        minV = inp[i];
      }
    return minNum;  
}

//histogram(dvbg, nbins = 100, min = 0, max = dvbg_mean*3, locations=vbg)
bool histogram(float* ar, size_t numel, int* har, int nbins, float minV, float maxV)
{ float astep = (maxV-minV)/(float)nbins;
  for (int i=0; i<nbins; i++)
    har[i] = 0;
  for (size_t i=0; i<numel; i++)
  { int cind = (ar[i]-minV)/astep;
    if (cind>=0 && cind<nbins)
      har[cind]++;
  }
  return true;
}

double MaD2 (double A1, double A2)
{ return (A1>A2?A1:A2); }

double MiD2 (double A1, double A2)
{ return (A1<A2?A1:A2); }

int MaxD3 (double A1, double A2, double A3)
{ if (A1>A2) if (A1>A3) return 0;
             else       return 2;
  else if (A2>A3) return 1;
             else return 2;
}

int MinD3 (double A1, double A2, double A3)
{ if (A1<A2) if (A1<A3) return 0;
             else       return 2;
  else if (A2<A3) return 1;
             else return 2;
}

double IntegralMy(double _a, double _b, double _n, double* _Data)
{                                                       // �������� �� �������
  double Sum = 0;
  double dx = fabs(_b - _a) / _n;
  for (int i=0; i<=_n-2; i+=2)
  { Sum = Sum + (_Data[i] + 4*_Data[i+1] + _Data[i+2])*dx/3;
  };
//  for (int i=1; i<=N; i++)
//  { sum = sum + dx*fabs(func[i] + func[i+1])/2;};
  return Sum;
};


//--------------------------------- something else ----------------------------

double gaussian_noise(double expected, double stddev)
{
	double x1, x2, noise;
	/* A uniformly distributed random number between 0 and 1 */
	x1 = (double)(rand()+1.)/((double)RAND_MAX+1.);                  // My mod - to avoid log(0) in noise
	x2 = (double)rand()/(double)RAND_MAX;
	noise = sqrt(-2.0*log(x1)) * cos(2.0*M_PI*x2);
	return expected + noise*stddev;
}

int poisson_noise(double expected)
{
	double L;
	int k = 0;
	double p = 1.0;
	/* For large values of the mean, we get big problems with arithmetic.
	 * In such cases, fall back on a Gaussian with the right variance. */
//?	if ( expected > 100.0 ) return (int)gaussian_noise(expected, sqrt(expected));;
	L = exp(-expected);
	do {
		double r;
		k++;
		r = (double)rand()/(double)RAND_MAX;
		p *= r;
	} while ( p > L );
	return k - 1;
}

float GaussDistrib(float dist, float sigma)
{ float _sig = (sigma>MinVal?1./sigma:0.);
  return exp(-0.5*(SQR(dist*_sig)));
}

void MakeGauss3D(float** gvol, float sigmaX, float sigmaY, float sigmaZ, int dX, int dY, int dZ)
{
  float* _gvol = *gvol;
  float Ampl = 1.;
  float siX = (sigmaX>MinVal?1./sigmaX:0.);
  float siY = (sigmaY>MinVal?1./sigmaY:0.);
  float siZ = (sigmaZ>MinVal?1./sigmaZ:0.);
  for(int i = 0; i<dZ; i++)
	for(int j = 0; j<dY; j++)
	  for(int k = 0; k<dX; k++)
		_gvol[k+dX*(j+dY*i)] = Ampl*exp(-0.5*(SQR((k-0.5*dX)*siX)+SQR((j-0.5*dY)*siY)+SQR((i-0.5*dZ)*siZ)));
};

void Smooth3D(float* gVol, int dX, int dY, int dZ, int smoC)
{ if (dX<1) dX = 1;
  if (dY<1) dY = 1;
  if (dZ<1) dZ = 1;
  int numel = (dX>0?dX:1)*(dY>0?dY:1)*(dZ>0?dZ:1);
  float* _gvol = new float[numel];
  int* _gvolN = new int[numel];
  for(int i = 0; i<numel; i++)
  {	_gvol[i] = gVol[i];
	gVol[i] = 0.;
	_gvolN[i] = 0;
  }
  for(int i = 0; i<dZ; i++)
   for(int j = 0; j<dY; j++)
	for(int k = 0; k<dX; k++)
	  for(int si = -smoC; si<=smoC; si++)
	   for(int sj = -smoC; sj<=smoC; sj++)
		for(int sk = -smoC; sk<=smoC; sk++)
		  if (i+si>=0 && i+si<dZ && j+sj>=0 && j+sj<dY && k+sk>=0 && k+sk<dX)
		  { if (si*si+sj*sj+sk*sk > smoC*smoC) continue;
			gVol[k+dX*(j+dY*i)] += _gvol[k+sk+dX*(j+sj+dY*(i+si))];
			_gvolN[k+dX*(j+dY*i)]++;
		  }
  for(int i = 0; i<numel; i++)
	if (_gvolN[i]>0) gVol[i] /= _gvolN[i];
  delete _gvol;
  delete _gvolN;
};

struct TPoints
{ float x;
  float y;
  float z;
  float I;
  TPoints* ne;
};

int Interpolate3DVolume(TPoints* allPo, float* vol1, float* vol2, int* dim, float fillWith,
                   float phi, float tet, float psi, float radiusF, float radiusE)
{ int size2 = dim[0]*dim[1]*dim[2];
  float cen2[3] = {0.5*(dim[0]-1),0.5*(dim[1]-1),0.5*(dim[2]-1)};
  for (size_t i=0; i<size2; i++) vol2[i] = fillWith;

  TPoints** AllPo = new TPoints*[size2];
  for (int i=0; i<size2; i++) AllPo[i] = NULL;

  TPoints** cPo;

  for (int iz=0; iz<dim[2]; iz++)
   for (int iy=0; iy<dim[1]; iy++)
    for (int ix=0; ix<dim[0]; ix++)
    { double cof[3];
      int coi[3] = {0,0,0};
      RotateEuler(phi, tet, psi, ix-cen2[0], iy-cen2[1], iz-cen2[2],
                      &cof[0], &cof[1], &cof[2]);
      bool outOf = false;
      for (int i=0; i<3; i++)
      { coi[i] = roundMy(cof[i]+cen2[i]);
        if (coi[i]<0 || coi[i]>=dim[i]) outOf = true;
      }
      if (outOf) continue;

//      int _ind = coi[2]+dim[2]*(coi[1]+dim[1]*coi[0]);
      int _ind = coi[0]+dim[0]*(coi[1]+dim[1]*coi[2]);
//-      if (vol1[_ind]<=fillWith) continue;
      float valI = vol1[ix+dim[0]*(iy+dim[1]*iz)];
//      float valI = vol1[iz+dim[2]*(iy+dim[1]*ix)];
      if (valI<=fillWith) continue;

      cPo = &AllPo[_ind];
      while (*cPo != NULL)
        cPo = &(*cPo)->ne;

      *cPo = new TPoints;
      (*cPo)->ne = NULL;
      (*cPo)->z = cof[2]+cen2[2];
      (*cPo)->y = cof[1]+cen2[1];
      (*cPo)->x = cof[0]+cen2[0];
      (*cPo)->I = valI;
//?      (*cPo)->I = vol1[ix+dim[0]*(iy+dim[1]*iz)];
//          vol2[coi[0]+dim[0]*(coi[1]+dim[1]*coi[2])] = vol1[ix+dim[0]*(iy+dim[1]*iz)];
    }

// actual interpolation

  /// radiusE - different radius for empty points!
  float radiusF2 = radiusF*radiusF;
  float radiusE2 = radiusE*radiusE;
  float _timeS = time(0);
//  int* vol3Dnum = new int[size2];
//  for (int i=0; i<size2; i++) vol3Dnum[i]=0;
#pragma omp parallel for
//private(_dist, wi, sumWi, _cI, _cIi, _cIr, cPo)
  for (int zi=0; zi<dim[2]; zi++)
  {
    double _cI, sumWi, wi, _dist;

#ifdef _OPENMP
	printf("Thread %d executes zi=%d of %d at %d\n", omp_get_thread_num(),zi,dim[2],(int)(time(0)-_timeS));
#endif
   for (int yi=0; yi<dim[1]; yi++)
	for (int xi=0; xi<dim[0]; xi++)
	{ sumWi = 0.;
	  _cI = 0.;

	  bool bigRad = false;

      float radius;
      float radius2;
	  if (AllPo[xi+dim[0]*(yi+dim[1]*zi)]==NULL)
      { radius = radiusE;
        radius2 = radiusE2;
	  } else
      { radius = radiusF;
        radius2 = radiusF2;
      }

	  for (int xii=xi-radius; xii<=xi+radius; xii++)
	   for (int yii=yi-radius; yii<=yi+radius; yii++)
		for (int zii=zi-radius; zii<=zi+radius; zii++)
		{ if (xii<0 || yii<0 || zii<0 || xii>=dim[0] || yii>=dim[1] || zii>=dim[2]) continue;
		  cPo = &AllPo[xii+dim[0]*(yii+dim[1]*zii)];
		  while (*cPo != NULL)
		  { _dist = (SQR((*cPo)->x-xi)+SQR((*cPo)->y-yi)+SQR((*cPo)->z-zi));
			if (_dist > radius2)
			{ cPo = &(*cPo)->ne;
			  continue;
			}
			if (_dist<MinVal) // a point at this location found!
			{ sumWi = 1.;
			  _cI = (*cPo)->I;
//			  vol3Dnum[xi+dim[0]*(yi+dim[1]*zi)]=1;
			  xii=dim[0]; yii=dim[1]; zii=dim[2];
			  break; // must be up to sumWi=0
			}
			_dist = sqrt(_dist);

			wi = SQR((radius-_dist)/(radius*_dist));
			sumWi += wi;
			_cI += wi * (*cPo)->I;
//			vol3Dnum[xi+dim[0]*(yi+dim[1]*zi)]++;
			cPo = &(*cPo)->ne;
		  }
		}
	  if (sumWi>MinVal) //???
        vol2[xi+dim[0]*(yi+dim[1]*zi)] = _cI/sumWi;
	  else
        vol2[xi+dim[0]*(yi+dim[1]*zi)] = fillWith;
	}
  }

//  delete vol3Dnum;
  delete AllPo;
  return 1;
}

