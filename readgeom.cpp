#include "MyMath.h"


float convertAnalyzedNumber(char* astr)
{   float x = 0;
    if (strlen(astr)<1) x = 1;
    else if (strlen(astr)==1 && astr[0]=='-') x = -1;
    else if (strlen(astr)==1 && astr[0]=='+') x = 1;
    else x = atof(astr);
    return x;
}

bool analyzeNumber(char* astr, float* x, float* y, float* z)
{ if (astr==NULL) return false;
  if (strchr(astr,'x'))
  { strchr(astr,'x')[0] = 0;
    *x = convertAnalyzedNumber(astr);
  } else if (strchr(astr,'y'))
  { strchr(astr,'y')[0] = 0;
    *y = convertAnalyzedNumber(astr);
  } else if (strchr(astr,'z'))
  { strchr(astr,'z')[0] = 0;
    *z = convertAnalyzedNumber(astr);
  } else
    return false;
  return true;
}
bool processFsSs(char* str, float* x, float* y, float* z)
{
  analyzeNumber(strtok(str," \n"), x, y, z);
  analyzeNumber(strtok(NULL," \n"), x, y, z);
  analyzeNumber(strtok(NULL," \n"), x, y, z);
  return true;
}

struct badregion
{
	char name[256];
	double min_x;
	double max_x;
	double min_y;
	double max_y;
	double min_fs;
	double max_fs;
	double min_ss;
	double max_ss;
	char local_coord; //==1 - use ss,fs, otherwise x,y
	int is_fsss;

	char *panel;

};

struct panel
{
  char     name[1024];  /* Name for this panel */
  int      min_fs;  /* Smallest FS value considered to be in the panel */
  int      max_fs;  /* Largest FS value considered to be in this panel */
  int      min_ss;  /* ... and so on */
  int      max_ss;
  double   cnx;       /* Location of corner (min_fs,min_ss) in pixels */
  double   cny;
  double   coffset;
  double   clen;     /* Camera length in metres */
  char    *clen_from;
  double   res;      /* Resolution in pixels per metre */
  double   adu_per_eV;   /* Number of ADU per eV */
  double   max_adu;  /* Treat pixel as unreliable if higher than this */

  double fsx;
  double fsy;
  double ssx;
  double ssy;

  char *mask;
  char *data;    //like this /entry_1/data_1/data

  int w;  // Width, calculated as max_fs-min_fs+1
  int h;  // Height, calculated as max_ss-min_ss+1

};

struct detector
{
  struct panel     *panels;
  int               n_panels;

  int               max_fs;
  int               max_ss;  /* Size of overall array needed, minus 1 */

  struct badregion *bad;
  int               n_bad;

  char              *data;    //like this /entry_1/data_1/data
  char              *mask;

  char              *dim0;
  char              *dim1;
  char              *dim2;
  char              *dim3;
  int               *dimS;

  unsigned int       mask_bad;
  unsigned int       mask_good;

  panel       defaults;

  float              coffsetD;
  float              clenD;
  char              *clen_fromD;
  char              *photon_energy;
  float              max_aduD;
  float              adu_per_eVD;        
  int                det_type;
};


bool ReadGeometry(char* FilNam, float** arxi, float** aryi, float** arzi,
       float* istep, int* numfs, int* numss, int* numcomp, float* Kwp, char* datafield, char* maskfield, 
       int* frameaxis, float* lessthan, float* morethan, float* equalto)
{ // I need to copy from crystfel a function for conversion fs,ss -> x,y
  // I can use NumOfDifs for indicating which geometry array to use foe a diffr pattern
  // Get number of components in geometry file and dimentions

  int maxFS = 0;
  int maxSS = 0;
  FILE* stGeom = fopen(FilNam, "rt");
  if (stGeom==NULL) 
  { printf("File %s not oppened!\n", FilNam); 
    return false;
  }
  char str1[MaxStrLen];
  char str2[MaxStrLen];
  int minfs, maxfs, minss, maxss;
  minfs = maxfs = minss = maxss = 0;
  int numpix = 0;
  int numpanels = 0;
  int numbads = 0;
  char cpanel[MaxStrLen] = "";
  char cclen[MaxStrLen] = "";
  char eventdata[MaxStrLen] = "";
  char eventmask[MaxStrLen] = "";
  char cbadreg[MaxStrLen] = "qq";
  float aduev = 0;
//  float lessthan = 0;
//  float morethan = 0;
//  float equalto = 0;


  while (!feof(stGeom))
  { fgets(str1,MaxStrLen,stGeom);
    if ((str1[0]==';' || str1[0]=='#' || str1[0]=='\%' || str1[0]=='\n' || (str1[0]=='/' && str1[1]=='/') || str1[0]==0) && !feof(stGeom)) 
      continue;
    if (strncmp(str1,"adu_per_",8)==0)
      aduev = atof(strchr(str1,'=')+1);
    if (!feof(stGeom))
    { if (strchr(str1,'/') == NULL ) continue;
      strcpy(str2, strchr(str1,'/')+1);
      *strchr(str1,'/') = '\0';
      if (strchr(str1,'=') != NULL ) continue;
    }
    if ((strcmp(str1,cpanel)!=0 || feof(stGeom)) && strncmp(str1,"bad",3)!=0)
    { strcpy(cpanel,str1);
      if (maxfs-minfs>0 || maxss-minss>0)
      { for (int iss=0; iss <= maxss-minss; iss++)
          for (int ifs=0; ifs <= maxfs-minfs; ifs++)
            numpix++;
        if (maxfs>maxFS) maxFS = maxfs;
        if (maxss>maxSS) maxSS = maxss;
        numpanels++;
      }
      minfs = maxfs = minss = maxss = 0;
    }
    if (strncmp(str2,"min_fs",6)==0)
    minfs = atoi(strchr(str2,'=')+1);
    else if (strncmp(str2,"min_ss",6)==0)
    minss = atoi(strchr(str2,'=')+1);
    else if (strncmp(str2,"max_fs",6)==0)
    maxfs = atoi(strchr(str2,'=')+1);
    else if (strncmp(str2,"max_ss",6)==0)
    maxss = atoi(strchr(str2,'=')+1);
    // count bad regions
    if (strncmp(str1,"bad",3)==0)
    { if (strcmp(str1,cbadreg)!=0)
      { strcpy(cbadreg,str1);
        numbads++;
      }
      minfs = maxfs = minss = maxss = 0;
      continue;
    }
  }
  rewind(stGeom);
  maxFS+=1;
  maxSS+=1;
  *numfs = maxFS;
  *numss = maxSS;
  if (maxSS*maxFS != numpix)
  { printf("Something wrong in geometry: MaxSS*MaxFS != NumPix, %d*%d != %d\n",maxSS,maxFS,numpix);
    if (maxSS*maxFS > numpix)
    { printf("  Setting NumPix = MaxSS*MaxFS\n");
      numpix = maxSS*maxFS;
    }
  }
  strcpy(cpanel,"");

  // Creating Tom's det structure
  detector* det = new detector;
  det->n_panels = numpanels;
  det->max_fs = maxFS-1;
  det->max_ss = maxSS-1;
  det->panels = new panel[numpanels];
  det->n_bad = numbads;                       //??? maybe divide by 4?
  det->bad = new badregion[det->n_bad];
  det->mask_good = 0;
  det->mask_bad = 0;
  for (int bi=0; bi<det->n_bad; bi++)
  { strcpy(det->bad[bi].name,"");
    det->bad[bi].min_x = 0;
    det->bad[bi].min_y = 0;
    det->bad[bi].max_x = 0;
    det->bad[bi].max_y = 0;
    det->bad[bi].min_fs = 0;
    det->bad[bi].min_ss = 0;
    det->bad[bi].max_fs = 0;
    det->bad[bi].max_ss = 0;
    det->bad[bi].local_coord = 0;
    det->bad[bi].panel = NULL;
  }
  numbads = 0;

  *numcomp = numpix;
  *arxi = new float[numpix];
  *aryi = new float[numpix];
  *arzi = new float[numpix];

  float* arx = *arxi;
  float* ary = *aryi;
  float* arz = *arzi;
  minfs = maxfs = minss = maxss = 0;
  float pixel, fsx, fsy, fsz, ssx, ssy, ssz, cornx, corny, cx, cy, coffset;
  coffset = 0;
//?  float peaksep = 0;
  int noindex = 0;
  char cbadraw = '-';
  float globpixel = 1.;
  pixel = globpixel;
  fsx = fsy = fsz = ssx = ssy = ssz = cornx = corny = coffset = 0.;
  numpix = 0;
  numpanels = 0;
  det->coffsetD = 0;
  det->adu_per_eVD = 0;
  det->max_aduD = 0;
  det->clen_fromD = new char[MaxStrLen];
  det->clen_fromD[0] = 0;
  det->data = new char[MaxStrLen];
  det->data[0] = 0;
  det->mask = new char[MaxStrLen];
  det->mask[0] = 0;
  det->dim0 = new char[MaxStrLen];
  det->dim0[0] = 0;
  det->dim1 = new char[MaxStrLen];
  det->dim1[0] = 0;
  det->dim2 = new char[MaxStrLen];
  det->dim2[0] = 0;
  det->dim3 = new char[MaxStrLen];
  det->dim3[0] = 0;
  det->dimS = new int[4];
  for (int i=0; i<4; i++)
    det->dimS[i] = i-1;          // this is correct only for a single h5 file
  det->dimS[4] = -1;
  det->photon_energy = new char[MaxStrLen];
  det->photon_energy[0] = 0;
  det->clenD = 0;
  strcpy(det->clen_fromD, "");
  strcpy(str1, "");
  strcpy(str2, "");
  int cnp;
  while (true)//(!feof(stGeom))
  { fgets(str1,MaxStrLen,stGeom);
//?    TrimNoCh(str1);
    if ((str1[0]==';' || str1[0]=='#' || str1[0]=='\%' || str1[0]=='\n' || (str1[0]=='/' && str1[1]=='/') || str1[0]==0) && !feof(stGeom)) 
      continue;
    if (!feof(stGeom))
    {
      if (strncmp(str1,"coffset",7)==0) //Global coffset
      { det->coffsetD = atof(strchr(str1,'=')+1);
        continue;
      }
      if (strncmp(str1,"max_adu",7)==0) //Global max adu
      { det->max_aduD = atof(strchr(str1,'=')+1);
        continue;
      }
      if (strncmp(str1,"adu_per_eV",9)==0) //Global adu per ev
      { det->adu_per_eVD = atof(strchr(str1,'=')+1);
        continue;
      }
      if (strncmp(str1,"clen",4)==0) //Global clen
      { strcpy(det->clen_fromD, strchr(str1,'=')+1);
        TrimNoCh(det->clen_fromD);
        if (isdigit(det->clen_fromD[0]) || det->clen_fromD[0]=='-' || det->clen_fromD[0]=='.' || det->clen_fromD[0]=='+')
          det->clenD = atof(det->clen_fromD);
        continue;
      }
      if (strncmp(str1,"res ",4)==0) //Global resolution
      { globpixel = atof(strchr(str1,'=')+1);
        pixel = globpixel;
        continue;
      }
      if (strncmp(str1,"data ",5)==0) //Global data
      { strcpy(eventdata,strchr(str1,'=')+1);
        TrimNoCh(eventdata);
        strcpy(det->data,eventdata);
        continue;
      }
      if (strncmp(str1,"dim0 ",5)==0) //dimentions
      { strcpy(eventdata,strchr(str1,'=')+1);
        TrimNoCh(eventdata);
        strcpy(det->dim0,eventdata);
        continue;
      }
      if (strncmp(str1,"dim1 ",5)==0) //dimentions
      { strcpy(eventdata,strchr(str1,'=')+1);
        TrimNoCh(eventdata);
        strcpy(det->dim1,eventdata);
        continue;
      }
      if (strncmp(str1,"dim2 ",5)==0) //dimentions
      { strcpy(eventdata,strchr(str1,'=')+1);
        TrimNoCh(eventdata);
        strcpy(det->dim2,eventdata);
        continue;
      }
      if (strncmp(str1,"dim3 ",5)==0) //dimentions
      { strcpy(eventdata,strchr(str1,'=')+1);
        TrimNoCh(eventdata);
        strcpy(det->dim3,eventdata);
        continue;
      }
      if (strncmp(str1,"mask ",5)==0) //Global mask
      { strcpy(det->mask,strchr(str1,'=')+1);
        TrimNoCh(det->mask);
        continue;
      }
      if (strncmp(str1,"photon_energy",12)==0) //Global max adu
      { strcpy(det->photon_energy,strchr(str1,'=')+1);
        TrimNoCh(det->photon_energy);
        continue;
      }
      if (strncmp(str1,"mask_good",9)==0) //Global max adu               //!!!!!!!!!
      { sscanf(strchr(str1,'=')+1,"%x",&det->mask_good);
        continue;
      }
      if (strncmp(str1,"mask_bad",8)==0) //Global max adu
      { sscanf(strchr(str1,'=')+1,"%x",&det->mask_bad);
        continue;
      }
      if (strncmp(str1,"rigid_group",11)==0) //Global max adu    rigid groups
      { 
        continue;
      }
      //flag_lessthan, flag_morethan, flag_equal
      if (strncmp(str1,"flag_lessthan",12)==0)
      { *lessthan = atof(strchr(str1,'=')+1);
        continue;
      }
      if (strncmp(str1,"flag_morethan",12)==0)
      { *morethan = atof(strchr(str1,'=')+1);
        continue;
      }
      if (strncmp(str1,"flag_equal",9)==0)
      { *equalto = atof(strchr(str1,'=')+1);
        continue;
      }


      if (strchr(str1,'/') != NULL ) //continue;
      { strcpy(str2, strchr(str1,'/')+1);
        *strchr(str1,'/') = '\0'; //now str1 - the first half of the string (before /) and str2 - the second
      } else strcpy(str1,"");
    }
    // take into account bad regions
    if (strncmp(str1,"bad",3)==0)
    { int _cbad = -1;
      for (int bi=0; bi<numbads; bi++)
        if (strncmp(str1,det->bad[bi].name,strlen(str1))==0)
        { _cbad = bi;
          break;
        }
      if (_cbad<0)
      { strcpy(det->bad[numbads].name,str1);
        numbads++;
        _cbad = numbads-1;
      }
      if (numbads>det->n_bad)
      { numbads=det->n_bad;
        if (!feof(stGeom)) continue;
      }
      // analysis of the second part
      if (strncmp(str2,"min_x",5)==0)
      det->bad[_cbad].min_x = atoi(strchr(str2,'=')+1);
      else if (strncmp(str2,"max_x",5)==0)
      det->bad[_cbad].max_x = atoi(strchr(str2,'=')+1);
      else if (strncmp(str2,"min_y",5)==0)
      det->bad[_cbad].min_y = atoi(strchr(str2,'=')+1);
      else if (strncmp(str2,"max_y",5)==0)
      det->bad[_cbad].max_y = atoi(strchr(str2,'=')+1);
      else if (strncmp(str2,"min_fs",6)==0)
      det->bad[_cbad].min_fs = atoi(strchr(str2,'=')+1);
      else if (strncmp(str2,"max_fs",6)==0)
      det->bad[_cbad].max_fs = atoi(strchr(str2,'=')+1);
      else if (strncmp(str2,"min_ss",6)==0)
      det->bad[_cbad].min_ss = atoi(strchr(str2,'=')+1);
      else if (strncmp(str2,"max_ss",6)==0)
      det->bad[_cbad].max_ss = atoi(strchr(str2,'=')+1);
      if (det->bad[_cbad].max_ss>0 || det->bad[_cbad].max_fs>0 ||
          det->bad[_cbad].min_ss>0 || det->bad[_cbad].min_fs>0)
        det->bad[_cbad].local_coord = 1;

      if (!feof(stGeom)) continue;
    }

    if (strcmp(str1,cpanel)!=0 || feof(stGeom))
    { if (maxfs-minfs>0 || maxss-minss>0)
      { numpanels++;
        cnp = numpanels-1;
        strcpy(det->panels[cnp].name, cpanel);
        TrimNoCh(eventdata);
        TrimNoCh(eventmask);
        TrimNoCh(cclen);
        det->panels[cnp].clen_from = new char[strlen(cclen)];
        det->panels[cnp].data = new char[strlen(eventdata)];
        det->panels[cnp].mask = new char[strlen(eventmask)];
        strcpy(det->panels[cnp].clen_from, cclen);
        if (isalnum(det->panels[cnp].clen_from[0]))
          det->panels[cnp].clen = atof(det->panels[cnp].clen_from);
        else det->panels[cnp].clen = 0;
        strcpy(cclen,"");
        strcpy(det->panels[cnp].data, eventdata);
        strcpy(det->panels[cnp].mask, eventmask);
        det->panels[cnp].min_fs = minfs;
        det->panels[cnp].min_ss = minss;
        det->panels[cnp].max_fs = maxfs;
        det->panels[cnp].max_ss = maxss;
        det->panels[cnp].cnx = cornx;
        det->panels[cnp].cny = corny;
        det->panels[cnp].fsx = fsx;
        det->panels[cnp].fsy = fsy;
        det->panels[cnp].ssx = ssx;
        det->panels[cnp].ssy = ssy;
//                printf("fsx=%0.2f, fsy=%0.2f, ssx=%0.2f, ssy=%0.2f\n",fsx,fsy,ssx,ssy);
        det->panels[cnp].res = pixel;
        det->panels[cnp].adu_per_eV = aduev;
        det->panels[cnp].coffset = coffset;
        det->panels[cnp].w = maxfs-minfs+1;
        det->panels[cnp].h = maxss-minss+1;
        *istep = pixel;
        for (int iss=0; iss <= maxss-minss; iss++)
          for (int ifs=0; ifs <= maxfs-minfs; ifs++)
          { cx = (cornx+(ifs+0.0)*fsx+(iss+0.0)*ssx)/pixel; //0.5 to have coordinates of the center of each pixel
            cy = (corny+(ifs+0.0)*fsy+(iss+0.0)*ssy)/pixel;
//  float xx = (ifs*fsx + iss*ssx + cornx) / det->panels[cnp].res;
//  float yy = (ifs*fsy + iss*ssy + corny) / det->panels[cnp].res;
            arx[(iss+minss)*maxFS+ifs+minfs] = cx;
            ary[(iss+minss)*maxFS+ifs+minfs] = cy;
            arz[(iss+minss)*maxFS+ifs+minfs] = (fabs(det->panels[cnp].clen)>MinVal?det->panels[cnp].clen:(fabs(det->clenD)>MinVal?det->clenD:BstpReg-1));
            if (arz[(iss+minss)*maxFS+ifs+minfs] > BstpReg-1)
              arz[(iss+minss)*maxFS+ifs+minfs] += (fabs(coffset)>MinVal?coffset:det->coffsetD);
//test_CSPAD
//            arco[(iss+minss)*maxFS+ifs+minfs] = cnp;
            numpix++;
          }
      }
      strcpy(cpanel,str1);
      minfs = maxfs = minss = maxss = 0;
      pixel = globpixel;
      noindex = 0;
      fsx = fsy = fsz = ssx = ssy = ssz = cornx = corny = coffset = 0.;
    }
    if (strncmp(str2,"min_fs",6)==0)
    minfs = atoi(strchr(str2,'=')+1);
    else if (strncmp(str2,"min_ss",6)==0)
    minss = atoi(strchr(str2,'=')+1);
    else if (strncmp(str2,"max_fs",6)==0)
    maxfs = atoi(strchr(str2,'=')+1);
    else if (strncmp(str2,"max_ss",6)==0)
    maxss = atoi(strchr(str2,'=')+1);
    else if (strncmp(str2,"no_index",8)==0)
    noindex = atoi(strchr(str2,'=')+1);
    else if (strncmp(str2,"clen",4)==0)
    strcpy(cclen,strchr(str2,'=')+1);
    else if (strncmp(str2,"data",4)==0)
    strcpy(eventdata,strchr(str2,'=')+1);
    else if (strncmp(str2,"mask",4)==0)
    strcpy(eventmask,strchr(str2,'=')+1);
    else if (strncmp(str2,"badraw",6)==0)
    cbadraw = strchr(str2,'=')[2];
    else if (strncmp(str2,"adu_per_",8)==0)
    aduev = atof(strchr(str2,'=')+1);
    else if (strncmp(str2,"res ",4)==0)
    { pixel = atof(strchr(str2,'=')+1);
      if (fabs(pixel)<1e-10) 
        pixel = globpixel;
    }
    else if (strncmp(str2,"corner_x",8)==0)
    cornx = atof(strchr(str2,'=')+1);
    else if (strncmp(str2,"corner_y",8)==0)
    corny = atof(strchr(str2,'=')+1);
    else if (strncmp(str2,"coffset",7)==0)
    coffset = atof(strchr(str2,'=')+1);
    else if (strncmp(str2,"fs ",3)==0)
    { strcpy(str2, strchr(str2,'=')+1);
      TrimNoCh(str2);
      processFsSs(str2, &fsx, &fsy, &fsz);
      if (fsx*fsx+fsy*fsy>1.1 || fsx*fsx+fsy*fsy<0.9) printf("Wrong at %s fsx=%0.2f and fsy=%0.2f\n",cpanel,fsx,fsy);
    }
    else if (strncmp(str2,"ss ",3)==0)
    { strcpy(str2, strchr(str2,'=')+1);
      TrimNoCh(str2);
      processFsSs(str2, &ssx, &ssy, &ssz);
      if (ssx*ssx+ssy*ssy>1.1 || ssx*ssx+ssy*ssy<0.9) printf("Wrong at %s ssx=%0.2f and ssy=%0.2f\n",cpanel,ssx,ssy);
    }
    if (feof(stGeom)) break;
  }
  fclose(stGeom);

  //correct dims:
  det->dimS[0] = (det->dim0[0]=='\%'?0:(det->dim1[0]=='\%'?1:(det->dim2[0]=='\%'?2:-1)));
  det->dimS[1] = (strstr(det->dim0,"ss")!=NULL?0:(strstr(det->dim1,"ss")!=NULL?1:(strstr(det->dim2,"ss")!=NULL?2:0)));
  det->dimS[2] = (strstr(det->dim0,"fs")!=NULL?0:(strstr(det->dim1,"fs")!=NULL?1:(strstr(det->dim2,"fs")!=NULL?2:1)));
//  printf("dims are: %d,%d,%d\n", det->dimS[0], det->dimS[1], det->dimS[2]);

  if (det->adu_per_eVD<MinVal) det->adu_per_eVD = aduev;
  for (int i=0; i<numpix; i++)
    if (arz[i]<MinVal)
    { arz[0] = BstpReg-1;
      break;
    }

  //change units for bad regions
  for (int _cbad = 0; _cbad < det->n_bad; _cbad++)
  { det->bad[_cbad].min_x /= *istep;
    det->bad[_cbad].max_x /= *istep;
    det->bad[_cbad].min_y /= *istep;
    det->bad[_cbad].max_y /= *istep;
  }
  // I don't know what is it for, taken from crystfel
  //* Calculate matrix inverse *
  for (int i=0; i<det->n_panels; i++)
  { struct panel *p;
    double d;
    p = &det->panels[i];
    if ( p->fsx*p->ssy == p->ssx*p->fsy )
      printf("Panel %i transformation singular.\n", i);
  }

  float energ = atof(det->photon_energy);
  *Kwp = energ/12398.4;
  strcpy(datafield, det->data);
  strcpy(maskfield, det->mask);
  if (det->dim0[0]=='%') *frameaxis = 0;
  else if (det->dim1[0]=='%') *frameaxis = 1;
  else if (det->dim2[0]=='%') *frameaxis = 2;
  else if (det->dim3[0]=='%') *frameaxis = 3;

  delete[] det;
  return true;
}
