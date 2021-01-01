#include <math.h>
#include <stdio.h>
#include <stdlib.h>
/*void dump_sat_();
  To avoid numer. problems it should have been compiled with flags like  --dbl --long

  In one example the reduction of ~100 in the storage is achieved as
  compared to the full data stored if the line with != below is
  used. Care must be taken in implementing this storage.

 */
void dump_sat_(double *bgy_sat, long *N, long *npf, long *nmu, long *mtl_m1, long *ntae, double *wtae, double *qpf, double *pmu, double *dwp_sat, double *dwv_sat, double *vres_sat, double *fdist_sat, double *gam_sat){
  long i=1, j=2;
  FILE *fpn;
/* use vres.dat to store the matrixes */
  char flnm[12]="Out/vres.dat";
  fpn=fopen(flnm,"w");
  fwrite(npf,sizeof(long),i,fpn);
  fwrite(nmu,sizeof(long),i,fpn);
  fwrite(mtl_m1,sizeof(long),i,fpn);
  fwrite(&j,sizeof(long),i,fpn);
  fwrite(ntae,sizeof(long),i,fpn);
  fwrite(wtae,sizeof(double),i,fpn);
  fwrite(qpf,sizeof(double),*npf,fpn);
  fwrite(pmu,sizeof(double),*nmu,fpn);
/* Here N is npf*nmu*(2*mtl_m1+1)*2*2 */
  for(j=0; j<(*N); j++){
    /*    if(gam_sat[j]!=0){*/
      fwrite(&j,sizeof(long),i,fpn);
      fwrite(dwv_sat+j,sizeof(double),i,fpn);
      fwrite(vres_sat+j,sizeof(double),i,fpn);
      fwrite(fdist_sat+j,sizeof(double),i,fpn);
      fwrite(gam_sat+j,sizeof(double),i,fpn);
      fwrite(dwp_sat+j,sizeof(double),i,fpn);
      fwrite(bgy_sat+j,sizeof(double),i,fpn);
      /*    }*/
  }
  fclose(fpn);
}
