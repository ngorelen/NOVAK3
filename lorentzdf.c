#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define npl 200
/* number of Legender polynomials to expand the source */

void sal_(), lorentzdf_();
double fun4_eq_(), gfl_();
/*
************************************************************************
   Compute the distribution function at one point in pitch angle. This is not efficient, but selfcontained. 
*/
void lorentzdf_(double *g1cs, double *gL0, double *gc, double *sfL, double *d1gfL, double *d1gfgt, double *gtend, int *i1, double *sa)
{
  double hhh, rnorm, rnormI, gfL_(), vsfL[npl], d1sfL[npl], d1sfgt[npl], sgc;
  int i, i0; 
      /* 
	 Make normalization coefficients for continuous solution. 
      */
      if(*gtend > 0.01&& (*g1cs)>fabs(*gL0)){
	hhh=*g1cs;
	rnormI=0.5*gfL_(i1, d1gfL, d1gfgt, gtend, &hhh, sa);
	rnorm=rnormI;
	hhh=-(*g1cs);
	rnorm+=0.5*gfL_(i1, d1gfL, d1gfgt, gtend, &hhh, sa);
	rnormI=rnorm/rnormI;
      }else {
	rnorm=0.5;
	rnormI=1;
      }
      for(i=0; i<npl; i++){
	sgc=-1.+2.*i/(npl-1);
	hhh=sgc;
	if((*g1cs) < fabs(*gL0)){
	if(sgc>(*g1cs)){
	  vsfL[i]=gfL_(i1, d1sfL+i, d1sfgt+i, gtend, &hhh, sa);
	  hhh=2.*(*g1cs)-sgc;
	  vsfL[i]-=0.5*gfL_(i1, sfL, d1gfgt, gtend, &hhh, sa);
	  d1sfL[i]+=0.5*(*sfL);
	  d1sfgt[i]-=0.5*(*d1gfgt);
	  d1sfL[i]*=rnormI;
	  d1sfgt[i]*=rnormI;
	  vsfL[i]*=rnormI;
	}
	else if(fabs(sgc) < (*g1cs)){
	  vsfL[i]=0.5*gfL_(i1, sfL, d1sfgt+i, gtend, &hhh, sa); 
	  hhh=-sgc;
	  vsfL[i]+=0.5*gfL_(i1, d1sfL+i, d1gfgt, gtend, &hhh, sa);
	  d1sfgt[i]+=(*d1gfgt);
	  d1sfgt[i]*=0.5;
	  d1sfL[i]-=(*sfL);
	  d1sfL[i]*=-0.5;
	} 
	else{
/*
  Counter-passing fast ions.
*/
	  hhh=sgc+2*(*g1cs);
	  vsfL[i]=rnormI*0.5*gfL_(i1, d1sfL+i, d1sfgt+i, gtend, &hhh, sa);
	  d1sfgt[i]*=0.5*rnormI;
	  d1sfL[i]*=0.5*rnormI;
	}}
	else{
	  vsfL[i]=gfL_(i1, d1sfL+i, d1sfgt+i, gtend, &hhh, sa);
	}
      }	
      hhh=0.;
      for(i=0; i<npl; i++){
	hhh+=vsfL[i];
      }
      hhh=1./hhh;
      hhh*=npl*0.5;
      for(i=0; i<npl; i++){
	vsfL[i]*=hhh;
	d1sfL[i]*=hhh;
	d1sfgt[i]*=hhh;
      }
/* 
   interpolate in time to isotropic from anisotropic solution 
*/
      i0=(npl-1)*((*gc+1.000001)*0.5)+1-1;
      i=(i0<1)?1:i0;
      i=(i>npl-3)?(npl-3):i; 
      hhh=0.5*(*gc+1.)*(npl-1.)-i;
      fun4_eq_(vsfL+i-1, &hhh, sfL);
      fun4_eq_(d1sfL+i-1, &hhh, d1gfL);
      fun4_eq_(d1sfgt+i-1, &hhh, d1gfgt);
      hhh=(*gtend)*(*gtend)*(*gtend);
      *sfL=(0.5*hhh+(*sfL))/(hhh+1);
      *d1gfL=(*d1gfL)/(hhh+1);
      *d1gfgt=((-(*sfL)+
	0.5)*3.*(*gtend)*(*gtend)+(*d1gfgt))/(hhh+1);
/*
      printf(" chis %d %g %g %g \n",i,*gc,-1.+i*2./(npl-1.),hhh);
    printf(" %d %g %g \n",i,sfL[i],sPl);
      sf[i]*=Nprt/50.;
*/
}
/* 
************************************************************ 
Summ over the Legendre polynomial for the solution at a time gtend.

It is based on the recurrence equation: 
l*P_l = x*(2l-1)*P_{l-1} - (l-1)*P_{l-2}

and its derivative:
l*P_l' = (2l-1)*(P_{l-1}+x*P_{l-1}') - (l-1)*P_{l-2}'

*/
double gfL_(int *N, double *d1gfL, double *d1gfgt, double *gt, double *chi, double *sa){
  int i1;
  double sfL;
  double sPl, sP1, sP0, d1sPl, d1sP1, d1sP0, hhh;
  sfL=0.;
  *d1gfL=0.;
  *d1gfgt=0.;
  for(i1=0; i1<(*N); i1++){
    if(i1==0){
      sPl=1;
      d1sPl=0;}
    else if(i1==1){
      sPl=(*chi);
      d1sPl=1.;}
    else {
      d1sPl=1./i1;
      sPl=((2*i1-1)*(*chi)*sP1-(i1-1)*sP0)*d1sPl;
      d1sPl*=((2*i1-1)*((*chi)*d1sP1+sP1)-(i1-1)*d1sP0);
    }
    sP0=sP1;
    sP1=sPl;
    d1sP0=d1sP1;
    d1sP1=d1sPl;
    hhh=sa[i1]*exp(-i1*(i1+1)*(*gt)*0.5);
    *d1gfL+=d1sPl*hhh;
    hhh*=sPl;
    sfL+=hhh;
    *d1gfgt+=-i1*(i1+1)*0.5*hhh;
    /*
    printf(" %d %g %g %g %g %g %g \n",i1,sfL,sPl,sP1,sP0,*gt,*chi);
    */
  }
  return sfL;
}
/* 
************************************************************ 
*/
double gf(double t,double x,double x0){
  return exp(-(x-x0)*(x-x0)/(2*t))/sqrt(t*6.283185307); 
}
/* 
************************************************************ 
*/
void sal_(int *N, double *sal, double *gt0, double *gL0){
  int i, j, Nint=20001;
  double chi, gdchi=2./(Nint-1), sPl, sP0, sP1;
  for(j=0; j<(*N); j++){sal[j]=0;}
  for(i=0; i<Nint; i++){
    chi=-1.+i*gdchi;
    for(j=0; j<(*N); j++){
      if(j==0){sPl=1;}
      else if(j==1){sPl=chi;}
      else {
	sPl=((2*j-1)*chi*sP1-(j-1)*sP0)/j;
      }
      sP0=sP1;
      sP1=sPl;
      sal[j]+=gf(*gt0, chi, *gL0)*sPl;
      /*
      sal[j]+=1*sPl;
      if(j==2){printf(" %g %g \n",chi,sPl);}
      */
    }
    /*
    printf(" %g ",chi);
    */
  }
  for(j=0; j<(*N); j++){
    sal[j]*=(j+0.5)*gdchi;
    /*
    printf(" %g \n",sal[j]);
    */
  }
}
