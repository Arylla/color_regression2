#include <stdio.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>

#define STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_RESIZE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image.h"
#include "stb_image_write.h"
#include "stb_image_resize.h"

/* -------------------------------------------------------------------------- */

void mstd (int *vns, int L, int R, float p_mu, float p_std, float p_n, 
    float *mu_out, float *std_out) {
  float totv  = 0 + p_mu*p_n;
  float totv2 = 0 + (p_std*p_std + p_mu*p_mu)*p_n;
  float totn = 0 + p_n;
  for (int v=L;v<R;v++) {
    int n = vns[v];
    totv += v*n;
    totv2 += v*v*n;
    totn += n;
  }
  float mu = totv/((float) totn);
  float va = totv2/((float) totn) - mu*mu;
  *mu_out = mu;
  *std_out = pow(va,0.5);
}

float gauss (float x, float m, float s) {
  float z = (x-m)/s;
  float e = exp(-z*z/2.0);
  e /= pow(2*3.141592653589793*s,0.5);
  if (e < 1e-16) e = 1e-16;
  return e;
}

void opp5 (float A, float B, int is_BW, float *out) {
  float K = 1.0;
  if (is_BW) K = 2.0;
  float mAB = A < B ? A : B;
  float mid = K*mAB;
  A -= mAB;
  B -= mAB;
  mAB = mid;
  
  float WAm = A < mAB ? A : mAB;
  float Am = 2*WAm;
  A -= WAm;
  mAB -= WAm;
  
  float WBm = B < mAB ? B : mAB;
  float Bm = 2*WBm;
  B -= WBm;
  mAB -= WBm;
 
  out[0] = A;
  out[1] = Am;
  out[2] = mAB;
  out[3] = Bm;
  out[4] = B;
}

void prism (float R, float G, float B, float *out) {
  for (int i=0;i<18;i++) out[i] = 0.0;

  float WH = R < G ? R : G;
  WH = WH < B ? WH : B;
  float BK = R > G ? R : G;
  BK = BK > B ? BK : B;
  BK = 1.0 - BK;
  opp5(WH,BK,1,out);

  R -= WH;
  G -= WH;
  B -= WH;

  if (B < 1e-6) opp5(R,G,0,out+5);
  else if (R < 1e-6) opp5(G,B,0,out+9);
  else {
    opp5(B,R,0,out+13);
    out[5] = out[17];
    out[17] = 0.0;
  } 
}

void otsu (unsigned char *rgb, int w, int h, float* out) {
  // prepare the gamma map
  float gamma_expand[256] = {}; int gamma_compress[256] = {};
  for (int i=0;i<256;i++) {
    gamma_expand[i] = pow(((float)i)/255.0,2.2);
    gamma_compress[i] = 0; //(int) (255.0*pow(((int)i)/255.0,1.0/2.2) + 0.5);
  }

  // calculate the luma map
  int luma_counts[1001] = {};
  for (int i=0;i<h;i++) {
    for (int j=0;j<w;j++) {
      unsigned char *p = rgb + 3*(i*w + j);
      float Rf = gamma_expand[p[0]];
      float Gf = gamma_expand[p[1]];
      float Bf = gamma_expand[p[2]];
      int luma = (int) ((100*(2*Rf + 7*Gf + 1*Bf))+0.5);
      luma_counts[luma] += 1;
    }
  } 

  // calculate global statistics
  float global_m = 0, global_s = 0;
  mstd(luma_counts,0,1001,0.0,1.0,0, &global_m, &global_s);

  // find best split 
  float best_score = -1;
  float best_mL = -1.0, best_mR = -1.0, best_sL = -1.0, best_sR = -1.0;
  int nL = 0, nR = 0, best_M = -1, best_nL = -1, best_nR = -1;
  for (int i=0;i<1001;i++) nR += luma_counts[i];
  for (int M=1;M<1001;M++) {
    float mL = 0, sL = 0, mR = 0, sR = 0;
    nL += luma_counts[M-1];
    nR -= luma_counts[M-1];
    mstd(luma_counts,0,M,global_m,global_s,100.0, &mL, &sL);
    mstd(luma_counts,M,1001,global_m,global_s,100.0, &mR, &sR);
    float score = nL*sL*sL + nR*sR*sR;
    if (best_score < 0 || score < best_score) {
      best_score = score;
      best_M = M;
      best_mL = mL;
      best_mR = mR;
      best_sL = sL;
      best_sR = sR;
      best_nL = nL;
      best_nR = nR;
    } 
  }

  float brightness = best_mL/1000.0, contrast = (best_mR-best_mL)/1000.0;
  float pri = best_nL/((float) best_nL + best_nR);

  float WLS[1001] = {};
  float WRS[1001] = {};
  for (int i=0;i<1001;i++) {
    float pL = pri*gauss(i,best_mL,best_sL);
    float pR = (1-pri)*gauss(i,best_mR,best_sR);
    float WR = 0.0, WL = 0.0;
    if (pR > 9*pL) WR = (pR-9*pL)/(pR+9*pL);
    if (pL > 9*pR) WL = (pL-9*pR)/(pL+9*pR);
    WLS[i] = WL;
    WRS[i] = WR;
  }  

  float WpriL = 0.0, WpriR = 0.0;
  float prism_L[18] = {}, prism_R[18] = {};
  for (int oi=0;oi<w*h;oi++) {
    unsigned char *px = rgb + 3*oi;
    float Rf = gamma_expand[px[0]];
    float Gf = gamma_expand[px[1]];
    float Bf = gamma_expand[px[2]];   
    int luma = (int) ((100*(2*Rf + 7*Gf + 1*Bf))+0.5);
    float WL = WLS[luma];
    float WR = WRS[luma];
    float curr_prism[18] = {};
    prism(Rf,Gf,Bf,curr_prism);
    for (int j=0;j<18; j++) {
      prism_L[j] += WL*curr_prism[j];
      prism_R[j] += WR*curr_prism[j];
    }
    WpriL += WL;
    WpriR += WR;
  }

  WpriL = WpriL > 1e-4 ? WpriL : 1e-4;
  WpriR = WpriR > 1e-4 ? WpriR : 1e-4;

  for (int i=0;i<18;i++) {
    prism_L[i] /= WpriL;
    prism_R[i] /= WpriR;
  }
  prism_L[17] = brightness;
  prism_R[17] = contrast;

  for (int i=0;i<18;i++) {
    out[i] = prism_L[i];
    out[i+18] = prism_R[i];
  }
}

float analyze (uint8_t *pix, int w, int h, float *out) {
  float rat = 1.0;
  if (w > h) rat = w > 384 ? (((float) 384) / w) : 1.0;
  if (h > w) rat = h > 384 ? (((float) 384) / h) : 1.0;
  int w2 = (int) (rat*w + 0.5);
  int h2 = (int) (rat*h + 0.5);
  unsigned char* rgb384 = malloc(3*384*384);
  stbir_resize_uint8_srgb 
    (pix, w,h,0,rgb384,w2,h2,0,3,STBIR_ALPHA_CHANNEL_NONE, 0);

  otsu(rgb384,w2,h2,out);

  free(rgb384);
  return 0.0;
}

int bucketize (float v, float *buckL, float *buckR) {
  int i = 61/2;
  float L = buckL[i], R = buckR[i];
  if (L <= v && v <= R) return i;
  int D = v < L ? -1 : 1;
  while (0 < i && i < 60) {
    i += D;
    L = buckL[i];
    R = buckR[i];
    if (L <= v && v <= R) return i;
  } 
  return i;
}

// ovec is a 36-slot vector from otsu analysis
float score (float *ovec, float* model) {
  int i=0;
  float tot = 0.0;
  for (int K0=0;K0<36;K0++) {
    float v0 = ovec[K0];
    for (int K1=K0+1;K1<36;K1++) {
      float v1 = ovec[K1];
      float *off = model + 65*61*i;
      int b0 = bucketize(v0, off , off + 61);
      int b1 = bucketize(v1, off + 61*2, off + 61*3);
      off += 61*4;
      tot += off[61*b0 + b1];
      //printf("%i:%f:%i:%i:%f:%i:%f\n",K0,v0,b0,K1,v1,b1,off[61*b0+b1]);
      i++;
    }
  }  
  return tot;
}

#ifndef NOMAIN

int main (int argc, char *argv[]) {
//  FILE *f = fopen(argv[1],"rb");

  const char* labs[18] = {
    "WH", "LG", "MG", "DG", "BK", "R", "O", "Y", "SG", "G", "TU", "CY", "AZ",
    "B", "PU","MA","PI","BC"
  };

  // load the image, copy it to local memory, convert to linear color
  int bpp = 0, w = 0, h = 0;
  unsigned char* stb_rgb = stbi_load(argv[1],&w,&h,&bpp,3);
  float vec[36] = {};
  analyze(stb_rgb,w,h,vec); 

  for (int i=0;i<18;i++) printf("%4s %6.4f %6.4f\n", labs[i], vec[i], vec[i+18]);

  stbi_image_free(stb_rgb); 
}

#endif
