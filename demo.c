#define OTSU_NOMAIN
#include "model_data.c"

#define NOMAIN
#include "analyze.c"

int main (int argc, char *argv[]) {
  int bpp = 0, w = 0, h = 0;
  unsigned char* stb_rgb = stbi_load(argv[1],&w,&h,&bpp,3);
  // ... or use stbi_load_from_memory if jpg data is in ram

  float vec[36] = {};
  analyze(stb_rgb,w,h,vec);
  printf("score %f\n", score(vec, model));

  stbi_image_free(stb_rgb);
}


