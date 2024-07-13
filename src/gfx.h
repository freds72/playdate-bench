#ifndef lib3d_gfx_h
#define lib3d_gfx_h

#include "3dmath.h"

void gfx_init(PlaydateAPI* playdate);
void polyfill(const Point3duv* verts, const int n, uint32_t* dither, uint32_t* bitmap);
void texfill(const Point3duv* verts, const int n, Texture* texture, uint8_t* bitmap);
void texfill_baseline(const Point3duv* verts, const int n, Texture* texture, uint8_t* bitmap);

#endif
