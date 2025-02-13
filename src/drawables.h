#ifndef _drawables_h
#define _drawables_h

#include <pd_api.h>
#include "3dmath.h"

#define MAX_DRAWABLES 2048

typedef struct {
    // original flags
    int flags;
    // number of points
    int n;
    // texture
    Texture* texture;
    // clipped points in camera space
    Point3duv pts[5];
} DrawableFace;

typedef struct {    
    union {
        struct {
            uint16_t alpha;
            uint16_t color;
        };
        uint32_t material;
    };
    float radius;
    float angle;
    Point3d pos;
} DrawableParticle;

struct Drawable_s;
typedef void(*draw_drawable)(struct Drawable_s* drawable, uint8_t* bitmap);

// generic drawable thingy
// cache-friendlyness???
typedef struct Drawable_s {
    float key;
    draw_drawable draw;
    union {
        DrawableFace face;
        DrawableParticle particle;
    };
} Drawable;

void drawables_init(PlaydateAPI* playdate);
void reset_drawables();
Drawable* pop_drawable(const float sortkey);
void draw_drawables(uint8_t* bitmap);

#endif
