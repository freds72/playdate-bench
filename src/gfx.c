#include <pd_api.h>
#include <float.h>
#include "simd.h"
#include "gfx.h"

// float32 display ptr width
#define LCD_ROWSIZE32 (LCD_ROWSIZE/4)

static PlaydateAPI* pd = NULL;

void gfx_init(PlaydateAPI* playdate) {
    pd = playdate;
}


#if TARGET_PLAYDATE
static __attribute__((always_inline))
#else
static __forceinline
#endif
inline void _drawMaskPattern(uint32_t* p, uint32_t mask, uint32_t color)
{
    *p = (*p & ~mask) | (color & mask);
}

#if TARGET_PLAYDATE
static __attribute__((always_inline))
#else
static __forceinline
#endif
inline void _drawMaskPatternOpaque(uint32_t* p, uint32_t color)
{
    *p = color;
}

static void drawFragment(uint32_t* row, int x1, int x2, uint32_t color)
{
    if (x2 < 0 || x1 >= LCD_COLUMNS)
        return;

    if (x1 < 0)
        x1 = 0;

    if (x2 > LCD_COLUMNS)
        x2 = LCD_COLUMNS;

    if (x1 > x2)
        return;

    // Operate on 32 bits at a time

    const int startbit = x1 & 31;
    const uint32_t startmask = swap((1 << (32 - startbit)) - 1);
    const int endbit = x2 & 31;
    const uint32_t endmask = swap(((1 << endbit) - 1) << (32 - endbit));

    const int col = x1 >> 5;
    uint32_t* p = row + col;

    if (col == x2 >> 5)
    {
        uint32_t mask = 0;

        if (startbit > 0 && endbit > 0)
            mask = startmask & endmask;
        else if (startbit > 0)
            mask = startmask;
        else if (endbit > 0)
            mask = endmask;

        _drawMaskPattern(p, mask, color);
    }
    else
    {
        int x = x1;

        if (startbit > 0)
        {
            _drawMaskPattern(p++, startmask, color);
            x += (32 - startbit);
        }

        x2 -= 32;
        while (x <= x2)
        {
            _drawMaskPatternOpaque(p++, color);
            x += 32;
        }

        if (endbit > 0)
            _drawMaskPattern(p, endmask, color);
    }
}

// x1/x2 are fixed 16:16!!
static void drawTextureFragment(uint8_t* row, int x1, int x2, int lu, int lv, int lw, int ru, int rv, int rw, uint8_t* texture)
{
    if (x2 < 0 || x1 >= LCD_COLUMNS<<16)
        return;

    int dx = x2 - x1;
    float frac = 1.0f - ((x1 & 0xffff) / (float)0xffff);
    if (dx < 0xffff) return;
    // uvw source is fixed point already
    dx >>= 16;
    int du = (ru - lu) / dx;
    int dv = (rv - lv) / dx;
    int dw = (rw - lw) / dx;
    x1 >>= 16;
    x2 >>= 16;
    
    if (x1 < 0) {
        lu -= x1 * du;
        lv -= x1 * dv;
        lw -= x1 * dw;
        x1 = 0;
        frac = 0.f;
    }
    // sub-pix shift
    lu += (int)(frac * du);
    lv += (int)(frac * dv);
    lw += (int)(frac * dw);
    
    if (x2 > LCD_COLUMNS)
        x2 = LCD_COLUMNS;

    if (x1 > x2)
        return;

    // Operate on 8 bits at a time

    const int startbit = x1 & 7;
    const uint8_t startmask = (1 << (8 - startbit)) - 1;
    const int endbit = x2 & 7;
    const uint8_t endmask = ((1 << endbit) - 1) << (8 - endbit);
    
    const int col = x1 >> 3;
    uint8_t* p = row + col;
    
    if (col == x2 >> 3)
    {
        uint8_t mask = 0;

        if (startbit > 0 && endbit > 0)
            mask = startmask & endmask;
        else if (startbit > 0)
            mask = startmask;
        else if (endbit > 0)
            mask = endmask;
        
        // build pixels
        uint8_t src = 0;
        for (int x = x1; x < x2; x++) {
            int u = (lu / lw)&31, v = (lv / lw)&31;
            // texture encoded as 1 pixel per byte
            src |= texture[u + v * 32] >> (x&7);
            // src |= 0x80 >> (x & 7);
            lu += du;
            lv += dv;
            lw += dw;
        }
        *p = (*p & ~mask) | src;
    }
    else
    {
        int x = x1;

        if (startbit > 0)
        {
            uint8_t src = 0;
            for (int i = startbit; i < 8; i++) {
                int u = (lu / lw) &31, v = (lv / lw) &31;
                // texture encoded as 1 pixel per byte
                src |= texture[u + v * 32] >> i;
                // src |= 0x80 >> i;
                lu += du;
                lv += dv;
                lw += dw;
            }
            *p = (*p & ~startmask) | src;
            p++;
            x += (8 - startbit);
        }

        x2 -= 8;
        while (x <= x2)
        {
            uint8_t src = 0;
            for (int i = 0; i < 8; i++) {
                int u = (lu / lw) &31, v = (lv / lw) &31;
                // texture encoded as 1 pixel per byte
                src |= texture[u + v * 32] >> i;
                // src |= 0x80 >> i;
                lu += du;
                lv += dv;
                lw += dw;
            }

            *(p++) = src;
            x += 8;
        }

        if (endbit > 0) {
            uint8_t src = 0;
            for (int i = 0; i < endbit; i++) {
                int u = (lu / lw) &31, v = (lv / lw) & 31;
                // texture encoded as 1 pixel per byte
                src |= texture[u + v * 32] >> i;
                // src |= 0x80 >> i;
                lu += du;
                lv += dv;
                lw += dw;
            }
            *p = (*p & ~endmask) | src;
        }
    }
}

void polyfill(const Point3duv* verts, const int n, uint32_t* dither, uint32_t* bitmap) {
	float miny = FLT_MAX, maxy = -FLT_MAX;
	int mini = -1;
	// find extent
	for (int i = 0; i < n; ++i) {
		float y = verts[i].y;
		if (y < miny) miny = y, mini = i;
		if (y > maxy) maxy = y;
	}
    // out of screen?
    if (miny > LCD_ROWS || maxy < 0.f) {
        return;
    }

	// data for left& right edges :
	int lj = mini, rj = mini;
    int ly = -1, ry = -1;
    int lx = 0, ldx = 0, rx = 0, rdx = 0;
    int ystart = (int)ceilf(miny), yend = (int)ceilf(maxy);
    if (yend > LCD_ROWS) yend = LCD_ROWS;
    if (ystart < 0) ystart = 0;
    bitmap += ystart * LCD_ROWSIZE32;
    for (int y = ystart; y < yend; y++, bitmap += LCD_ROWSIZE32, lx += ldx, rx += rdx) {
        // maybe update to next vert
        while (ly < y) {
            const Point3duv* p0 = &verts[lj];
            // lj = (lj + 1) % n;
            lj++;
            if (lj >= n) lj = 0;
            const Point3duv* p1 = &verts[lj];
            const float y0 = p0->y, y1 = p1->y;
            ly = (int)y1;
            ldx = __TOFIXED16((p1->x - p0->x) / (y1 - y0));
            //sub - pixel correction
            lx = __TOFIXED16(p0->x) + (int)((y - y0) * ldx);
        }
        while (ry < y) {
            const Point3duv* p0 = &verts[rj];
            // rj = (rj + n - 1) % n;
            rj--;
            if (rj < 0) rj = n - 1;
            const Point3duv* p1 = &verts[rj];
            const float y0 = p0->y, y1 = p1->y;
            ry = (int)y1;
            rdx = __TOFIXED16((p1->x - p0->x) / (y1 - y0));
            //sub - pixel correction
            rx = __TOFIXED16(p0->x) + (int)((y - y0) * rdx);
        }

        drawFragment(bitmap, lx>>16, rx>>16, dither[y&31]);
    } 
}

// perspective correct texturing
// z contains dither color
void texfill(const Point3duv* verts, const int n, uint8_t* texture, uint8_t* bitmap) {

    float miny = FLT_MAX, maxy = -FLT_MAX;
    int mini = -1;
    // find extent
    for (int i = 0; i < n; ++i) {
        float y = verts[i].y;
        if (y < miny) miny = y, mini = i;
        if (y > maxy) maxy = y;
    }
    // out of screen?
    if (miny > LCD_ROWS || maxy < 0) {
        return;
    }

    // data for left& right edges :
    int lj = mini, rj = mini;
    int ly = -1, ry = -1;
    int lx = 0, ldx = 0, rx = 0, rdx = 0;
    int lu = 0, ldu = 0, ru = 0, rdu = 0;   
    int lv = 0, ldv = 0, rv = 0, rdv = 0;
    int lw = 0, ldw = 0, rw = 0, rdw = 0;
    int ystart = (int)ceilf(miny), yend = (int)ceilf(maxy);
    if (yend > LCD_ROWS) yend = LCD_ROWS;
    if (ystart < 0) ystart = 0;
    bitmap += ystart * LCD_ROWSIZE;
    for (int y = ystart; y < yend; y++, bitmap += LCD_ROWSIZE, lx += ldx, rx += rdx, lu += ldu, ru += rdu, lv += ldv, rv += rdv, lw += ldw, rw += rdw) {
        // maybe update to next vert
        while (ly < y) {
            const Point3duv* p0 = &verts[lj];
            lj++;
            if (lj >= n) lj = 0;
            const Point3duv* p1 = &verts[lj];
            const float y0 = p0->y, y1 = p1->y;
            const float w0 = p0->z, w1 = p1->z;
            const float u0 = p0->u * w0, v0 = p0->v * w0;
            const float dy = y1 - y0;
            ly = (int)y1;
            ldx = __TOFIXED16((p1->x - p0->x) / dy);
            ldw = __TOFIXED16((w1 - w0) / dy);
            ldu = __TOFIXED16((p1->u * w1 - u0) / dy);
            ldv = __TOFIXED16((p1->v * w1 - v0) / dy);
            //sub - pixel correction
            const float cy = y - y0;
            lx = __TOFIXED16(p0->x) + (int)(cy * ldx);
            lw = __TOFIXED16(w0) + (int)(cy * ldw);
            lu = __TOFIXED16(u0) + (int)(cy * ldu);
            lv = __TOFIXED16(v0) + (int)(cy * ldv);
        }
        while (ry < y) {
            const Point3duv* p0 = &verts[rj];
            rj--;
            if (rj < 0) rj = n - 1;
            const Point3duv* p1 = &verts[rj];
            const float y0 = p0->y, y1 = p1->y;
            const float w0 = p0->z, w1 = p1->z;
            const float u0 = p0->u * w0, v0 = p0->v * w0;
            const float dy = y1 - y0;
            ry = (int)y1;
            rdx = __TOFIXED16((p1->x - p0->x) / dy);
            rdw = __TOFIXED16((w1 - w0) / dy);
            rdu = __TOFIXED16((p1->u * w1 - u0) / dy);
            rdv = __TOFIXED16((p1->v * w1 - v0) / dy);
            //sub - pixel correction
            const float cy = y - y0;
            rx = __TOFIXED16(p0->x) + (int)(cy * rdx);
            rw = __TOFIXED16(w0) + (int)(cy * rdw);
            ru = __TOFIXED16(u0) + (int)(cy * rdu);
            rv = __TOFIXED16(v0) + (int)(cy * rdv);
        }
        drawTextureFragment(bitmap, lx, rx, lu, lv, lw, ru, rv, rw, texture);
    }
}
