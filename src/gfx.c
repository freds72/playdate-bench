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
static void drawTextureFragment(uint8_t* row, float lx, float rx, float lu, float lv, float lw, float ru, float rv, float rw, uint8_t* texture, const int tw, const int tmask)
{
    int x1 = (int)lx, x2 = (int)rx;
    if (x2 < 0 || x1 >= LCD_COLUMNS)
        return;

    if (x1 > x2)
        return;

    float frac = x1 - lx;
    float dx = rx - lx;

    if (x2 > LCD_COLUMNS)
        x2 = LCD_COLUMNS;

    if ((int)dx == 0) return;

    // 
    const float du = (ru - lu) / dx;
    const float dv = (rv - lv) / dx;
    const float dw = (rw - lw) / dx;
    
    if (x1 < 0) {
        lu -= x1 * du;
        lv -= x1 * dv;
        lw -= x1 * dw;
        x1 = 0;
        frac = 0.f;
    }
    // sub-pix shift
    lu += frac * du;
    lv += frac * dv;
    lw += frac * dw;
    
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
            int u = (int)(lu / lw)&tmask, v = (int)(lv / lw)& tmask;
            // texture encoded as 1 pixel per byte
            src |= texture[u + v * tw] >> (x&7);
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
                int u0 = (int)(lu / lw) & tmask, v0 = (int)(lv / lw) & tmask;
                // texture encoded as 1 pixel per byte
                src |= texture[u0 + v0 * tw] >> i;
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
        if (x <= x2) {
            // starting point
            float u0 = (lu / lw), v0 = (lv / lw);
            const float du_strip = 8 * du, dv_strip = 8 * dv, dw_strip = 8 * dw;
            while (x <= x2)
            {
                // next
                lu += du_strip;
                lv += dv_strip;
                lw += dw_strip;
                const float u1 = lu / lw, v1 = lv / lw;
                const float ddu = (u1 - u0) / 8.0f, ddv = (v1 - v0) / 8.0f;
                uint8_t src = 0;
                for (int i = 0; i < 8; i++, u0 += ddu, v0 += ddv) {
                    // texture encoded as 1 pixel per byte
                    src |= texture[(((int)u0) & tmask) + ((((int)v0)) & tmask) * tw]>>i;
                    // src |= 0x80 >> i;
                }

                *(p++) = src;
                x += 8;
                u0 = u1;
                v0 = v1;
            }
        }
        if (endbit > 0) {
            uint8_t src = 0;
            for (int i = 0; i < endbit; i++) {
                int u0 = (int)(lu / lw) & tmask, v0 = (int)(lv / lw) & tmask;
                // texture encoded as 1 pixel per byte
                src |= texture[u0 + v0 * tw] >> i;
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
void texfill(const Point3duv* verts, const int n, Texture* texture, uint8_t* bitmap) {

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
    float lx = 0.f, ldx = 0.f, rx = 0.f, rdx = 0.f;
    float lu = 0.f, ldu = 0.f, ru = 0.f, rdu = 0.f; 
    float lv = 0.f, ldv = 0.f, rv = 0.f, rdv = 0.f;
    float lw = 0.f, ldw = 0.f, rw = 0.f, rdw = 0.f;
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
            ldx = (p1->x - p0->x) / dy;
            ldw = (w1 - w0) / dy;
            ldu = (p1->u * w1 - u0) / dy;
            ldv = (p1->v * w1 - v0) / dy;
            //sub - pixel correction
            const float cy = y - y0;
            lx = p0->x + (cy * ldx);
            lw = w0 + (cy * ldw);
            lu = u0 + (cy * ldu);
            lv = v0 + (cy * ldv);
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
            rdx = (p1->x - p0->x) / dy;
            rdw = (w1 - w0) / dy;
            rdu = (p1->u * w1 - u0) / dy;
            rdv = (p1->v * w1 - v0) / dy;
            //sub - pixel correction
            const float cy = y - y0;
            rx = p0->x + (cy * rdx);
            rw = w0 + (cy * rdw);
            ru = u0 + (cy * rdu);
            rv = v0 + (cy * rdv);
        }
        drawTextureFragment(bitmap, lx, rx, lu, lv, lw, ru, rv, rw, texture->data, texture->size, texture->size-1);
        // drawFragment((uint32_t*)bitmap, lx >> 16, rx >> 16, 0x0);
    }
}

static int32_t fixed_div_16(int32_t x, int32_t y)
{
    return ((int64_t)x * (1 << FIXED16_SHIFT)) / y;
}

// perspective correct texturing - base version
// x1/x2 are fixed 16:16!!
// tw: texture width
// tsize: texture modulo mask
static void drawTextureFragment_baseline(uint8_t* row, int x1, int x2, int lu, int lv, int lw, int ru, int rv, int rw, uint8_t* texture, const int tw, const int tmask)
{
    if (x2 < 0 || x1 >= LCD_COLUMNS << 16)
        return;

    if (x1 > x2)
        return;

    float frac = 1.0f - (x1 & 0xffff) / (float)0xffff;
    // keep fixed point precision
    const float dx = x2 - x1;
    // convert to screen units
    x1 >>= 16;
    x2 >>= 16;

    if (x2 > LCD_COLUMNS)
        x2 = LCD_COLUMNS;

    // int dx = x2 - x1;
    if (x2-x1 == 0) return;

    // delta have a precision of 1/16
    const int du = __TOFIXED16((ru - lu) / dx);
    const int dv = __TOFIXED16((rv - lv) / dx);
    const int dw = __TOFIXED16((rw - lw) / dx);

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
            int u = (lu / lw) & tmask, v = (lv / lw) & tmask;
            // texture encoded as 1 pixel per byte
            src |= texture[u + v * tw] >> (x & 7);
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
                int u = (lu / lw) & tmask, v = (lv / lw) & tmask;
                // texture encoded as 1 pixel per byte
                src |= texture[u + v * tw] >> i;
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

        if (x <= x2) {
            // starting point (with additional precision)
            int u0 = ((64 * lu) / lw), v0 = ((64 * lv) / lw) ;
            const int du_strip = 8 * du, dv_strip = 8 * dv, dw_strip = 8 * dw;
            while (x <= x2)
            {
                // next
                lu += du_strip;
                lv += dv_strip;
                lw += dw_strip;
                // end point (with additional precision)
                const int u1 = ((64 * lu) / lw)  , v1 = ((64 * lv) / lw) ;
                const int ddu = (u1 - u0) / 8, ddv = (v1 - v0) / 8;
                uint8_t src = 0;
                for (int i = 0; i < 8; i++, u0 += ddu, v0 += ddv) {
                    int u = (u0>>6) & tmask, v = (v0>>6) & tmask;
                    // texture encoded as 1 pixel per byte
                    src |= texture[u + v * tw] >> i;
                    // src |= 0x80 >> i;
                }

                *(p++) = src;
                x += 8;
                u0 = u1;
                v0 = v1;
            }
        }

        if (endbit > 0) {
            uint8_t src = 0;
            for (int i = 0; i < endbit; i++) {
                int u = (lu / lw) & tmask, v = (lv / lw) & tmask;
                // texture encoded as 1 pixel per byte
                src |= texture[u + v * tw] >> i;
                // src |= 0x80 >> i;
                lu += du;
                lv += dv;
                lw += dw;
            }
            *p = (*p & ~endmask) | src;
        }
    }
}

// z contains dither color
void texfill_baseline(const Point3duv* verts, const int n, Texture* texture, uint8_t* bitmap) {

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
        drawTextureFragment_baseline(bitmap, lx, rx, lu, lv, lw, ru, rv, rw, texture->data, texture->size, texture->size-1);
        // drawFragment((uint32_t*)bitmap, lx >> 16, rx >> 16, 0x0);
    }
}


void polyfill_fixed(const Point3duv* verts, const int n, uint32_t* bitmap) {
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
    Fixed x = { 0 }, dx = { 0 };
    int ystart = (int)ceilf(miny), yend = (int)ceilf(maxy);
    if (yend > LCD_ROWS) yend = LCD_ROWS;
    if (ystart < 0) ystart = 0;
    bitmap += ystart * LCD_ROWSIZE32;
    for (int y = ystart; y < yend; y++, bitmap += LCD_ROWSIZE32, x.i32 += dx.i32) {
        // maybe update to next vert
        while (ly < y) {
            const Point3duv* p0 = &verts[lj];
            // lj = (lj + 1) % n;
            lj++;
            if (lj >= n) lj = 0;
            const Point3duv* p1 = &verts[lj];
            const float y0 = p0->y, y1 = p1->y;
            ly = (int)y1;
            dx.q0 = __TOFIXED8((p1->x - p0->x) / (y1 - y0));
            //sub - pixel correction
            x.q0 = __TOFIXED8(p0->x) + (int16_t)((y - y0) * dx.q0);
        }
        while (ry < y) {
            const Point3duv* p0 = &verts[rj];
            // rj = (rj + n - 1) % n;
            rj--;
            if (rj < 0) rj = n - 1;
            const Point3duv* p1 = &verts[rj];
            const float y0 = p0->y, y1 = p1->y;
            ry = (int)y1;
            dx.q1 = __TOFIXED8((p1->x - p0->x) / (y1 - y0));
            //sub - pixel correction
            x.q1 = __TOFIXED8(p0->x) + (int16_t)((y - y0) * dx.q1);
        }

        drawFragment(bitmap, x.q0 >> FIXED8_SHIFT, x.q1 >> FIXED8_SHIFT, 0);
    }
}

// perspective correct texturing - fixed 12:4 version
// tw: texture width
// tsize: texture modulo mask
static void drawTextureFragment_fixed(uint8_t* row, Fixed x, Fixed luv, Fixed ruv, Fixed w, uint8_t* texture, const int tw, const int tmask)
{
    if (x.q0 < 0 || x.q1 >= __TOFIXED8(LCD_COLUMNS))
        return;

    if (x.q0 > x.q1)
        return;

    // keep dx as fixed
    int16_t dx = x.q1 - x.q0;
    // convert to screen units
    int x1 = x.q0 >> FIXED8_SHIFT;
    int x2 = x.q1 >> FIXED8_SHIFT;

    if (x2 > LCD_COLUMNS)
        x2 = LCD_COLUMNS;

    // int dx = x2 - x1;
    if (x2 - x1 == 0) return;
    // uvw source is fixed point already
    const int16_t dw = ((1<<FIXED8_SHIFT)*(w.q1 - w.q0)) / dx;
    const Fixed duv = {
        .q0 = ((1<<FIXED8_SHIFT)*(ruv.q0 - luv.q0)) / dx,
        .q1 = ((1<<FIXED8_SHIFT)*(ruv.q1 - luv.q1)) / dx
    };

    if (x1 < 0) {
        luv.q0 -= x1 * duv.q0;
        luv.q1 -= x1 * duv.q1;
        w.q0 -= x1 * dw;
        x1 = 0;
    }
    // sub-pix shift
    // lu += (int)(frac * du);
    // lv += (int)(frac * dv);
    // lw += (int)(frac * dw);

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
            int16_t u = (luv.q0 / w.q0) & tmask, v = (luv.q1 / w.q0) & tmask;
            // texture encoded as 1 pixel per byte
            src |= texture[u + v * tw] >> (x & 7);
            // src |= 0x80 >> (x & 7);
            luv.i32 += duv.i32;
            w.q0 += dw;
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
                int u = (luv.q0 / w.q0) & tmask, v = (luv.q1 / w.q0) & tmask;
                // texture encoded as 1 pixel per byte
                src |= texture[u + v * tw] >> i;
                // src |= 0x80 >> i;
                luv.i32 += duv.i32;
                w.q0 += dw;
            }
            *p = (*p & ~startmask) | src;
            p++;
            x += (8 - startbit);
        }

        x2 -= 8;

        Fixed uv0 = {
            .q0 = ((1<<FIXED8_SHIFT) * luv.q0) / w.q0,
            .q1 = ((1<<FIXED8_SHIFT) * luv.q1) / w.q0
        };
        const Fixed duv_stride = {
        .q0 = 8 * duv.q0,
        .q1 = 8 * duv.q1
        };
        const int16_t dw_stride = 8 * dw;
        while (x <= x2)
        {
            luv.i32 += duv_stride.i32;
            w.q0 += dw_stride;
            const Fixed uv1 = {
                .q0 = ((1<<FIXED8_SHIFT) * luv.q0) / w.q0,
                .q1 = ((1<<FIXED8_SHIFT) * luv.q1) / w.q0
            };
            const Fixed dduv = {
                .q0 = (uv1.q0 - uv0.q0) / 8,
                .q1 = (uv1.q1 - uv0.q1) / 8
            };
            uint8_t src = 0;
            for (int i = 0; i < 8; i++, uv0.i32 += dduv.i32) {
                int16_t u = (uv0.q0>>FIXED8_SHIFT) & tmask, v = (uv0.q1>> FIXED8_SHIFT) & tmask;
                // texture encoded as 1 pixel per byte
                src |= texture[u + v * tw] >> i;
            }

            *(p++) = src;
            x += 8;

            uv0 = uv1;
        }

        if (endbit > 0) {
            uint8_t src = 0;
            for (int i = 0; i < endbit; i++) {
                int u = (luv.q0 / w.q0) & tmask, v = (luv.q1 / w.q0) & tmask;
                // texture encoded as 1 pixel per byte
                src |= texture[u + v * tw] >> i;
                // src |= 0x80 >> i;
                luv.i32 += duv.i32;
                w.q0 += dw;
            }
            *p = (*p & ~endmask) | src;
        }
    }
}

void texfill_fixed(const Point3duv* verts, const int n, Texture* texture, uint8_t* bitmap) {

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
    Fixed x = { 0 }, dx = { 0 };
    Fixed luv = { 0 }, lduv = { 0 };
    Fixed ruv = { 0 }, rduv = { 0 };
    Fixed w = { 0 }, dw = { 0 };
    int ystart = (int)ceilf(miny), yend = (int)ceilf(maxy);
    if (yend > LCD_ROWS) yend = LCD_ROWS;
    if (ystart < 0) ystart = 0;
    bitmap += ystart * LCD_ROWSIZE;
    for (int y = ystart; y < yend; y++, bitmap += LCD_ROWSIZE, x.i32 += dx.i32, luv.i32 += lduv.i32, ruv.i32 += rduv.i32, w.i32 += dw.i32) {
        // maybe update to next vert
        while (ly < y) {
            const Point3duv* p0 = &verts[lj];
            lj++;
            if (lj >= n) lj = 0;
            const Point3duv* p1 = &verts[lj];
            const float y0 = p0->y, y1 = p1->y;
            const float w0 = 8 * p0->z, w1 = 8 * p1->z;
            const float u0 = p0->u * w0, v0 = p0->v * w0;
            const float dy = y1 - y0;
            ly = (int)y1;
            dx.q0 = __TOFIXED8((p1->x - p0->x) / dy);
            dw.q0 = __TOFIXED8((w1 - w0) / dy);
            lduv.q0 = __TOFIXED8((p1->u * w1 - u0) / dy);
            lduv.q1 = __TOFIXED8((p1->v * w1 - v0) / dy);
            //sub - pixel correction
            const float cy = y - y0;
            x.q0 = __TOFIXED8(p0->x) + (int16_t)(cy * dx.q0);
            w.q0 = __TOFIXED8(w0) + (int16_t)(cy * dw.q0);
            luv.q0 = __TOFIXED8(u0) + (int16_t)(cy * lduv.q0);
            luv.q1 = __TOFIXED8(v0) + (int16_t)(cy * lduv.q1);
        }
        while (ry < y) {
            const Point3duv* p0 = &verts[rj];
            rj--;
            if (rj < 0) rj = n - 1;
            const Point3duv* p1 = &verts[rj];
            const float y0 = p0->y, y1 = p1->y;
            const float w0 = 8 * p0->z, w1 = 8 * p1->z;
            const float u0 = p0->u * w0, v0 = p0->v * w0;
            const float dy = y1 - y0;
            ry = (int)y1;
            dx.q1 = __TOFIXED8((p1->x - p0->x) / dy);
            dw.q1 = __TOFIXED8((w1 - w0) / dy);
            rduv.q0 = __TOFIXED8((p1->u * w1 - u0) / dy);
            rduv.q1 = __TOFIXED8((p1->v * w1 - v0) / dy);
            //sub - pixel correction
            const float cy = y - y0;
            x.q1 = __TOFIXED8(p0->x) + (int16_t)(cy * dx.q1);
            w.q1 = __TOFIXED8(w0) + (int16_t)(cy * dw.q1);
            ruv.q0 = __TOFIXED8(u0) + (int16_t)(cy * rduv.q0);
            ruv.q1 = __TOFIXED8(v0) + (int16_t)(cy * rduv.q1);
        }
        drawTextureFragment_fixed(bitmap, x, luv, ruv, w, texture->data, texture->size, texture->size - 1);        
    }
}
