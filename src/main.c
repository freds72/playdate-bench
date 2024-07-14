//
//  main.c
//  Extension
//

// compile
// cmake .. -G "NMake Makefiles" --toolchain="%PLAYDATE_SDK_PATH%\C_API\buildsupport\arm.cmake" -DCMAKE_BUILD_TYPE=Release
// nmake

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <pd_api.h>

#include "3dmath.h"
#include "drawables.h"
#include "gfx.h"

static int update(void* userdata);
static void init(void* userdata);
const char* fontpath = "/System/Fonts/Asheville-Sans-14-Bold.pft";
LCDFont* font = NULL;
static PlaydateAPI* pd;

void* (*lib3d_realloc)(void* ptr, size_t size);

#define lib3d_malloc(s) lib3d_realloc(NULL, (s))
#define lib3d_free(ptr) lib3d_realloc((ptr), 0)

static void lib3d_setRealloc(void* (*realloc)(void* ptr, size_t size)) {
	lib3d_realloc = realloc;
}

#ifdef _WINDLL
__declspec(dllexport)
#endif
int eventHandler(PlaydateAPI* playdate, PDSystemEvent event, uint32_t arg)
{
	(void)arg; // arg is currently only used for event = kEventKeyPressed

	if ( event == kEventInit )
	{
		pd = playdate;
		const char* err;
		font = pd->graphics->loadFont(fontpath, &err);
		
		if ( font == NULL )
			pd->system->error("%s:%i Couldn't load font %s: %s", __FILE__, __LINE__, fontpath, err);

		// Note: If you set an update callback in the kEventInit handler, the system assumes the game is pure C and doesn't run any Lua code in the game
		pd->system->setUpdateCallback(update, pd);

		// default font
		pd->graphics->setFont(font);

		// init sub modules
		init(pd);		
	}
	
	return 0;
}

// 3d stuff
#define OUTCODE_IN 0
#define OUTCODE_FAR 1
#define OUTCODE_NEAR 2
#define OUTCODE_RIGHT 4
#define OUTCODE_LEFT 8

#define Z_NEAR 0.5f

typedef struct {
	float u;
	float v;
} UV;

typedef struct {
	int vi[4];
	UV uv[4];
	Point3d n;
	Texture* texture;
	float cp;	
} ThreedFace;

typedef struct {
	Point3d* vertices;
	ThreedFace* faces;
	int nfaces;
} ThreedModel;

// 3d models
static Point3d cube_verts[8] = {
	{ .v = {0.f,0.f,0.f} },
	{ .v = {1.f,0.f,0.f} },
	{ .v = {1.f,0.f,1.f} },
	{ .v = {0.f,0.f,1.f} },
	{ .v = {0.f,1.f,0.f} },
	{ .v = {1.f,1.f,0.f} },
	{ .v = {1.f,1.f,1.f} },
	{ .v = {0.f,1.f,1.f} },
};

static Texture side_texture = { 0 };
static Texture top_texture = { 0 };

static ThreedFace cube_faces[6]={
	{.vi = {0,3,2,1}, .texture = &top_texture, .uv = { {.u = 0.f, .v = 0.f }, {.u = 0.f, .v = 63.9f},{.u = 63.9f, .v = 63.9f }, {.u = 63.9f, .v = 0.f } }},
	{.vi = {0,1,5,4}, .texture = &side_texture, .uv = { {.u = 0.f, .v = 0.f }, {.u = 63.9f,.v = 0.f}, {.u = 63.9f, .v = 63.9f }, {.u = 0.f,  .v = 63.9f} }},
	{.vi = {1,2,6,5}, .texture = &side_texture, .uv = { {.u = 0.f, .v = 0.f }, {.u = 63.9f,.v = 0.f}, {.u = 63.9f, .v = 63.9f }, {.u = 0.f,  .v = 63.9f} }},
	{.vi = {2,3,7,6}, .texture = &side_texture, .uv = { {.u = 0.f, .v = 0.f }, {.u = 63.9f,.v = 0.f}, {.u = 63.9f, .v = 63.9f }, {.u = 0.f,  .v = 63.9f} }},
	{.vi = {3,0,4,7}, .texture = &side_texture, .uv = { {.u = 63.9f,.v =0.f },{.u = 0.f, .v = 0.f}, {.u = 0.f,  .v = 63.9f }, {.u = 63.9f, .v = 63.9f} }},
	{.vi = {4,5,6,7}, .texture = &top_texture, .uv=  { {.u = 0.f, .v = 0.f},  {.u = 63.9f,.v = 0.f}, {.u = 63.9f, .v = 63.9f }, {.u = 0.f, .v = 63.9f} }}
};

static ThreedModel cube = {
	.vertices = cube_verts,
	.faces = cube_faces,
	.nfaces = 6
};

// clip polygon against near-z
static int z_poly_clip(const float z, const float flip, Point3duv* in, int n, Point3duv* out) {
    Point3duv v0 = in[n - 1];
    float d0 = flip * (v0.z - z);
    int nout = 0;
    for (int i = 0; i < n; i++) {
        Point3duv v1 = in[i];
        int side = d0 > 0;
        if (side) out[nout++] = (Point3duv){ .pos = { v0.x, v0.y, v0.z }, .u = v0.u, .v = v0.v };
        const float d1 = flip * (v1.z - z);
        if ((d1 > 0) != side) {
            // clip!
            const float t = d0 / (d0 - d1);
            out[nout++] = (Point3duv){ .pos = {
                lerpf(v0.x,v1.x,t),
                lerpf(v0.y,v1.y,t),
                z},
                .u = lerpf(v0.u,v1.u,t),
                .v = lerpf(v0.v,v1.v,t)
            };
        }
        v0 = v1;
        d0 = d1;
    }
    return nout;
}

static int _mode = 0;

static void draw_face(Drawable* drawable, uint8_t* bitmap) {
    DrawableFace* face = &drawable->face;

    const int n = face->n;
    Point3duv* pts = face->pts;
    for (int i = 0; i < n; ++i) {
        // project 
        const float w = 1.0f / pts[i].z;
        pts[i].x = 199.5f +  199.5f * w * pts[i].x;
        pts[i].y = 119.5f -  199.5f * w * pts[i].y;
		pts[i].z = w;
    }

  	//polyfill(pts, n, _ordered_dithers + (int)(face->material * (1.f - shading)) * 32, (uint32_t*)bitmap);
	if (_mode == 0) {
		texfill_baseline(pts, n, face->texture, bitmap);
	}
	else if (_mode == 1) {
		texfill(pts, n, face->texture, bitmap);
	}
	else if (_mode == 2) {
		texfill_fixed(pts, n, face->texture, bitmap);
	}

    // Point3duv* p0 = &pts[n - 1];
	// for (int i = 0; i < n; ++i) {
	// 		Point3duv* p1 = &pts[i];
	// 		if (p0->u) {
	// 				pd->graphics->drawLine((int)p0->x, (int)p0->y, (int)p1->x, (int)p1->y, 1, kColorBlack);
	// 		}
	// 		p0 = p1;
	// }
}

static void push_threeD_model(const ThreedModel* model, const Point3d cv, const Mat4 m) {
    Point3duv tmp[4];		
	for (int i = 0; i < model->nfaces; i++) {
		ThreedFace* f = &model->faces[i];
		// visible?
        if ( v_dot(f->n, cv) - f->cp > 0.01f) {
            // transform
            int outcode = 0xfffffff, is_clipped_near = 0;
            float min_key = FLT_MAX;
            float max_key = -FLT_MAX;
			int n = 4;
			for(int i=0;i<4;i++) {
                Point3duv* res = &tmp[i];

                // project using active matrix
                m_x_v(m, model->vertices[f->vi[i]], &res->p);

                int code =
                    (((Flint) { .f = res->z - Z_NEAR }.i >> 30) & OUTCODE_NEAR) |
                    (((Flint) { .f = res->z - res->x }.i >> 29) & OUTCODE_RIGHT) |
                    (((Flint) { .f = res->z + res->x }.i >> 28) & OUTCODE_LEFT);

                if (res->z < min_key) min_key = res->z;
                if (res->z > max_key) max_key = res->z;
                outcode &= code;
                is_clipped_near |= code;

				// uv's
				res->u = f->uv[i].u;
				res->v = f->uv[i].v;
            }

            // visible?
            if (outcode == 0) {
                Drawable* drawable = pop_drawable(min_key);
                drawable->draw = draw_face;
                drawable->key = min_key;
				DrawableFace* face = &drawable->face;
				face->texture = f->texture;
                if (is_clipped_near & OUTCODE_NEAR) {
                    face->n = z_poly_clip(Z_NEAR, 1.0f, tmp, n, face->pts);
                }
                else {
                    face->n = n;
                    memcpy(face->pts, tmp, n * sizeof(Point3duv));
                }
            }       
        }
    }
}

static void load_texture(const char* name, Texture* texture) {
	const char* err;
	char* path = NULL;

	// read dither table
	pd->system->formatString(&path, "images/%s", name);
	LCDBitmap* bitmap = pd->graphics->loadBitmap(path, &err);

	if (!bitmap)
		pd->system->logToConsole("Failed to load: %s, %s", path, err);

	int w = 0, h = 0, r = 0;
	uint8_t* mask = NULL;
	uint8_t* data = NULL;
	pd->graphics->getBitmapData(bitmap, &w, &h, &r, &mask, &data);
	if (w != h)
		pd->system->logToConsole("Invalid image format: %dx%d", w, h);

	// allocate byte buffer
	uint8_t* tex = lib3d_malloc(w * h);

	for (uint8_t j = 0; j < h; j++, data += w/8) {
		for (uint8_t i = 0; i < w; i++) {
			tex[i + j * w] = (data[i/8] & (0x80>>(i&7)))?0x80:0;
		}
	}
	texture->size = w;
	texture->data = tex;

	pd->graphics->freeBitmap(bitmap);
}

static void init(void* userdata) {
	PlaydateAPI* pd = userdata;
	
	lib3d_setRealloc(pd->system->realloc);

	drawables_init(pd);
	gfx_init(pd);

	// init models
	ThreedModel* model = &cube;

	for(int i=0;i<model->nfaces;i++) {
		ThreedFace* f = &model->faces[i];
		Point3d n0;
		Point3d n1;
		Point3d v0 = model->vertices[f->vi[0]];
		make_v(v0,model->vertices[f->vi[3]],&n0);
		make_v(v0,model->vertices[f->vi[1]],&n1);
		
		v_cross(n0,n1,&f->n);
		v_normz(&f->n);
		f->cp = v_dot(f->n, v0);
	}

	// load & prep textures
	load_texture("crate_side", &side_texture);
	load_texture("crate_top", &top_texture);
}

float total = 0;
int runs = 0;
float cam_dist = 2.f;
float cam_angle = 2.5f;
static int update(void* userdata)
{
	static char* modes[] = {
	"8-pixel strip + fixed16 + q16q16",
	"8-pixel strip + floats",
	"8-pixel strip + fixed16 + unrolled"};
	

	pd->graphics->clear(kColorWhite);
	uint8_t* screen = pd->graphics->getFrame(); // working buffer

	PDButtons pressed;
	PDButtons pushed;
	PDButtons released;
	pd->system->getButtonState(&pressed, &pushed, &released);
	if (released & kButtonA) {
		total = 0;
		runs = 0;
		_mode = (_mode + 1) % 3;
	}
	if (pressed & kButtonUp) {
		cam_dist += 0.1f;
	}
	if (pressed & kButtonDown) {
		cam_dist -= 0.1f;
	}
	if (cam_dist < 0.2f) cam_dist = 0.2f;

	if (pressed & kButtonLeft) {
		cam_angle -= 0.1f;
	}
	if (pressed & kButtonRight) {
		cam_angle += 0.1f;
	}

	const float c = cosf(cam_angle), s = sinf(cam_angle);
	Point3d cam_pos = { .v = {cam_dist * c, cam_dist/2.f, cam_dist * s }};
	Point3d look_at = { .v = {0.f, 0.f, 0.f }};

	// camera matrix
	Mat4 m = {0};
	m_look_at(cam_pos, look_at, m);
	m_inv(m);
	m_inv_translate(m, cam_pos);

	// model matrix
	Mat4 mvv = {0};
	float crank_angle = pd->system->getCrankAngle() / 16.0f;
	const float cc = cosf(crank_angle), ss = sinf(crank_angle);
	Mat4 cube_m = { 0 };
	m_x_m((Mat4) {
		1.f, 0.f, 0.f, 0.f,
			0.f, cc, -ss, 0.f,
			0.f, ss, cc, 0.f,
			0.f, 0.f,0.f, 1.f
	}, (Mat4) {
		1.f, 0.f, 0.f, 0.f,
			0.f, 1.0f, 0.f, 0.f,
			0.f, 0.f, 1.0f, 0.f,
			-0.5f, -0.5f, -0.5f, 1.f
		}, cube_m);
	
	m_x_m(m, cube_m, mvv);

	// cam pos in 3d model space
	Point3d inv_cam_pos;
	m_inv_x_v(cube_m, cam_pos, &inv_cam_pos);

	//
	reset_drawables();

	push_threeD_model(&cube, inv_cam_pos, mvv);

	const int N = 1;

	const float t0 = pd->system->getElapsedTime();
	draw_drawables(screen);
	const float t1 = pd->system->getElapsedTime();

	total += t1 - t0;
	runs += 1;
	const float avg = total / runs;
	char* buf;
	int len = pd->system->formatString(&buf,"Mode: %s\nElapsed time: %f s\n(%i calls/s)\nraw: %f s",modes[_mode], avg, (int)(N / avg), t1 - t0);
	pd->graphics->drawText(buf, len, kASCIIEncoding, 2, 2);
	pd->system->realloc(buf,0);
		
	pd->system->drawFPS(0, 230);

	pd->graphics->markUpdatedRows(0, LCD_ROWS - 1);
	return 1;
}

