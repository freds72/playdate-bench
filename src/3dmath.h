#ifndef _lib3d_3dmath_h
#define _lib3d_3dmath_h

#include <stdint.h>
#include <stddef.h>

#define PI 3.1415927410125732421875f
#define MAT4x4 16
#define VEC3 3

#ifndef max
    #define max(a,b) (((a) > (b)) ? (a) : (b))
    #define min(a,b) (((a) < (b)) ? (a) : (b))
#endif // !max

// misc. math typedefs
typedef struct Point3d {
  union {
    // named access
    struct {
      float x;
      float y;
      float z;
    };
    // values array
    float v[VEC3];
  };
} Point3d;

typedef struct {
    union {
        struct {
            float x;
            float y;
        };
        float v[2];
    };
} Point2d;


// 3d point with uv coordinate
typedef struct {
    union {
        struct Point3d p;
        struct {
            float x;
            float y;
            float z;
        };
        // values array
        float pos[VEC3];
    };
    union {
        struct {
            float u;
            float v;
        };
        float uv[2];
    };
} Point3duv;

// range helpers
typedef struct {
    int min;
    int max;
} IntRange;

typedef struct {
    float min;
    float max;
} FloatRange;

// aliases a float to a 32bits memory address
typedef struct {
    union {
        float f;
        uint32_t i;
    };
} Flint;

// basic texture definition
typedef struct {
    uint8_t size;
    uint8_t* data;
} Texture;

// fixed type helper
typedef struct {
    union {
        struct {
            int8_t b0;
            int8_t b1;
            int8_t b2;
            int8_t b3;
        };
        struct {
            int16_t q0;
            int16_t q1;
        };
        uint32_t i32;
    };
} Fixed;

// matrix struct
typedef float Mat4[MAT4x4];

#define __STATIC_INLINE static inline

// convert a tau angle [0;1] into a radian angle
#if TARGET_PLAYDATE
__attribute__((always_inline)) __STATIC_INLINE
#else
static __forceinline
#endif
float detauify(const float tau) {
    return tau * 2 * PI;
}

// convert a rad angle [0;2*PI] into a tau angle [0;1]
#if TARGET_PLAYDATE
__attribute__((always_inline)) __STATIC_INLINE
#else
static __forceinline
#endif
float tauify(const float rad) {
    return rad / (2 * PI);
}

// returns a random number between 0-1
float randf();

// returns a rand number between [0;max[
int randi(const int max);

// https://benpfaff.org/writings/clc/shuffle.html
// in-place shuffling of int array
void shuffle(int* array, const size_t n);

// lerp between 2 float values
#if TARGET_PLAYDATE
__attribute__((always_inline)) __STATIC_INLINE
#else
static __forceinline
#endif
float lerpf(const float a, const float b, const float t) {
  return a + (b - a) * t;
}

// lerp between 2 int values
#if TARGET_PLAYDATE
__attribute__((always_inline)) __STATIC_INLINE
#else
static __forceinline
#endif
int lerpi(const int a, const int b, const float t) {
    return (int)lerpf((float)a, (float)b, t);
}

#if TARGET_PLAYDATE
__attribute__((always_inline)) __STATIC_INLINE
#else
static __forceinline
#endif
float v_dot(const Point3d v0, const Point3d v1) {
    // faster??
    float acc = 0;
    acc += v0.v[0] * v1.v[0];
    acc += v0.v[1] * v1.v[1];
    acc += v0.v[2] * v1.v[2];
    return acc;
}

#if TARGET_PLAYDATE
__attribute__((always_inline)) __STATIC_INLINE
#else
static __forceinline
#endif
void m_x_v(const Mat4 m, const Point3d v, Point3d* out) {
    const float x = v.x, y = v.y, z = v.z;
    out->v[0] = m[0] * x + m[4] * y + m[8] * z + m[12];
    out->v[1] = m[1] * x + m[5] * y + m[9] * z + m[13];
    out->v[2] = m[2] * x + m[6] * y + m[10] * z + m[14];
}

void make_v(const Point3d a, Point3d b, Point3d* out);

// returns 1/len
float v_normz(Point3d* a);
void v_cross(const Point3d a, const Point3d b, Point3d* out);

// matrix multiply
void m_x_m(const Mat4  a, const Mat4 b, Mat4 out);
// translate matrix by vector v
void m_x_translate(const Mat4 a, const Point3d v, Mat4 out);
// multiply by y rotation
void m_x_y_rot(const Mat4 a, const float angle, Mat4 out);
// matrix vector multiply invert
// inc.position
void m_inv_x_v(const Mat4 m, const Point3d v, Point3d* out);
// interpolate points
void v_lerp(const Point3d a, const Point3d b, const float t, Point3d* out);

void m_look_at(const Point3d from, const Point3d to, Mat4 out);

void m_inv(Mat4 m);
void m_inv_translate(Mat4 m, const Point3d v);

#endif