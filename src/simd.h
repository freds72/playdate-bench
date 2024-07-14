
#ifndef _lib3d_simd_h
#define _lib3d_simd_h

#include <stdint.h>

typedef int32_t q31_t;

#if TARGET_PLAYDATE
#define INLINE static __attribute__((always_inline)) inline
#else
#define INLINE static __forceinline
#endif

INLINE uint32_t swap(uint32_t n)
{
#if TARGET_PLAYDATE
    //return __REV(n);
    uint32_t result;

    __asm volatile ("rev %0, %1" : "=l" (result) : "l" (n));
    return(result);
#else
    return ((n & 0xff000000) >> 24) | ((n & 0xff0000) >> 8) | ((n & 0xff00) << 8) | (n << 24);
#endif
}

INLINE uint32_t __SADD16(uint32_t op1, uint32_t op2)
{
#if TARGET_PLAYDATE
    uint32_t result;

    __asm volatile ("sadd16 %0, %1, %2" : "=r" (result) : "r" (op1), "r" (op2));
    return(result);
#else
  return  
    (uint32_t)((int16_t)(op1&0xffff) + ((int16_t)(op2&0xffff))) |
    ((uint32_t)(((int16_t)(op1>>8) + (int16_t)(op2>>8))))<<8;
#endif
}

typedef struct {
    union {
        struct {
            int16_t q0;
            int16_t q1;
        };
        uint32_t i32;
    };
} q16_t;

INLINE uint32_t __SSUB16(uint32_t op1, uint32_t op2)
{
#if TARGET_PLAYDATE
    uint32_t result;

    __asm volatile ("ssub16 %0, %1, %2" : "=r" (result) : "r" (op1), "r" (op2));
    return(result);
#else
    q16_t a = { .i32 = op1 }, b = { .i32 = op2 };
    q16_t res = {
        .q0 = a.q0 - b.q0,
        .q1 = a.q1 - b.q1 
    };
    return res.i32;
#endif
}

INLINE uint32_t __SMLAD(uint32_t x, uint32_t y, uint32_t sum)
{
#if TARGET_PLAYDATE
  uint32_t result;

  __asm volatile ("smlad %0, %1, %2, %3" : "=r" (result) : "r" (x), "r" (y), "r" (sum) );
  return(result);
#else
 return ((uint32_t)(((((q31_t)x << 16) >> 16) * (((q31_t)y << 16) >> 16)) +
                       ((((q31_t)x      ) >> 16) * (((q31_t)y      ) >> 16)) +
                       ( ((q31_t)sum    )                                  )   ));
#endif
}

#define FIXED16_SHIFT 16

// convert the given float into a 16:16 fixed point (int32)
INLINE int32_t __TOFIXED16(float x)
{
    // will corectly generate a vcvt asm instruction
    return (int32_t)(x * (1 << FIXED16_SHIFT));
}

#define FIXED8_SHIFT 6

// convert float to 12:4 fixed (int16)
INLINE int16_t __TOFIXED8(float x)
{
    // will corectly generate a vcvt asm instruction
    return (int16_t)(x * (1 << FIXED8_SHIFT));
}

#endif
