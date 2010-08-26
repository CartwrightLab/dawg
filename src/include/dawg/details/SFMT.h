/** 
 * @file SFMT.h 
 *
 * @brief SIMD oriented Fast Mersenne Twister(SFMT) pseudorandom
 * number generator
 *
 * @author Mutsuo Saito (Hiroshima University)
 * @author Makoto Matsumoto (Hiroshima University)
 *
 * Copyright (C) 2006, 2007 Mutsuo Saito, Makoto Matsumoto and Hiroshima
 * University. All rights reserved.
 *
 * The new BSD License is applied to this software.
 * see LICENSE.txt
 *
 * @note We assume that your system has inttypes.h.  If your system
 * doesn't have inttypes.h, you have to typedef uint32_t and uint64_t,
 * and you have to define PRIu64 and PRIx64 in this file as follows:
 * @verbatim
 typedef unsigned int uint32_t
 typedef unsigned long long uint64_t  
 #define PRIu64 "llu"
 #define PRIx64 "llx"
@endverbatim
 * uint32_t must be exactly 32-bit unsigned integer type (no more, no
 * less), and uint64_t must be exactly 64-bit unsigned integer type.
 * PRIu64 and PRIx64 are used for printf function to print 64-bit
 * unsigned int and 64-bit unsigned int in hexadecimal format.
 */

#ifndef SFMT_H
#define SFMT_H

#include <stdio.h>

/*-----------------
  BASIC DEFINITIONS
  -----------------*/
/** Mersenne Exponent. The period of the sequence 
 *  is a multiple of 2^MEXP-1. */
#define SFMT_MEXP 19937
/** SFMT generator has an internal state array of 128-bit integers,
 * and N is its size. */
#define SFMT_N (SFMT_MEXP / 128 + 1)
/** N32 is the size of internal state array when regarded as an array
 * of 32-bit integers.*/
#define SFMT_N32 (SFMT_N * 4)
/** N64 is the size of internal state array when regarded as an array
 * of 64-bit integers.*/
#define SFMT_N64 (SFMT_N * 2)

#if defined(__STDC_VERSION__) && (__STDC_VERSION__ >= 199901L)
  #include <inttypes.h>
#elif defined(_MSC_VER) || defined(__BORLANDC__)
  typedef unsigned int uint32_t;
  typedef unsigned __int64 uint64_t;
  #define inline __inline
#else
  #include <inttypes.h>
  #if defined(__GNUC__)
    #define inline __inline__
  #endif
#endif

#ifndef PRIu64
  #if defined(_MSC_VER) || defined(__BORLANDC__)
    #define PRIu64 "I64u"
    #define PRIx64 "I64x"
  #else
    #define PRIu64 "llu"
    #define PRIx64 "llx"
  #endif
#endif

#if defined(__GNUC__)
#define ALWAYSINLINE __attribute__((always_inline))
#else
#define ALWAYSINLINE
#endif

#if defined(_MSC_VER)
  #if _MSC_VER >= 1200
    #define PRE_ALWAYS __forceinline
  #else
    #define PRE_ALWAYS inline
  #endif
#else
  #define PRE_ALWAYS inline
#endif

#if defined(__BIG_ENDIAN__) && !defined(__amd64) && !defined(BIG_ENDIAN64)
#define BIG_ENDIAN64 1
#endif
#if defined(HAVE_ALTIVEC) && !defined(BIG_ENDIAN64)
#define BIG_ENDIAN64 1
#endif
#if defined(ONLY64) && !defined(BIG_ENDIAN64)
  #if defined(__GNUC__)
    #error "-DONLY64 must be specified with -DBIG_ENDIAN64"
  #endif
#undef ONLY64
#endif
/*------------------------------------------------------
  128-bit SIMD data type for Altivec, SSE2 or standard C
  ------------------------------------------------------*/
#if defined(HAVE_ALTIVEC)
  #if !defined(__APPLE__)
    #include <altivec.h>
  #endif
/** 128-bit data structure */
union W128_T {
    vector unsigned int s;
    uint32_t u[4];
};
/** 128-bit data type */
typedef union W128_T w128_t;

#elif defined(HAVE_SSE2)
  #include <emmintrin.h>

/** 128-bit data structure */
union W128_T {
    __m128i si;
    uint32_t u[4];
};
/** 128-bit data type */
typedef union W128_T w128_t;

#else

/** 128-bit data structure */
struct W128_T {
    uint32_t u[4];
};
/** 128-bit data type */
typedef struct W128_T w128_t;

#endif

/** the 128-bit internal state array */
struct SFMT_T {
	/** the 128-bit internal state array */
	w128_t status[SFMT_N];
	/** index counter to the 32-bit internal state array */
	int idx;
};
typedef struct SFMT_T sfmt_t;

/** the 32bit integer pointer to the 128-bit internal state array */
inline uint32_t * psfmt32(sfmt_t *sfmt) { return &sfmt->status[0].u[0]; }
#if !defined(BIG_ENDIAN64) || defined(ONLY64)
/** the 64bit integer pointer to the 128-bit internal state array */
inline uint64_t * psfmt64(sfmt_t *sfmt) { return (uint64_t *)&sfmt->status[0].u[0]; }
#endif


uint32_t sfmt_gen_rand32(sfmt_t *sfmt);
uint64_t sfmt_gen_rand64(sfmt_t *sfmt);
void sfmt_fill_array32(sfmt_t *sfmt, uint32_t *array, int size);
void sfmt_fill_array64(sfmt_t *sfmt, uint64_t *array, int size);
void sfmt_init_gen_rand(sfmt_t *sfmt, uint32_t seed);
void sfmt_init_by_array(sfmt_t *sfmt, uint32_t *init_key, int key_length);
const char *sfmt_get_idstring(void);
int sfmt_get_min_array_size32(void);
int sfmt_get_min_array_size64(void);

/* These real versions are due to Isaku Wada */
/** generates a random number on [0,1]-real-interval */
inline static double to_real1(uint32_t v)
{
    return v * (1.0/4294967295.0); 
    /* divided by 2^32-1 */ 
}

/** generates a random number on [0,1]-real-interval */
inline static double sfmt_genrand_real1(sfmt_t *sfmt)
{
    return to_real1(sfmt_gen_rand32(sfmt));
}

/** generates a random number on [0,1)-real-interval */
inline static double to_real2(uint32_t v)
{
    return v * (1.0/4294967296.0); 
    /* divided by 2^32 */
}

/** generates a random number on [0,1)-real-interval */
inline static double sfmt_genrand_real2(sfmt_t *sfmt)
{
    return to_real2(sfmt_gen_rand32(sfmt));
}

/** generates a random number on (0,1)-real-interval */
inline static double to_real3(uint32_t v)
{
    return (((double)v) + 0.5)*(1.0/4294967296.0); 
    /* divided by 2^32 */
}

/** generates a random number on (0,1)-real-interval */
inline static double sfmt_genrand_real3(sfmt_t *sfmt)
{
    return to_real3(sfmt_gen_rand32(sfmt));
}
/** These real versions are due to Isaku Wada */

/** generates a random number on [0,1) with 53-bit resolution*/
inline static double to_res53(uint64_t v) 
{ 
    return v * (1.0/18446744073709551616.0L);
}

/** generates a random number on [0,1) with 53-bit resolution from two
 * 32 bit integers */
inline static double to_res53_mix(uint32_t x, uint32_t y) 
{ 
    return to_res53(x | ((uint64_t)y << 32));
}

/** generates a random number on [0,1) with 53-bit resolution
 */
inline static double sfmt_genrand_res53(sfmt_t *sfmt) 
{ 
    return to_res53(sfmt_gen_rand64(sfmt));
} 

/** generates a random number on [0,1) with 53-bit resolution
    using 32bit integer.
 */
inline static double sfmt_genrand_res53_mix(sfmt_t *sfmt) 
{ 
    uint32_t x, y;

    x = sfmt_gen_rand32(sfmt);
    y = sfmt_gen_rand32(sfmt);
    return to_res53_mix(x, y);
} 
#endif
