#define BOOST_TEST_MODULE Dawg::testu01
#define BOOST_TEST_DYN_LINK

#include "../boost_test_helper.h"

// Test for uniformity
#define DAWG_NO_CONFIG_H
#ifndef TEST_GEN
#	define TEST_GEN <dawg/details/xorshift_64.h>
#endif

#define RANDOM_GEN_HEADER TEST_GEN

#define XSTR(s) XSTR_(s)
#define XSTR_(s) #s
#define GEN_NAME XSTR(RANDOM_GEN_HEADER)

#include <iostream>
#include <stdint.h>
#include <cstdio>

extern "C" {
#if defined(__FreeBSD__)
#include <TestU01/unif01.h>
#include <TestU01/bbattery.h>
#include <TestU01/gdefconf.h>
#else
#include <unif01.h>
#include <bbattery.h>
#endif // defined
}

#include <dawg/details/mutt.h>

using namespace dawg;
using namespace dawg::details;

uint64_t revbits(uint64_t x) {
	uint64_t y = 0;
	for(int i=0;i<8;++i) {
		unsigned char b = (unsigned char)x;
		b = ((b * 0x80200802ULL) & 0x0884422110ULL) * 0x0101010101ULL >> 32;
		y = (y << 8) + b;
		x =  x >> 8;
	}
	return y;
}

uint32_t to_32(uint32_t x) {
	return x;
}
uint32_t to_32(uint64_t x) {
#ifdef TEST_REV
	x = revbits(x) ;//<< 4;
#endif
#ifdef TEST_LOWER
	return (x & 0xFFFFFFFFUL);
#else
	return ((x >> 32) & 0xFFFFFFFFUL);
#endif
}

unsigned long get_bits(void *params, void *state) {
	mutt_gen_default *g = static_cast<mutt_gen_default*>(state);
#ifndef TEST_MIX
	return to_32(g->rand_native());
#else
	unsigned long u=0;
	for(int i=0;i<8;++i)
		u = u << 8 | (g->rand_native() & 255);
	return u;
#endif
}

double get_u01(void *params, void *state) {
	mutt_gen_default *g = static_cast<mutt_gen_default*>(state);
	return g->rand_real();
}

void write_gen(void *state) {
	printf("N/A");
}

unif01_Gen *create_gen(unsigned int u) {
	static char name[] = GEN_NAME
#ifdef DAWG_DISABLE_WEYL_GENERATOR
		"-now"
#endif
#ifdef TEST_LOWER
		"-low"
#endif
#ifdef TEST_REV
		"-rev"
#endif
#ifdef TEST_MIX
		"-mix"
#endif
	;
	unif01_Gen *gen = new unif01_Gen;
	mutt_gen_default *g = new mutt_gen_default;
	g->seed(u);
	gen->state = g;
	gen->name = &name[0];
	gen->param = NULL;
	gen->GetU01 =  &get_u01;
	gen->GetBits = &get_bits;
	gen->Write = &write_gen;
	return gen;
}

void delete_gen(unif01_Gen *gen) {
	if(NULL == gen)
		return;
	delete static_cast<mutt_gen_default*>(gen->state);
	delete gen;
}

BOOST_AUTO_TEST_CASE(test_small_crush)
{
    unsigned int u = 1276590013;
    unif01_Gen *g = create_gen(u);
    bbattery_SmallCrush(g);
    delete_gen(g);
}

//BOOST_AUTO_TEST_CASE(test_crush)
//{
//    unsigned int u = 1276590013;
//    unif01_Gen *g = create_gen(u);
//    bbattery_Crush(g);
//    delete_gen(g);
//}

//BOOST_AUTO_TEST_CASE(test_big_crush)
//{
//    unsigned int u = 1276590013;
//    unif01_Gen *g = create_gen(u);
//    bbattery_BigCrush(g);
//    delete_gen(g);
//}

