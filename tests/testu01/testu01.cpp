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
#include <unif01.h>
#include <bbattery.h>
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
	x = revbits(x) << 4;
#endif
#ifdef TEST_LOWER
	return (x & 0xFFFFFFFFUL);	
#else
	return ((x >> 32) & 0xFFFFFFFFUL);
#endif
}

unsigned long get_bits(void *params, void *state) {
	mutt_gen_default *g = static_cast<mutt_gen_default*>(state);
	return to_32(g->rand_native());
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

int main(int argc, const char *argv[]) {
	unsigned int u = 1276590013;
	unif01_Gen *g = create_gen(u);
	if(argc < 1 || argv[1][0] == 's')
		bbattery_SmallCrush(g);
	else if(argv[1][0] == 'm')
		bbattery_Crush(g);
	else if(argv[1][0] == 'b')
		bbattery_BigCrush(g);
	else {
		std::cerr << "Unknown Test" << std::endl;
		delete_gen(g);
		return 1;
	}
	delete_gen(g);
	return 0;
}

