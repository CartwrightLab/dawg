// Test for uniformity
#ifndef GEN
#	define GEN <dawg/details/xorshift_64.h>
#endif

#define RANDOM_GEN_HEADER GEN

#include <iostream>
#include <stdint.h>

#include <unif01.h>
#include <bbattery.h>
#include <dawg/details/mutt.h>

unsigned long get_bits(void *params, void *state) {
	mutt_gen_default *g = static_cast<mutt_gen_default*>(state);
	uint64_t u = g->rand_uint64();
#ifdef LOWER
	u >>= 32;
#endif
	return (u & 0xFFFFFFFFUL);
}

double get_u01(void *params, void *state) {
	mutt_gen_default *g = static_cast<mutt_gen_default*>(state);
	return g->rand_real();
}

void write_xorshift(void *state) {
	printf("hidden");
}

unif01_Gen *create_gen(unsigned int u) {
	uinf01_Gen *gen = new uinf01_gen;
	mutt_gen_default *g = new mutt_gen_default;
	g->seed(u);
	gen->state = g;
	gen->name = #GEN;
	gen->param = NULL;
	gen->GetU01 =  &get_u01;
	gen->GetBits = &get_bits;
	gen->Write = &write_gen;
	return gen;
}

void delete_gen(unif01_Gen *gen) {
	if(NULL == gen)
		return;
	delete gen->state;
	delete gen;
}

int main(int argc, char *argv[]) {
	unsigned int u = 1276590013;
	uinf01_Gen *g = create_gen(u);
	if(argc < 1 || argv[1][0] == "s")
		bbattery_SmallCrush(g);
	else if(argv[1][0] == "m")
		bbattery_Crush(g);
	else if(argv[1][0] == "b"
		bbattery_BigCrush(g);
	else {
		std::cerr << "Unknown Test" << std::endl;
		delete_gen(g);
		return 1;
	}
	delete_gen(g);
	return 0;
}

