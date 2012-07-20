// Test for uniformity
#define DAWG_NO_CONFIG_H
#ifndef GEN
#	define GEN <dawg/details/xorshift_64.h>
#endif

#define RANDOM_GEN_HEADER GEN

#define STRINGIZE(x) #x
#define GEN_NAME STRINGIZE(GEN) 

#include <iostream>
#include <stdint.h>

#include <unif01.h>
#include <bbattery.h>
#include <dawg/details/mutt.h>

using namespace dawg;
using namespace dawg::details;

unsigned long get_bits(void *params, void *state) {
	mutt_gen_default *g = static_cast<mutt_gen_default*>(state);
	uint64_t u = g->rand_uint64();
#ifndef LOWER
	u >>= 32;
#endif
	return (u & 0xFFFFFFFFUL);
}

double get_u01(void *params, void *state) {
	mutt_gen_default *g = static_cast<mutt_gen_default*>(state);
	return g->rand_real();
}

void write_gen(void *state) {
	printf("hidden");
}

unif01_Gen *create_gen(unsigned int u) {
	static char name[] = GEN_NAME;
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

