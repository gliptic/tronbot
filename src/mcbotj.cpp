#include "mcbotj.hpp"

static bool ucb_tab_init = false;
f32 ucb_tab[1024];

void init_ucb_tab() {
	if (!ucb_tab_init) {
		ucb_tab_init = true;

		ucb_tab[0] = INFINITY;
		for (u32 i = 1; i < 1024; ++i) {
			ucb_tab[i] = CounterF::ucb(c2, f32(i));
		}
	}
}