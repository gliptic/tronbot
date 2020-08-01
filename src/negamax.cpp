#include "negamax.hpp"

namespace negamax {

u32 ADJACENT_CELLS[5 * 5];

void init() {
	for (int y = 0; y < 5; ++y)
	for (int x = 0; x < 5; ++x) {
		auto p = y * 5 + x;
		u32 move_mask = 0;
		for (int y2 = y - 1; y2 <= y + 1; ++y2)
		for (int x2 = x - 1; x2 <= x + 1; ++x2) {
			if (x2 >= 0 && x2 < 5 && y2 >= 0 && y2 < 5) {
				auto p2 = y2 * 5 + x2;
				if (p != p2) {
					move_mask |= 1 << p2;
				}
			}
		}

		ADJACENT_CELLS[p] = move_mask;
	}
}

}