#ifndef RAND_HPP
#define RAND_HPP 1

#include "game_enums.hpp"

struct TinyEntropy {
	u32 x;

	TinyEntropy(u32 x) : x(x) {
	}

	u32 get_u32(u32 max) {
		u64 v = u64(this->x) * max;
		this->x = u32(v);
		return u32(v >> 32);
	}
};

struct LcgPair {
	u32 s[2];

	LcgPair()
		: LcgPair(0, 0) {
	}

	LcgPair(u32 a_init, u32 b_init) {
		s[0] = a_init;
		s[1] = b_init;
	}

	static LcgPair from_bits(u32 a, u32 b) {
		a = a * 2654435761 ^ b;
		b = b * 2654435761 ^ a;
		return LcgPair(a, b);
	}

	u32 next() {
		u32 x = s[0];
		s[0] = s[0] * 29943829 + 0xffff;
		return x;
	}

	TinyEntropy tiny_entropy() {
		return TinyEntropy(next());
	}

	static u32 lim(u32 v, u32 max) {
		return u32(u64(v) * max >> 32);
	}

	u32 get_u32(u32 max) {
		return lim(this->next(), max);
	}
};

#endif
