#ifndef SILLY_BOT_HPP
#define SILLY_BOT_HPP 1

#include "bot.hpp"

struct SillyBot : Bot {
	int variant;

	SillyBot(LcgPair& rng, int variant);

	// Action
	virtual BoardMoves move(int time);
};

#endif // !SILLY_BOT_HPP
