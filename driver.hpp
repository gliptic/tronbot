#ifndef DRIVER_HPP
#define DRIVER_HPP 1

#include <utility>
#include <memory>
#include <vector>
#include <cmath>
#include <cassert>

#include "src/bot.hpp"
#include "src/board.hpp"
#include "src/rand.hpp"

using std::unique_ptr;
using std::vector;

struct Game : GameState {
	unique_ptr<BaseBot> bot[2];
	u32 step;

	Game(LcgPair& rng, unique_ptr<BaseBot> bot1, unique_ptr<BaseBot> bot2)
		: GameState(GameState::standard()), step(0) {
		bot[0] = move(bot1);
		bot[1] = move(bot2);

		for (int i = 0; i < 2; ++i) {
			bot[i]->field_width(16);
			bot[i]->field_height(16);
			bot[i]->your_bot_id(Player(i));
		}
	}

	BoardMovesJ update_bots() {
		BoardMovesJ new_moves(NONE, NONE);

		for (int i = 0; i < 2; ++i) {
			bot[i]->update(this->board);
			new_moves.p[i] = bot[i]->move(200);
		}

		return new_moves;
	}

	u32 update_(BoardMovesJ new_moves) {
		
		++step;

		return update(new_moves);
	}
};

#endif // !DRIVER_HPP
