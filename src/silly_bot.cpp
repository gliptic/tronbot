#include "silly_bot.hpp"

SillyBot::SillyBot(LcgPair& rng, int variant)
	: Bot(rng), variant(variant) {
}

BoardMoves SillyBot::move(int time) {
	//auto moves = state.board.legal_moves(player_id);
	//auto moves = state.board.best_fill_moves(player_id);
	auto moves = best_dist_moves(state.board, player_id, this->rng);
	//auto moves = variant ? state.board.legal_moves2(player_id) : state.board.legal_moves(player_id);
	if (moves == NONE) {
		return UP;
	}

	while (true) {
		u32 move = this->rng.next() >> (32 - 2);
		if (moves & (1 << move)) {
			return BoardMoves(1 << move);
		}
	}
}
