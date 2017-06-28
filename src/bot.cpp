#include "bot.hpp"
using namespace std;

#ifdef _WIN32
#include <windows.h>
u32 ticks() {
	return GetTickCount();
}
#else
#define _POSIX_C_SOURCE 200809L
#include <time.h>

static bool ticks_hasstart = false;
static time_t start_sec;

u32 ticks() {
	struct timespec spec;
	if (!ticks_hasstart) {
		clock_gettime(CLOCK_MONOTONIC, &spec);
		start_sec = spec.tv_sec;
		ticks_hasstart = true;
	}

	clock_gettime(CLOCK_MONOTONIC, &spec);
	return u32(spec.tv_sec - start_sec) * 1000 + spec.tv_nsec / 1000000;
}
#endif

Bot::Bot(LcgPair& rng)
	: rng(rng) {
}

void Bot::round(int time) {  };

void Bot::update(Board board) {
	auto prev_move0 = board.move_from(0, this->state.board);
	auto prev_move1 = board.move_from(1, this->state.board);
	this->state.board = board;
	this->state.board.headings[0].prev_move = prev_move0;
	this->state.board.headings[1].prev_move = prev_move1;
}

void Bot::timebank(int time) { };
void Bot::time_per_move(int time) { };
void Bot::your_bot(string name) { };
void Bot::your_bot_id(Player player_id) { this->player_id = player_id;  };
void Bot::field_width(int width) { this->width = width; };
void Bot::field_height(int height) { this->height = height; };
void Bot::player_names(string player1, string player2) { };

BaseBot::~BaseBot() {
}
