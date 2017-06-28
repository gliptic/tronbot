#ifndef BOT_h
#define BOT_h

#include <string>
#include "board.hpp"
#include "rand.hpp"

u32 ticks();

struct BaseBot {
	BaseBot() {}
	// Action
	virtual BoardMoves move(int time) = 0;
	// Update
	virtual void round(int time) = 0;
	virtual void update(Board board) = 0;
	// Settings
	virtual void timebank(int time) = 0;
	virtual void time_per_move(int time) = 0;
	virtual void your_bot(std::string name) = 0;
	virtual void your_bot_id(Player player_id) = 0;
	virtual void field_width(int width) = 0;
	virtual void field_height(int height) = 0;
	virtual void player_names(std::string player1, std::string player2) = 0;
	virtual void print_diag(double (&diag)[16][16]) {}

	virtual ~BaseBot();
};

struct Bot : BaseBot {
	Bot(LcgPair& pair);
	// Action
	virtual BoardMoves move(int time) = 0;
	// Update
	virtual void round(int time);
	virtual void update(Board board);
	// Settings
	virtual void timebank(int time);
	virtual void time_per_move(int time);
	virtual void your_bot(std::string name);
	virtual void your_bot_id(Player player_id);
	virtual void field_width(int width);
	virtual void field_height(int height);
	virtual void player_names(std::string player1, std::string player2);

protected:
	Player player_id;
	int width, height;
	GameState state;
	LcgPair& rng;
};

#endif
