#include "parser.hpp"
using namespace std;

Parser::Parser(Bot &bot) : bot(bot) {
}

string Parser::next_cmd() {
	string k;
	getline(cmdline, k, ' ');
	return k;
}

void Parser::parse() {
	string line;
	while (getline(cin, line)) {
		cmdline.clear();
		cmdline.str(line);
		// Process command
		string cmd = next_cmd();
		if      (cmd == "action")   process_action();
		else if (cmd == "update")   process_update();
		else if (cmd == "settings") process_settings();
	}
}

void Parser::process_action() {
	string cmd = next_cmd();
	if (cmd == "move") {
		auto boardMove = bot.move(stoi(next_cmd()));

		switch (boardMove) {
		case UP:    printf("up\n"); break;
		case DOWN:  printf("down\n"); break;
		case LEFT:  printf("left\n"); break;
		case RIGHT: printf("right\n"); break;
		}
	}
}

void Parser::process_update() {
	next_cmd();
	string cmd = next_cmd();
	if      (cmd == "round") { bot.round(stoi(next_cmd())); }
	else if (cmd == "field") { // Potentially replace with your own boardstate parser
		stringstream ss(next_cmd());
		Board board(ss);
		bot.update(board);
	}
}

void Parser::process_settings() {
	string cmd = next_cmd();
	if      (cmd == "timebank")      bot.timebank(stoi(next_cmd()));
	else if (cmd == "time_per_move") bot.time_per_move(stoi(next_cmd()));
	else if (cmd == "your_bot")      bot.your_bot(next_cmd());
	else if (cmd == "your_botid")   {
		int id = stoi(next_cmd());
		bot.your_bot_id((id == 0) ? Pl1 : Pl2);
	}
	else if (cmd == "field_width")   { width  = stoi(next_cmd()); bot.field_width(width);}
	else if (cmd == "field_height")  { height = stoi(next_cmd()); bot.field_height(height);}
	else if (cmd == "player_names") {
		stringstream args(next_cmd());
		string player1, player2;
		getline(args, player1, ',');
		getline(args, player2, ',');
		bot.player_names(player1,player1);}
}
