#ifndef PARSER_h
#define PARSER_h

#include <cstdint>
#include <iostream>
#include <sstream>
#include <queue>
#include "bot.hpp"
#include "game_enums.hpp"

class Parser {
public:
	Parser(Bot &bot);
	void parse();
private:
	Bot& bot;
	int width = 0;
	int height = 0;
	void process_command();
	void process_action();
	void process_update();
	void process_settings();
	std::string next_cmd();
	std::stringstream cmdline;
};

#endif
