#ifndef GFX_HPP
#include "src/tigr.h"
#include "src/board.hpp"

extern Tigr* scr;
extern bool debug_board;

void render_board(usize board_count, BitBoard** boards, TPixel const* colors, PlayerHeading* headings, f64(*diag)[16][16]);

#endif
