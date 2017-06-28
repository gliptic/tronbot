#include "gfx.hpp"

Tigr* scr;
bool debug_board = false;

void render_board(usize board_count, BitBoard** boards, TPixel const* colors, PlayerHeading* headings, f64(*diag)[16][16]) {
	tigrClear(scr, tigrRGB(0, 0, 0));
	for (u32 y = 0; y < 16; ++y)
		for (u32 x = 0; x < 16; ++x) {
			u32 size = 48;
			u32 pshrink = 10;
			u32 rx = x * size, ry = y * size;

			auto p = Vec(x, y);
			if (headings && headings[0].pos == p) {
				tigrFill(scr, rx + pshrink, ry + pshrink, size - pshrink * 2, size - pshrink * 2, tigrRGB(200, 0, 0));
			} else if (headings && headings[1].pos == p) {
				tigrFill(scr, rx + pshrink, ry + pshrink, size - pshrink * 2, size - pshrink * 2, tigrRGB(0, 0, 200));
			} else {
				
				usize filled[5] = {};
				u32 count = 0;
				bool drew = false;
				for (usize i = 0; i < board_count; ++i) {
					
					if (boards[i]->is_wall(p)) {
						filled[count++] = i;
						drew = true;
					}
				}

				for (usize i = 0; i < count; ++i) {
					usize c = filled[i];
					u32 subx = i * size / count;
					u32 nextx = (i + 1) * size / count;
					tigrFill(scr, rx + subx, ry, nextx - subx, size, colors[c]);
				}

				if (count == 0) {
					u8 g = diag ? u8(255.0 * (*diag)[y][x]) : 0;
					tigrFill(scr, rx, ry, size, size, tigrRGB(0, g, 0));
				}
			}
		}
	tigrUpdate(scr);
}