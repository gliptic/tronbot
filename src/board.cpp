#include "board.hpp"
#include <cstring>
#include <algorithm>
#include <intrin.h>
#pragma intrinsic(_BitScanForward)  

#include "../gfx.hpp"

Board2::Board2() {
	this->clear();
}

Board2::Board2(std::stringstream &stream) {
	this->read(stream);
}

struct D2 {
	BitBoard closed;
	BitBoard opened[2];
	u32 score[2];

	D2(BitBoard const& b)
		: closed(b) {
		opened[0].clear_no_walls();
		opened[1].clear_no_walls();
		score[0] = 0;
		score[1] = 0;
	}

	u32 expand(BitBoard& opened_, u32& begin, u32& end) {
		u32 new_opened_from_above = 0;

		u32 b = begin, e = end;
		u32 diff = 0;

		if (opened_.get_row(b) & ~closed.get_row(b - 1)) {
			begin = --b;
			diff = 1;
		}

		for (u32 y = b; y < e; ++y) {
			u32 cur_closed = closed.get_row(y);
			u32 new_opened_from_below = opened_.get_row(y + 1);

			u32 cur_opened = opened_.get_row(y);
			u32 new_opened = ((cur_opened << 1) | (cur_opened >> 1) | new_opened_from_above | new_opened_from_below) /*& ~cur_opened*/ & ~cur_closed;

			diff |= new_opened & ~cur_opened;
			opened_.bits[y + 1] |= new_opened;

			new_opened_from_above = cur_opened;
		}

		if (new_opened_from_above) {
			// May need to extend down
			u32 new_opened = new_opened_from_above & ~closed.get_row(e);
			if (new_opened) {
				assert(opened_.get_row(e) == 0);
				opened_.set_row(e, new_opened);
				end = ++e;
				diff = 1;
			}
		}

		return diff;
	}

	void close(BitBoard& opened_, u32 b, u32 e) {
		for (u32 y = b; y < e; ++y) {
			closed.bits[y + 1] |= opened_.get_row(y);
		}
	}

	u32 count(BitBoard& opened_, u32 b, u32 e) {
		u32 s = 0;
		for (u32 y = b; y < e; ++y) {
			s += popcount(opened_.get_row(y));
		}
		return s;
	}

	void estimate(Vec pos1, Vec pos2) {
		opened[0].set(pos1);
		opened[1].set(pos2);

		u32 begin[2], end[2];
		begin[0] = pos1.v[1];
		end[0] = pos1.v[1] + 1;
		begin[1] = pos2.v[1];
		end[1] = pos2.v[1] + 1;

		u32 diff;
		for (;;) {
			diff = expand(opened[0], begin[0], end[0]);
			diff |= expand(opened[1], begin[1], end[1]);
			if (!diff) break;

			if (false && debug_board) {
				BitBoard* boards[3] = { &closed, opened + 0, opened + 1 };
				TPixel colors[3] = { tigrRGB(100, 100, 100), tigrRGB(200, 100, 100), tigrRGB(100, 100, 200) };
				render_board(3, boards, colors, 0, 0);
			}
			//if (debug_board) render_board(opened[0], 0, 0);
			close(opened[0], begin[0], end[0]);
			//if (debug_board) render_board(opened[1], 0, 0);
			close(opened[1], begin[1], end[1]);
			//if (debug_board) render_board(closed, 0, 0);

			if (false && debug_board) {
				BitBoard* boards[3] = { &closed, opened + 0, opened + 1 };
				TPixel colors[3] = { tigrRGB(100, 100, 100), tigrRGB(200, 100, 100), tigrRGB(100, 100, 200) };
				render_board(3, boards, colors, 0, 0);
			}
		}

		score[0] = count(opened[0], begin[0], end[0]);
		score[1] = count(opened[1], begin[1], end[1]);
	}
};

#define COMBINED 0
#define CLOSE_ALL 1

struct D3 {
	BitBoard closed;
	BitBoard opened[2];
	u32 score[2];

	D3(BitBoard const& b)
		: closed(b) {
#if COMBINED
		closed.remove_walls();
#endif
		opened[0].clear_no_walls();
		opened[1].clear_no_walls();
		score[0] = 0;
		score[1] = 0;
	}

	u32 expand(u32& begin, u32& end) {
#if COMBINED
		u32 new_opened_from_above = 0;
#else
		u32 new_opened_from_above0 = 0, new_opened_from_above1 = 0;
#endif

		u32 b = begin, e = end;
		u32 diff = 0;

		u32 before_b_mask = ~closed.get_row(b - 1);
#if COMBINED
		before_b_mask = (before_b_mask << 16) | before_b_mask;

		if (opened[0].get_row(b) & before_b_mask) {
			begin = --b;
			diff = 1;
		}
#else
		if ((opened[0].get_row(b) | opened[1].get_row(b)) & before_b_mask) {
			begin = --b;
			diff = 1;
		}
#endif

		for (u32 y = b; y < e; ++y) {

#if COMBINED
			u32 cur_closed = closed.get_row(y);
			cur_closed = (cur_closed << 16) | cur_closed;
			u32 new_opened_from_below = opened[0].get_row(y + 1);

			u32 cur_opened = opened[0].get_row(y);

			u32 mask_a = 0xfffeffff;
			u32 mask_b = 0xffff7fff;

			u32 new_opened = (((cur_opened << 1) & mask_a)
				| ((cur_opened >> 1) & mask_b)
				| new_opened_from_above
				| new_opened_from_below) & ~cur_closed;

			diff |= new_opened & ~cur_opened;
			opened[0].bits[y + 1] |= new_opened;

			closed.bits[y + 1] |= new_opened;
			new_opened_from_above = cur_opened;
#else
			u32 cur_closed = closed.get_row(y);
			u32 new_opened_from_below0 = opened[0].get_row(y + 1);
			u32 new_opened_from_below1 = opened[1].get_row(y + 1);

			u32 cur_opened0 = opened[0].get_row(y);
			u32 cur_opened1 = opened[1].get_row(y);
			u32 new_opened0 = ((cur_opened0 << 1) | (cur_opened0 >> 1) | new_opened_from_above0 | new_opened_from_below0) & ~cur_closed;
			u32 new_opened1 = ((cur_opened1 << 1) | (cur_opened1 >> 1) | new_opened_from_above1 | new_opened_from_below1) & ~cur_closed;

			diff |= new_opened0 & ~cur_opened0;
			diff |= new_opened1 & ~cur_opened1;
			opened[0].bits[y + 1] |= new_opened0;
			opened[1].bits[y + 1] |= new_opened1;

			closed.bits[y + 1] |= new_opened0 | new_opened1;
			new_opened_from_above0 = cur_opened0;
			new_opened_from_above1 = cur_opened1;
#endif
		}

#if COMBINED
		if (new_opened_from_above) {
			// May need to extend down
			u32 cur_closed = closed.get_row(e);
			cur_closed = (cur_closed << 16) | cur_closed;
			u32 new_opened = new_opened_from_above & ~cur_closed;
			if (new_opened) {
				opened[0].set_row(e, new_opened);
				end = ++e;
				diff = 1;
			}
		}
#else
		if (new_opened_from_above0 || new_opened_from_above1) {
			// May need to extend down
			u32 cur_closed = closed.get_row(e);
			u32 new_opened0 = new_opened_from_above0 & ~cur_closed;
			u32 new_opened1 = new_opened_from_above1 & ~cur_closed;
			if (new_opened0 || new_opened1) {
				opened[0].set_row(e, new_opened0);
				opened[1].set_row(e, new_opened1);
				closed.bits[e + 1] |= new_opened0 | new_opened1;
				end = ++e;
				diff = 1;
			}
		}
#endif

		return diff;
	}

	u32 count(BitBoard& opened_, u32 b, u32 e) {
		u32 s = 0;
		for (u32 y = b; y < e; ++y) {
			s += popcount(opened_.get_row(y));
		}
		return s;
	}

	u32 count2(BitBoard& opened_, u32 b, u32 e) {
		u32 s0 = 0, s1 = 0;
		u32 pattern0 = 0xAAAA, pattern1 = 0xAAAA >> 1;
		for (u32 y = b; y < e; y += 2) {
			u32 row = opened_.get_row(y);
			s0 += popcount(row & pattern0);
			s1 += popcount(row & pattern1);
			if (y + 1 >= e) break;
			row = opened_.get_row(y + 1);
			s0 += popcount(row & pattern1);
			s1 += popcount(row & pattern0);
		}

		return 2 * std::min(s0, s1);
	}

	Vec count3(BitBoard& opened_, u32 b, u32 e) {
		u32 s0 = 0, s1 = 0;
		for (u32 y = b; y < e; ++y) {
			u32 row = opened_.get_row(y);
			s0 += popcount(row & 0xffff);
			s1 += popcount(row >> 16);
		}

		return Vec(s0, s1);
	}

	void estimate(Vec pos1, Vec pos2, bool patterned) {
		opened[0].set(pos1);
		opened[1].set(pos2);

		u32 begin = std::min(pos1.v[1], pos2.v[1]);
		u32 end = std::max(pos1.v[1] + 1, pos2.v[1] + 1);

		for (;;) {
			u32 diff = expand(begin, end);
			if (!diff) break;

			if (false && debug_board) {
				BitBoard* boards[3] = { &closed, opened + 0, opened + 1 };
				TPixel colors[3] = { tigrRGB(100, 100, 100), tigrRGB(200, 100, 100), tigrRGB(100, 100, 200) };
				render_board(3, boards, colors, 0, 0);
			}
		}

#if COMBINED
		auto c = count3(opened[0], begin, end);
		score[0] = c.v[0];
		score[1] = c.v[1];
#else
		if (patterned) {
			score[0] = count2(opened[0], begin, end);
			score[1] = count2(opened[1], begin, end);
		} else {
			score[0] = count(opened[0], begin, end);
			score[1] = count(opened[1], begin, end);
		}
#endif
	}
};

template<bool CloseAll>
struct D4 {
	BitBoard closed;
	BitBoard opened;
	u8 dist[2][16][16];
	u32 score[2];
	i32 min[2];

	D4(BitBoard const& b)
		: closed(b) {
		closed.dup();
		opened.clear_no_walls();
		score[0] = score[1] = 0;
		min[0] = min[1] = -1;

		memset(dist, 0xff, sizeof(dist));
	}

	inline void reg_new_bits(u32 y, u32 d, u32 new_opened) {
		if (!CloseAll) {
			unsigned long mask = new_opened;
			unsigned long offs;
			while (_BitScanForward(&offs, mask)) {
				mask ^= 1 << offs;
				dist[offs >> 4][y][offs & 15] = d;
			}
		}
	}

	u32 expand(u32& begin, u32& end, u8 d) {
		u32 new_opened_from_above = 0;

		u32 b = begin, e = end;
		u32 diff = 0;

		{
			u32 before_b_mask = ~closed.get_row(b - 1);
			u32 new_opened_from_below = opened.get_row(b) & before_b_mask;
			if (new_opened_from_below) {
				u32 new_opened = new_opened_from_below;
				reg_new_bits(b - 1, d, new_opened);
				opened.or_row(b - 1, new_opened);
				if (CloseAll)
					closed.or_row(b - 1, new_opened | ((new_opened >> 16) | (new_opened << 16)));
				diff = 1;
				--begin;
			}
		}

		for (u32 y = b; y < e; ++y) {

			u32 cur_closed = closed.get_row(y);
			u32 new_opened_from_below = opened.get_row(y + 1);
			u32 cur_opened = opened.get_row(y);

			u32 mask_a = 0xfffeffff;
			u32 mask_b = 0xffff7fff;

			u32 spread = (((cur_opened << 1) & mask_a)
				| ((cur_opened >> 1) & mask_b)
				| new_opened_from_above
				| new_opened_from_below);

			u32 new_opened = CloseAll ? (spread & ~cur_closed) : (spread & ~(cur_closed | cur_opened));

			reg_new_bits(y, d, new_opened);

			diff |= new_opened;
			opened.or_row(y, new_opened);
			if (CloseAll) {
				closed.or_row(y, new_opened | ((new_opened >> 16) | (new_opened << 16)));
			}

			new_opened_from_above = cur_opened;
		}

		if (new_opened_from_above) {
			// May need to extend down
			u32 cur_closed = closed.get_row(e);
			u32 new_opened = new_opened_from_above & ~cur_closed;
			if (new_opened) {
				reg_new_bits(e, d, new_opened);
				opened.set_row(e, new_opened);
				if (CloseAll) {
					closed.or_row(e, new_opened | ((new_opened >> 16) | (new_opened << 16)));
				}

				end = ++e;
				diff = 1;
			}
		}

		return diff;
	}

	u32 count(BitBoard& opened_, u32 b, u32 e) {
		u32 s = 0;
		for (u32 y = b; y < e; ++y) {
			s += popcount(opened_.get_row(y));
		}
		return s;
	}

	u32 count_pat(BitBoard& opened_, u32 pl, Vec start, u32 b, u32 e) {
		u32 s0 = 0, s1 = 0;
		u32 shift = 16 * pl;
		for (u32 y = b; y < e; y += 2) {
			u32 row = opened_.get_row(y) >> shift;
			u32 pattern0 = 0xAAAA >> (y & 1);
			u32 pattern1 = 0xAAAA >> ((y + 1) & 1);
			s0 += popcount(row & pattern0);
			s1 += popcount(row & pattern1);
			if (y + 1 >= e) break;
			row = opened_.get_row(y + 1) >> shift;
			s0 += popcount(row & pattern1);
			s1 += popcount(row & pattern0);
		}

		u32 start_color = ((start.x() ^ start.y()) + 1) & 1;
		if (start_color == 0) {
			return std::min(s0, s1 + 1) + std::min(s1, s0);
		} else {
			return std::min(s1, s0 + 1) + std::min(s0, s1);
		}
	}

	Vec count3(BitBoard& opened_, u32 b, u32 e) {
		if (CloseAll) {
			u32 s0 = 0, s1 = 0;
			for (u32 y = b; y < e; ++y) {
				u32 row = opened_.get_row(y);
				s0 += popcount(row & 0xffff);
				s1 += popcount(row >> 16);
			}

			return Vec(s0, s1);
		} else {
			u32 s0 = 0, s1 = 0;
			for (u32 y = b; y < e; ++y) {

				for (u32 x = 0; x < 16; ++x) {
					u32 d0 = dist[0][y][x];
					u32 d1 = dist[1][y][x];
					s0 += (d0 != 0xff) && d0 <= d1;
					s1 += (d1 != 0xff) && d1 <= d0;
				}
			}
			return Vec(s0, s1);
		}
	}

	void estimate(Vec pos1, Vec pos2, bool patterned, LcgPair& rng) {
		opened.set(pos1);
		opened.set(Vec(pos2.x() + 16, pos2.y()));

		u32 begin = std::min(pos1.v[1], pos2.v[1]);
		u32 end = std::max(pos1.v[1] + 1, pos2.v[1] + 1);

		Vec pos[2] = { pos1, pos2 };

		u32 diff = 1;
				
		for (u8 d = 1;; ++d) {
			if (true && debug_board) {
				BitBoard opened1;
				for (u32 y = 0; y < 16; ++y) {
					opened1.set_row(y, opened.get_row(y) >> 16);
				}
				BitBoard* boards[3] = { &closed, &opened, &opened1 };
				TPixel colors[3] = { tigrRGB(100, 100, 100), tigrRGB(200, 100, 100), tigrRGB(100, 100, 200) };

				render_board(3, boards, colors, 0, 0);
			}

			if (diff) {
				diff = expand(begin, end, d);
			}

			u32 end_state = 0;

			for (u32 i = 0; i < 2; ++i) {
				if (min[i] < 0) {
					BoardMoves moves = NONE;

					Vec candidates[] = {
						Vec(pos[i].x() + 1, pos[i].y()),
						Vec(pos[i].x(), pos[i].y() + 1),
						Vec(pos[i].x() - 1, pos[i].y()),
						Vec(pos[i].x(), pos[i].y() - 1),
					};

					for (u32 j = 0; j < 4; ++j) {
						if (candidates[j].x() < 16
						 && !closed.get(candidates[j])
						 && !opened.get(Vec(candidates[j].x() + (i == 0 ? 16 : 0), candidates[j].y()))) {
							moves = moves | BoardMoves(1 << j);
						}
					}

					if (moves != NONE) {
						while (true) {
							u32 move = rng.next() >> (32 - 2);
							if (moves & (1 << move)) {
								pos[i] = candidates[move];
								closed.set(pos[i]);
								break;
							}
						}
					} else {
						end_state |= 1 << i;
					}
				}
			}

			if (pos[0] == pos[1]) {
				end_state = 3;
			}

			for (u32 i = 0; i < 2; ++i) {
				if (min[i] < 0 && (end_state >> i) & 1) {
					// TODO: We can fill 'closed' with this players 'opened' to terminate BFS earlier
					min[i] = (i32)d - 1;
				}
			}

			if (!diff && min[0] >= 0 && min[1] >= 0) {
				break;
			}
		}

		auto c = patterned ? Vec(count_pat(opened, 0, pos1, begin, end), count_pat(opened, 1, pos2, begin, end)) : count3(opened, begin, end);
		score[0] = c.v[0];
		score[1] = c.v[1];
	}
};

struct D {
	static u32 const max_cells_per_dist = 128; // This is probably higher than necessary
	u32 coord[2][max_cells_per_dist];
	u32 next_count;
	
	BitBoard closed;
	BitBoard opened[2];
	u32 score[2];

	D(BitBoard const& b)
		: closed(b) {
		opened[0].clear();
		opened[1].clear();
		score[0] = 0;
		score[1] = 0;
	}

	void add_node_no_check(u32 dist, u32 p, Vec pos) {
		coord[dist & 1][next_count++] = (pos.x() + (pos.y() << 4)) + (p << 8);
		opened[p].set(pos);
		score[p] += 1;
	}

	void add_node(u32 dist, u32 p, Vec pos) {
		if (!closed.get(pos) && !opened[p].get(pos)) {
			add_node_no_check(dist, p, pos);
		}
	}

	void expand_node(u32 base_dist, u32 p, Vec pos) {
		u32 dist = base_dist + 1;
		add_node(dist, p, Vec(pos.x() - 1, pos.y()));
		add_node(dist, p, Vec(pos.x() + 1, pos.y()));
		add_node(dist, p, Vec(pos.x(), pos.y() - 1));
		add_node(dist, p, Vec(pos.x(), pos.y() + 1));
	}

	void estimate(Vec pos1, Vec pos2) {
		next_count = 0;

		add_node_no_check(0, 0, pos1);
		add_node_no_check(0, 1, pos2);

		for (u32 dist = 0; dist < 256; ++dist) {
			u32 c = next_count;
			if (!c) {
				break;
			}

			// Clear next distance slot
			next_count = 0;

			for (u32 i = 0; i < c; ++i) {
				u16 cell = coord[dist & 1][i];
				closed.set(Vec(cell & 0xf, (cell >> 4) & 0xf));
			}

			if (false && debug_board) {
				BitBoard* boards[3] = { &closed, opened + 0, opened + 1 };
				TPixel colors[3] = { tigrRGB(100, 100, 100), tigrRGB(200, 100, 100), tigrRGB(100, 100, 200) };
				render_board(3, boards, colors, 0, 0);
			}

			for (u32 i = 0; i < c; ++i) {
				u16 cell = coord[dist & 1][i];
				expand_node(dist, cell >> 8, Vec(cell & 0xf, (cell >> 4) & 0xf));
			}
			
			if (false && debug_board) {
				BitBoard* boards[3] = { &closed, opened + 0, opened + 1 };
				TPixel colors[3] = { tigrRGB(100, 100, 100), tigrRGB(200, 100, 100), tigrRGB(100, 100, 200) };
				render_board(3, boards, colors, 0, 0);
			}
		}
	}
};

i32 score_board(Board2& board) {
	D d(board);
	d.estimate(board.headings[0].pos, board.headings[1].pos);

	return d.score[0] - d.score[1];
}

i32 score_board2(Board2& board, Vec pos, Player pl) {
	D2 d(board);
	d.estimate(pl == 0 ? pos : board.headings[0].pos, pl == 1 ? pos : board.headings[1].pos);

	return d.score[pl] - d.score[pl ^ 1];
}

i32 score_board3(Board2& board, Vec pos, Player pl, bool patterned) {
	D3 d(board);
	d.estimate(pl == 0 ? pos : board.headings[0].pos, pl == 1 ? pos : board.headings[1].pos, patterned);

	return d.score[pl] - d.score[pl ^ 1];
}

void score_board4(Board2& board, bool patterned, LcgPair& rng, MinMax<> (&result)[2]) {
	D4<false> d(board);
	d.estimate(board.headings[0].pos, board.headings[1].pos, patterned, rng);

	result[0] = MinMax<>(d.min[0], d.score[0]);
	result[1] = MinMax<>(d.min[1], d.score[1]);
}

BoardMoves best_dist_moves(Board2& board, Player pl, LcgPair& rng) {
	auto pos = board.headings[pl].pos;
	i32 best = -0x7fffffff;
	BoardMoves best_moves = NONE;
	Vec candidates[] = {
		Vec(pos.x() + 1, pos.y()),
		Vec(pos.x(), pos.y() + 1),
		Vec(pos.x() - 1, pos.y()),
		Vec(pos.x(), pos.y() - 1),
	};

	for (u32 i = 0; i < 4; ++i) {
		if (!board.get(candidates[i])) {
			Board2 cur_board(board);
			cur_board.set(candidates[i]);

			Vec candidates2[] = {
				Vec(candidates[i].x() + 1, candidates[i].y()),
				Vec(candidates[i].x(), candidates[i].y() + 1),
				Vec(candidates[i].x() - 1, candidates[i].y()),
				Vec(candidates[i].x(), candidates[i].y() - 1),
			};

			for (u32 j = 0; j < 4; ++j) {
				if (!board.get(candidates2[j])) {
					Board2 cur_board2(cur_board);
					cur_board2.set(candidates2[j]);
					cur_board2.headings[pl].pos = candidates2[j];

					MinMax<> result[2];
					score_board4(cur_board, false, rng, result);
					i32 cand_best = (i32)result[0].max - (i32)result[1].max;
					if (cand_best > best) {
						best = cand_best;
						best_moves = BoardMoves(1 << i);
					} else if (cand_best == best) {
						best_moves = best_moves | BoardMoves(1 << i);
					}
				}
			}
		}
	}

	return best_moves;
}
