#ifndef BOARD_h
#define BOARD_h

#include <sstream>
#include <cassert>
#include "game_enums.hpp"
#include "rand.hpp"
#include <cstdlib>
#ifdef _WIN32
#include <nmmintrin.h>
#endif

void init_board_tables();

struct PlayerHeading {
	PlayerHeading() {
	}

	PlayerHeading(Vec pos, BoardMoves prev_move)
		: pos(pos), prev_move(prev_move) {
	}

	static Vec next_pos(Vec cur_pos, BoardMoves new_move) {
		Vec new_pos = cur_pos;
		switch (new_move) {
		case UP: --new_pos.v[1]; break;
		case DOWN: ++new_pos.v[1]; break;
		case LEFT: --new_pos.v[0]; break;
		case RIGHT: ++new_pos.v[0]; break;
		}

		return new_pos;
	}

	Vec make_move(BoardMoves new_move) {
		Vec new_pos = next_pos(this->pos, new_move);

		this->pos = new_pos;
		this->prev_move = new_move;
		return new_pos;
	}

	BoardMoves get_move(u32 child) {
		return move_from_child(this->prev_move, child);
	}

	Vec pos;
	BoardMoves prev_move;
};

struct ChildIndexJ : Vec {
	using Vec::Vec;

	u32 index() const {
		return y() * 3 + x();
	}
};

struct PlayerPositions {
	PlayerHeading headings[2];

	PlayerPositions() {
		headings[0] = PlayerHeading(Vec(3, 7), RIGHT);
		headings[1] = PlayerHeading(Vec(12, 7), LEFT);
	}

	BoardMovesJ get_moves(ChildIndexJ child) {
		auto next_move0 = headings[0].get_move(child.v[0]);
		auto next_move1 = headings[1].get_move(child.v[1]);
		return BoardMovesJ(next_move0, next_move1);
	}

	BoardMoves move_from(u32 player, PlayerPositions const& prev_board) {
		auto prev_pos = prev_board.headings[player].pos;
		auto cur_pos = this->headings[player].pos;
		if (cur_pos.x() == prev_pos.x() + 1 && cur_pos.y() == prev_pos.y())
			return RIGHT;
		if (cur_pos.x() == prev_pos.x() && cur_pos.y() == prev_pos.y() + 1)
			return DOWN;
		if (cur_pos.x() == prev_pos.x() - 1 && cur_pos.y() == prev_pos.y())
			return LEFT;
		if (cur_pos.x() == prev_pos.x() && cur_pos.y() == prev_pos.y() - 1)
			return UP;
		assert(cur_pos == prev_pos);
		return this->headings[player].prev_move;
	}
};

template<typename DerivedT>
struct BoardCommon : PlayerPositions {
	
	BoardCommon() {
	}

	void print() {
		DerivedT& derived = *static_cast<DerivedT*>(this);

		system("cls");
		for (u32 y = 0; y < 16; ++y) {
			for (u32 x = 0; x < 16; ++x) {
				auto p = Vec(x, y);
				if (player_pos[0] == p) {
					printf("1");
				} else if (player_pos[1] == p) {
					printf("2");
				} else {
					printf(derived.is_wall(p) ? "." : " ");
				}
			}

			printf("|\n");
		}

		for (u32 x = 0; x < 16; ++x) {
			printf("=");
		}
	}

	BoardMoves legal_moves(Vec pos) {
		DerivedT& derived = *static_cast<DerivedT*>(this);

		BoardMoves moves = NONE;
		if (!derived.is_wall(Vec(pos.x, pos.y - 1))) moves = moves | UP;
		if (!derived.is_wall(Vec(pos.x, pos.y + 1))) moves = moves | DOWN;
		if (!derived.is_wall(Vec(pos.x - 1, pos.y))) moves = moves | LEFT;
		if (!derived.is_wall(Vec(pos.x + 1, pos.y))) moves = moves | RIGHT;
		return moves;
	}

	BoardMoves legal_moves(Player pl) {
		DerivedT& derived = *static_cast<DerivedT*>(this);

		auto pos = this->headings[pl].pos;
		return legal_moves(pos);
	}

	BoardMoves legal_moves2(Player pl) {
		DerivedT& derived = *static_cast<DerivedT*>(this);

		auto pos = this->player_pos[pl];
		BoardMoves moves = NONE;
		if (!derived.is_wall(Vec(pos.x, pos.y - 1)) && legal_moves(Vec(pos.x, pos.y - 1)) != NONE) moves = moves | UP;
		if (!derived.is_wall(Vec(pos.x, pos.y + 1)) && legal_moves(Vec(pos.x, pos.y + 1)) != NONE) moves = moves | DOWN;
		if (!derived.is_wall(Vec(pos.x - 1, pos.y)) && legal_moves(Vec(pos.x - 1, pos.y)) != NONE) moves = moves | LEFT;
		if (!derived.is_wall(Vec(pos.x + 1, pos.y)) && legal_moves(Vec(pos.x + 1, pos.y)) != NONE) moves = moves | RIGHT;
		return moves;
	}

	BoardMoves best_fill_moves(Player pl) {
		DerivedT& derived = *static_cast<DerivedT*>(this);

		auto pos = this->player_pos[pl];
		u32 best = 0;
		BoardMoves best_moves = NONE;
		Vec candidates[] = {
			Vec(pos.x + 1, pos.y),
			Vec(pos.x, pos.y + 1),
			Vec(pos.x - 1, pos.y),
			Vec(pos.x, pos.y - 1),
		};

		for (u32 i = 0; i < 4; ++i) {
			u32 cand_best = derived.flood_acc(candidates[i], 15); // flood
			if (cand_best > best) {
				best = cand_best;
				best_moves = BoardMoves(1 << i);
			} else if (cand_best == best) {
				best_moves = best_moves | BoardMoves(1 << i);
			}
		}

		return best_moves;
	}

	bool is_outside(Vec pos) {
		return pos.x >= 16 || pos.y >= 16;
	}

	void read(std::stringstream &stream) {
		DerivedT& derived = *static_cast<DerivedT*>(this);

		derived.clear();

		u32 x = 0, y = 0;
		std::string line;
		while (std::getline(stream, line, ',')) {
			if (line == ".") derived.set(Vec(x, y));

			if (line == "0") {
				this->headings[Pl1].pos = Vec(x, y);
			}
			if (line == "1")
				this->headings[Pl1].pos = Vec(x, y);
			x = (x + 1) % 16;
			if (x == 0)
				y++;
		}
	}
};

inline u32 popcount(u32 v) {
#if _WIN32
	return _mm_popcnt_u32(v);
#else
	u32 c = 0;
	for (; v; c++) {
		v &= v - 1;
	}
	return c;
#endif
}

inline u32 popcount16(u32 v) {
	return popcount(v);
}

inline u32 hflood(u32 seed, u32 mask, u32 w) {
	u32 imask = ~mask;

	for (u32 i = 0; i < w; ++i) {
		seed &= imask;
		seed |= (seed << 1) | (seed >> 1);
	}

	return seed & imask;
}

inline u32 rotl(u32 v, u32 c) {
	c &= 31;
	return (v << c) | (v >> (32 - c));
}

inline u32 rotr(u32 v, u32 c) {
	c &= 31;
	return (v >> c) | (v << (32 - c));
}

struct BitBoard {
	u32 bits[1 + 16 + 1];

	void clear_no_walls() {
		memset(bits, 0, sizeof(bits));
	}

	void set(Vec pos) {
		assert(pos.y() + 1 < 18);
		bits[pos.y() + 1] |= (u32(1) << (pos.x() & 31));
	}

	bool get(Vec pos) const {
		assert(pos.y() + 1 < 18);
		return (bits[pos.y() + 1] >> (pos.x() & 31)) & 1;
	}

	u32 get_row(u32 y) const {
		return bits[y + 1];
	}

	void set_row(u32 y, u32 new_row) {
		assert(y + 1 < 18);
		bits[y + 1] = new_row;
	}

	void or_row(u32 y, u32 set_bits) {
		assert(y + 1 < 18);
		bits[y + 1] |= set_bits;
	}

	bool is_wall(Vec pos) const {
		return this->get(pos);
	}

	void clear() {
		bits[0] = 0xffffffff;
		bits[1 + 16] = 0xffffffff;

		for (u32 i = 1; i < 1 + 16; ++i) {
			bits[i] = 0xffff0000;
		}
	}

	void dup() {
		for (u32 i = 1; i < 1 + 16; ++i) {
			bits[i] = (bits[i] & 0xffff) | (bits[i] << 16);
		}
	}
};

template<typename T>
struct ShiftBoard {
	ShiftBoard(T const& a, i32 shift) : a(a), shift(shift) {
	}

	u32 get_row(u32 y) const {
		return (a.get_row(y) >> shift) | 0xffff0000;
	}

	T const& a;
	i32 shift;
};

template<typename A, typename B>
struct OrBoard {
	OrBoard(A const& a, B const& b) : a(a), b(b) {
	}

	u32 get_row(u32 y) const {
		return a.get_row(y) | b.get_row(y);
	}

	A const& a;
	B const& b;
};

template<typename A, typename B>
inline OrBoard<A, B> or_board(A const& a, B const& b) {
	return OrBoard<A, B>(a, b);
}

template<typename T>
inline ShiftBoard<T> shift_board(T const& a, i32 shift) {
	return ShiftBoard<T>(a, shift);
}

//  0
// 321
//76 54
// a98
//  b

template<typename B>
inline u32 diamond5(B const& board, Vec pos) {
	assert(pos.y() < 16);
	u32 r0 = pos.y() == 0 ? 0xffffffff : board.get_row(pos.y() - 2);
	u32 r1 = board.get_row(pos.y() - 1);
	u32 r2 = board.get_row(pos.y());
	u32 r3 = board.get_row(pos.y() + 1);
	u32 r4 = pos.y() == 15 ? 0xffffffff : board.get_row(pos.y() + 2);

	u32 d = (r0 >> pos.x()) & 1;
	d |= rotr(r1, pos.x() - 2) & (7 << 1);
	d |= rotr(r2, pos.x() - 6) & (3 << 4);
	d |= rotr(r2, pos.x() - 5) & (3 << 6);
	d |= rotr(r3, pos.x() - 9) & (7 << 8);
	d |= rotr(r4, pos.x() - 11) & (1 << 11);

	return d;
}

inline u32 index_diamond5(u32 d5, u32 first_dir) {
	switch (first_dir) {
	case RIGHT:
		return ((d5 >> 6) & 1) | (((d5 >> 10) & 1) << 1) | (((d5 >> 7) & 1) << 2) | (((d5 >> 3) & 1) << 3);
	case DOWN:
		return ((d5 >> 9) & 1) | (((d5 >> 8) & 1) << 1) | (((d5 >> 11) & 1) << 2) | (((d5 >> 10) & 1) << 3);
	case LEFT:
		return ((d5 >> 5) & 1) | (((d5 >> 1) & 1) << 1) | (((d5 >> 4) & 1) << 2) | (((d5 >> 8) & 1) << 3);
	case UP:
		return ((d5 >> 2) & 1) | (((d5 >> 3) & 1) << 1) | (((d5 >> 0) & 1) << 2) | (((d5 >> 1) & 1) << 3);
	default:
		return 0;
	}
}

struct Board2 : BoardCommon<Board2>, BitBoard {

	void compact(Vec pos, Vec opp_pos, u32 (&dest)[32], u32(&opp)[32]) {
		for (u32 y = 0; y < 32; ++y) {
			if (y + pos.y() < 15 || y + pos.y() > 15+15) {
				dest[y] = 0xffffffff;
			} else {
				dest[y] = rotr(bits[y + pos.y() - 15 + 1], pos.x() - 15);
			}

			if (y + pos.y() - 15 != opp_pos.y()) {
				opp[y] = 0;
			} else {
				opp[y] = rotr((1 << opp_pos.x()), pos.x() - 15);
			}
		}
	}

	Board2();
	Board2(std::stringstream &stream);
};

template<typename T = u8>
struct MinMax {
	T min, max;

	MinMax()
		: MinMax(0, 254) {
	}

	MinMax(T min, T max)
		: min(min), max(max) {
		assert(this->is_valid());
	}

	static MinMax zero() {
		return MinMax(0, 0);
	}

	static MinMax at_least(T m) {
		return MinMax(m, 254);
	}

	static MinMax at_most(T m) {
		return MinMax(0, m);
	}

	MinMax operator+(T x) const {
		MinMax r(this->min + x, this->max + x);
		return r;
	}

	MinMax operator-(T x) const {
		assert(this->min >= x && this->max >= x);
		return MinMax(this->min - x, this->max - x);
	}

	MinMax operator&(MinMax other) const {
		if (this->max <= other.min || other.max <= this->min) {
			// Previous values assumed no stupid moves were made, and
			// clearly we made stupid moves :(. Assume other is more accurate.
			return other;
		}
		return MinMax(std::max(this->min, other.min), std::min(this->max, other.max));
	}

	MinMax operator|(MinMax other) const {
		return MinMax(std::min(this->min, other.min), std::max(this->max, other.max));
	}

	MinMax operator^(MinMax other) const {
		// [1, 1] [1, 2]
		// 1 <= 2 -> true
		// 2 <= 1 -> false
		if (this->max <= other.min) {
			assert(other.is_valid());
			return other;
		} else if (other.max <= this->min) {
			assert(this->is_valid());
			return *this;
		}

		return MinMax(std::min(this->min, other.min), std::max(this->max, other.max));
	}

	bool and_(MinMax other) {
		bool changed = false;
		if (this->max <= other.min || other.max <= this->min) {
			// Previous values assumed no stupid moves were made, and
			// clearly someone made stupid moves. Assume other is more accurate.
			*this = other;
			return true;
		}

		if (other.min > this->min) {
			this->min = other.min;
			changed = true;
		}
		if (other.max < this->max) {
			this->max = other.max;
			changed = true;
		}
		assert(this->is_valid());

		return changed;
	}

	bool operator<(T other_min) const {
		return this->max < other_min || (this->max == other_min && this->min < other_min);
	}

	bool operator<(MinMax& other) const {
		return this->max < other.min || (this->max == other.min && this->min < other.min);
	}

	bool operator==(MinMax const& other) const {
		return this->min == other.min && this->max == other.max;
	}

	bool operator!=(MinMax const& other) const {
		return !operator==(other);
	}

	bool is_valid() const {
		return this->min <= this->max;
	}
};

i32 score_board(Board2& board);
i32 score_board2(Board2& board, Vec pos, Player pl);
i32 score_board3(Board2& board, Vec pos, Player pl, bool patterned);
void score_board4(Board2& board, bool patterned, LcgPair& rng, MinMax<>(&result)[2]);
void score_board_just_max(Board2& board, bool patterned, LcgPair& rng, u32(&result)[2]);

BoardMoves best_dist_moves(Board2& board, Player pl, LcgPair& rng);
BoardMoves best_lookup_moves(Board2& board, Player pl);

typedef Board2 Board;

struct GameState {
	Board board;

	GameState() {
	}

	static GameState standard() {
		GameState state;

		state.board.set(state.board.headings[0].pos);
		state.board.set(state.board.headings[1].pos);

		return state;
	}

	static GameState random(LcgPair& rng) {
		GameState state;

		u32 x = rng.get_u32(7), y = rng.get_u32(15);

		state.board.headings[0].pos = Vec(x, y);
		state.board.headings[1].pos = Vec(15 - x, 15 - y);

		state.board.set(state.board.headings[0].pos);
		state.board.set(state.board.headings[1].pos);

		return state;
	}

	GameState(Board const& board_init)
		: board(board_init) {
	}

	BoardMoves pick_move(LcgPair& rng, u32 player_id, BoardMoves moves) {
		if (moves == NONE) {
			return this->board.headings[player_id].prev_move;
		} else {
			while (true) {
				u32 move = rng.next() >> (32 - 2);
				if (moves & (1 << move)) {
					return BoardMoves(1 << move);
				}
			}
		}
	}

	u32 update(BoardMovesJ new_moves) {
		u32 hit = 0;

		for (int i = 0; i < 2; ++i) {
			auto new_move = new_moves.p[i];

			auto new_pos = board.headings[i].make_move(new_move);

			if (board.is_wall(new_pos)) {
				hit |= u32(1) << i;
			} else {
				board.set(new_pos);
			}
		}

		if (board.headings[0].pos == board.headings[1].pos) {
			hit = (u32(1) << 0) | (u32(1) << 1);
		}

		return hit;
	}
};

#endif
