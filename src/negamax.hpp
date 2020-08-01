#include <vector>
#include <algorithm>
#include "game_enums.hpp"
#include "rand.hpp"
#include <intrin.h>
#pragma intrinsic(_BitScanForward, _BitScanReverse)

using std::vector;
using std::min;
using std::max;

namespace negamax {

typedef u8 CellCoord;

struct Move {
	u16 v;

	Move() : v(0) {
	}

	Move(CellCoord worker_1, CellCoord worker_2)
		: v(0x8000 | ((u16)worker_2 << 5) | ((u16)worker_1)) {
	}

	Move(CellCoord from, CellCoord move_to, CellCoord build_on)
		: v(((u16)from << 10) | ((u16)build_on << 5) | ((u16)move_to)) {
	}

	bool is_initial() const {
		return this->v >= 0x8000;
	}

	CellCoord initial_worker_1() const {
		return (u8)(this->v & 31);
	}

	CellCoord initial_worker_2() const {
		return (u8)(this->v >> 5);
	}

	CellCoord from() const {
		return (u8)(this->v >> 10);
	}

	CellCoord move_to() const {
		return (u8)(this->v & 31);
	}

	CellCoord build_on() const {
		return (u8)((this->v >> 5) & 31);
	}
};

struct ZobHash {
	u64 levels[5 * 5][5];
	u64 player[5 * 5][2];

	ZobHash() {
		LcgPair rand = LcgPair::from_bits(1, 2);

		for (int c = 0; c < 5*5; ++c) {
			for (int i = 0; i < 2; ++i) {
				this->player[c][i] = rand.next_u64() & ~1;
			}

			for (int i = 0; i < 5; ++i) {
				this->levels[c][i] = rand.next_u64() & ~1;
			}
		}
	}
};

#define MAX_MOVES ((5 * 5) * (5 * 5 - 1) / 2)

extern u32 ADJACENT_CELLS[5 * 5];

void init();

inline u32 trailing_zeros(u32 x) {
	unsigned long v = 32;
	_BitScanForward(&v, x);
	return v;
}

inline u32 leading_zeros(u32 x) {
	unsigned long v = 0;
	_BitScanReverse(&v, x);
	return 31 - v;
}

struct Board {
	u32 cells_h0_, cells_h1_;
	u32 cells_b_[2];
	u64 hash;

	Board()
		: cells_h0_(0), cells_h1_(0)
		, cells_b_{0, 0}, hash(0) {
	}

	u8 current_player() {
		return u8(this->hash & 1);
	}

	u32 cells_b() const { return cells_b_[0] | cells_b_[1]; }

	bool is_blocked(CellCoord coord) const {
		return (this->cells_b() >> coord) & 1 != 0;
	}

	bool is_third(CellCoord coord) const {
		return ((this->cells_h0_ & this->cells_h1_) >> coord) & 1 != 0;
	}

	usize level_assume_0_to_3(CellCoord coord) const {
		auto l = (this->cells_h1_ >> coord << 1) & 2;
		l |= (this->cells_h0_ >> coord) & 1;
		return (usize)l;
	}

	void build(CellCoord coord, ZobHash& h) {

		auto mask = 1 << coord;
		auto pl = this->level_assume_0_to_3(coord);

		//let pcells_b_ = self.cells_b_ & mask;
		auto overflow = mask & this->cells_h1_ & this->cells_h0_;

		this->cells_b_[0] |= overflow; // Overflow to blocked
		this->cells_b_[1] |= overflow; // Overflow to blocked
		this->cells_h1_ ^= mask & this->cells_h0_;
		this->cells_h0_ ^= mask;

		//let ll = self.level(coord);
		auto l = pl + 1;

		/*
		if ll as usize != l {
		println!("building to {}. {}  at {}. prev l: {}", l, ll, coord.0, pl);
		panic!();
		}*/

		this->hash ^= h.levels[coord][l];
		//println!("hash level {} on pos {} -> {:x}", l, coord.0, self.hash);
	}

	void perform_move(Move m, ZobHash& h) {
        auto player = (usize)this->current_player();

        if (m.is_initial()) {
            
            auto w1 = m.initial_worker_1();
			auto w2 = m.initial_worker_2();
            this->cells_b_[player] = (1 << w1) | (1 << w2);

			this->hash ^= h.player[w1][player] ^ h.player[w2][player];
        } else {
            auto worker_pos = m.from();
            auto move_dest = m.move_to();

			this->cells_b_[player] ^= 1 << worker_pos;
			this->cells_b_[player] |= 1 << move_dest;

			this->hash ^= h.player[worker_pos][player] ^ h.player[move_dest][player];

            // No build if stepping to level 3
            if (!this->is_third(move_dest)) {
                auto build_dest = m.build_on();
				this->build(build_dest, h);
            }
        }

		this->hash ^= 1;
    }

	struct CoordPair { CellCoord c[2]; };

	u32 worker_mask(usize player) {
		return this->cells_b_[player] & ~this->cells_b_[1 - player];
	}

	inline u32 trailing_zeros(u32 x) {
		unsigned long v = 32;
		_BitScanForward(&v, x);
		return v;
	}

	CoordPair worker_pos(usize player) {
        auto x = this->worker_mask(player);
		CoordPair pair;
		pair.c[0] = trailing_zeros(x);
		pair.c[1] = 31 - leading_zeros(x);
		return pair;
    }

	usize valid_moves(Move* moves) {
		usize count = 0;
		if (this->cells_b_[1] == 0) {
			for (auto w1 = 0; w1 < 5 * 5 - 1; ++w1) {
				//let w1 = CellCoord::new(w1x, w1y);
				if (!this->is_blocked(w1)) {
					for (auto w2 = w1 + 1; w2 < 5*5; ++w2) {
						if (!this->is_blocked(w2)) {
							moves[count++] = Move(w1, w2);
						}
					}
				}
			}
		} else {
			auto player = (usize)this->current_player();
			auto other_player = 1 - player;

			//let mut opp_can_move_to_third = 0;
			auto cells_h0 = this->cells_h0_;
			auto cells_h1 = this->cells_h1_;
			auto cells_b = this->cells_b();

			/*
			for opp_worker in 0..2 {
			// Check opponent winning moves
			let opp_worker_pos = $board.worker_coord[other_player][opp_worker];

			//let from_cell_h1 = ((cells_h1 >> opp_worker_pos.0) & 1).wrapping_neg();
			//let within_level_and_third = from_cell_h1 & cells_h0 & cells_h1;
			//let opp_move_mask = reach << (opp_worker_pos.0 - 1*6 - 1);
			//opp_can_move_to_third |= within_level_and_third & opp_move_mask;

			if $board.is_third(opp_worker_pos) {
			// Opponent is already on third level, we lose
			$no_moves;
			}
			}*/

			auto thirds = this->cells_h0_ & this->cells_h1_;

			if ((this->cells_b_[other_player] & thirds) != 0) {
				return 0;
			}

			u32 can_move_tos[2] = { 0, 0 };
			auto worker_poses = this->worker_pos(player);

			for (int worker = 0; worker < 2; ++worker) {
				auto worker_pos = worker_poses.c[worker];
				auto from_cell_h0 = 0 - ((cells_h0 >> worker_pos) & 1);
				auto from_cell_h1 = 0 - ((cells_h1 >> worker_pos) & 1);

				auto without_level = ~from_cell_h1 & cells_h1 & (~from_cell_h0 | cells_h0);
				auto move_mask = ADJACENT_CELLS[worker_pos]; //reach << (worker_pos.0 - 1*6 - 1);
				auto can_move_to = ~(cells_b | without_level) & move_mask;

				auto win_moves = can_move_to & cells_h0 & cells_h1;

				if (win_moves != 0) {
					auto win_move = trailing_zeros(win_moves);
					moves[0] = Move(worker_pos, win_move, 0);
					return 1;
				}

				can_move_tos[worker] = can_move_to;
			}


			for (int worker = 0; worker < 2; ++worker) {

				auto worker_pos = worker_poses.c[worker];
				auto can_move_to = can_move_tos[worker];

				auto without_original = cells_b ^ (1 << worker_pos);
				auto can_build_on_all = ~without_original;

				while (can_move_to != 0) {
					auto move_to = trailing_zeros(can_move_to);
					can_move_to ^= 1 << move_to;

					auto build_mask = ADJACENT_CELLS[move_to]; //reach << (move_to - 1*6 - 1);
					auto can_build_on = can_build_on_all & build_mask;

					/*
					if $opp_win_check && opp_can_move_to_third != 0 {
					// If opponent can move to a third level, only consider building domes
					// on those.
					can_build_on &= opp_can_move_to_third;
					}
					*/

					while (can_build_on != 0) {
						auto build_on = trailing_zeros(can_build_on);
						can_build_on ^= 1 << build_on;

						moves[count++] = Move(worker_pos, move_to, build_on);
					}
				}
			}
		}

		return count;
	}
};

#define TT_TABLE_MASK (0xffffff)

struct TtEntry {
	u64 hash;
	i32 value;
	Move best_move;
	u8 depth;
	i8 flags;

	TtEntry()
		: hash(1), value(0), best_move(Move()), depth(0), flags(0) {
	}
};


struct Negamax {
		
	vector<TtEntry> tt;
	usize total_evals;
	ZobHash zob_hash;

	Negamax()
		: tt(TT_TABLE_MASK + 1, TtEntry()),
		  total_evals(0) {
	}

	void save(u64 h, u8 depth, i32 value, i8 flags, Move best_move) {
		auto& entry = this->tt[h & TT_TABLE_MASK];

		if (entry.hash == h && entry.depth > depth) { return; }

		entry.hash = h;
		entry.value = value;
		entry.flags = flags;
		entry.best_move = best_move;
		entry.depth = depth;
	}

	i32 negamax_3(Board board, u32 depth, i32 a, i32 b) {
		return this->negamax_3_(board, depth, a, b);
	}

	i32 negamax_3_(Board board, u32 depth, i32 a, i32 b) {
		if (depth == 0) {
			this->total_evals += 1;
			return 0;
		}

		auto alpha_orig = a;
        Move best_move;

        auto h = board.hash;
		auto entry = this->tt[h & TT_TABLE_MASK];
        if (entry.hash == h && entry.depth >= (u8)depth) {
            best_move = entry.best_move;
            if (entry.flags == 0) {
               return entry.value;
            } else if (entry.flags == 1) {
                a = max(a, entry.value);
            } else if (entry.flags == -1) {
                b = min(b, entry.value);
            }

            if (a >= b) {
                return entry.value;
            }
        } else {
            best_move = Move();
        }

        auto value = -10000;
        
        if (best_move.v != 0) {
            auto next_board = board;
            next_board.perform_move(best_move, this->zob_hash);
            auto v = -this->negamax_3_(next_board, depth - 1, -b, -a);
            value = max(v, value);
            a = max(a, value);
        }

        auto new_best_move = best_move;

        if (a < b) {

            Move actions[MAX_MOVES];
            auto count = board.valid_moves(actions);
            if (count == 0) {
                this->total_evals += 1;
                return board.current_player() == 0 ? -10000 : 10000;
            }

			for (int i = 0; i < count; ++i) {
				auto action = actions[i];
				if (action.v != best_move.v) {
					auto next_board = board;
					next_board.perform_move(action, this->zob_hash);
					auto v = -this->negamax_3_(next_board, depth - 1, -b, -a);
					a = max(a, v);
                    if (v > value) {
                        new_best_move = action;
                        value = v;
                        if (a >= b) {
                            break;
                        }
                    }
				}
			}
        }

		i8 flag;
        if (value <= alpha_orig) {
			flag = 1;
        } else if (value >= b) {
			flag = -1;
        } else {
			flag = 0;
        }

        this->save(board.hash, (u8)depth, value, flag, new_best_move);

		return value;
	}
};

}
