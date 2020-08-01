#ifndef MCBOTJ_HPP
#define MCBOTJ_HPP 1

#include <utility>
#include <memory>
#include <vector>
#include <cmath>
#include <cassert>

#include "gametree.hpp"
#include "bot.hpp"
#include "board.hpp"
#include "rand.hpp"

static bool const force_first_moves = false;

struct CounterI {
	u32 wins, total;

	static CounterI even() {
		CounterI z;
		z.wins = 1;
		z.total = 2;
		return z;
	}

	static CounterI unlikely() {
		CounterI z;
		z.wins = 0;
		z.total = 0x80000000;
		return z;
	}

	static CounterI empty() {
		CounterI z;
		z.wins = 0;
		z.total = 0;
		return z;
	}

	CounterI normalize(u32 wanted_total) {
		if (!total) {
			return *this;
		} else {
			int scale = (total + wanted_total - 1) / wanted_total;
			CounterI c;
			c.wins = wins / scale;
			c.total = total / scale;
			return c;
		}
	}

	f64 ucb(double c, int variant) const {
		// TODO: Table of c * sqrt(log(total + 1) / total)
		double adjusted_total = total * 0.5;
		return (double)wins / total + c * sqrt(log(adjusted_total + 1.0) / adjusted_total);
	}

	f64 ratio() const {
		return total ? (double)wins / total : 0.0;
	}

	void add(u32 wins, u32 total) {
		this->wins += wins;
		this->total += total;
	}

	void add(u32 wins) {
		this->wins += wins;
		this->total += 2;
	}

	bool less_visited(CounterI const& other) const {
		u32 f_total = this->total & 0x7fffffff;
		u32 f_other_total = other.total & 0x7fffffff;
		return f_total < f_other_total || (f_total == f_other_total && this->wins < other.wins);
	}

	bool operator==(CounterI const& other) const {
		return this->total == other.total && this->wins == other.wins;
	}
};

struct CounterF {
	f32 r;
	u32 total;

	struct SelectInfo {};

	CounterF() {}

	CounterF(f32 r, u32 total) : r(r), total(total) {
	}

	static CounterF even() {
		return CounterF(0.5f, 1);
	}

	static CounterF unlikely() {
		return CounterF(0.f, 0x80000000);
	}

	static CounterF empty() {
		return CounterF(0.f, 0);
	}

	static CounterF even_zero() {
		return CounterF(0.5f, 0);
	}

	static f32 ucb(f32 c, f32 adjusted_total) {
		return c * sqrtf(logf(adjusted_total + 1.0f) / adjusted_total);
	}

	template<int Variant>
	f32 ucb(f32 c, u32 total_total) const {
		f32 adjusted_total = f32(total);

		if (false) {
			return r + 1.15f * sqrtf((f32)total_total) / (adjusted_total + 1.f);
			//return r + 1.0f * sqrtf(total_total) / (adjusted_total + 1.f);
		} else {
			if (total < 1024) {
				return r + ucb_tab[total];
			} else {
				return r + c * sqrtf(logf(adjusted_total + 1.0f) / adjusted_total);
			}
		}
	}

	f32 ratio() const {
		return r;
	}

	template<int Variant>
	void add(f32 wins, u32 pl, SelectInfo const&) {
		if (!this->is_dead()) { // TODO: Check how often this triggers and try to minimize it
			auto prev_total = f32(this->total);

			this->r = this->r + (wins - this->r) / (prev_total + 1.f);
			++this->total;
		}
	}

	bool less_visited(CounterF const& other) const {
		u32 f_total = this->total & 0x7fffffff;
		u32 f_other_total = other.total & 0x7fffffff;
		return f_total < f_other_total || (f_total == f_other_total && this->r < other.r);
	}

	bool is_dead() const {
		return this->total == 0x80000000;
	}

	bool operator==(CounterF const& other) const {
		return this->total == other.total && this->r == other.r;
	}
};

typedef f64 WeightT;
static WeightT const gamma = 0.07;
static WeightT const neg_gamma = 1.0 - gamma;
static WeightT const adjusted_gamma = gamma / 3.0;
static WeightT const inv_max_int = 1.0 / 4294967296.0;
static WeightT const log_one_m_gamma = log(1.0 - gamma);

struct CounterExp3 {
	struct SelectInfo {
		WeightT dist[2];

		SelectInfo() {
		}

		SelectInfo(WeightT dist0, WeightT dist1) {
			dist[0] = dist0;
			dist[1] = dist1;
		}
	};

	WeightT lw;

	CounterExp3() {
	}

	CounterExp3(WeightT lw) : lw(lw) {
	}

	static CounterExp3 unlikely() {
		return CounterExp3(-100.0f);
	}

	static CounterExp3 even_zero() {
		return CounterExp3(0.0f);
	}

	template<int Variant>
	void add(f32 r, u32 pl, SelectInfo const& select_info) {
		WeightT estimated = r / select_info.dist[pl];
		auto le = estimated * adjusted_gamma;
		auto lnv = this->lw + le;
		this->lw = lnv;
	}

	static WeightT log_sum_of_exp(WeightT lw0, WeightT lw1, WeightT lw2) {
		if (lw1 >= lw0) {
			if (lw2 >= lw1) {
				return lw2 + log(exp(lw0 - lw2) + exp(lw1 - lw2) + 1.0f);
			} else {
				return lw1 + log(exp(lw0 - lw1) + 1.0 + exp(lw2 - lw1));
			}
		} else {
			if (lw2 >= lw0) {
				return lw2 + log(exp(lw0 - lw2) + exp(lw1 - lw2) + 1.0f);
			} else {
				return lw0 + log(1.0 + exp(lw1 - lw0) + exp(lw2 - lw0));
			}
		}
	}

	WeightT prob_(WeightT logsum) const {
		return exp(log_one_m_gamma + lw - logsum) + adjusted_gamma;
	}

	WeightT prob_2(WeightT logsum) const {
		return exp(log_one_m_gamma + lw - logsum) + adjusted_gamma;
	}

	/*
	WeightT prob_m_adjusted_gamma(WeightT logsum) const {
		//return expf(log_one_m_gamma + lw - logsum);

		return expf(log_one_m_gamma + lw) / sum;
	}*/

	f32 ratio() const {
		return lw;
	}
};

struct PlayerCounterExp3 {
	CounterExp3 move[3];

	template<int Variant>
	u32 pick_move(LcgPair& rng, WeightT& dist) {
		WeightT logsum = CounterExp3::log_sum_of_exp(move[0].lw, move[1].lw, move[2].lw);
		WeightT draw = WeightT(rng.next()) * inv_max_int;

		WeightT p0 = move[0].prob_(logsum) + adjusted_gamma;
		WeightT p1 = move[1].prob_(logsum) + adjusted_gamma;
		WeightT p2 = move[2].prob_(logsum) + adjusted_gamma;
		if (draw < p0) {
			dist = p0;
			return 0;
		}
		draw -= p0;
		
		if (draw < p1) {
			dist = p1;
			return 1;
		}

		dist = p2;
		return 2;
	}

	u32 pick_best(LcgPair& rng) {
		auto best = move[0].lw;
		u32 best_children[3];
		u32 best_children_count = 1;
		best_children[0] = 0;

		for (u32 m = 1; m < 3; ++m) {

			auto& cand_best = move[m];

			if (best < cand_best.lw) {
				best = cand_best.lw;
				best_children_count = 1;
				best_children[0] = m;
			} else if (best == cand_best.lw) {
				best_children[best_children_count++] = m;
			}
		}

		return best_children[rng.get_u32(best_children_count)];
	}
};

struct BaseGameNodeJ {
	static u32 const ChildCount = 3 * 3;
	NodeId children[ChildCount];

	BaseGameNodeJ() {
		memset(children, 0, sizeof(children));
	}

	static u32 child_index(BoardMovesJ prev_moves, BoardMovesJ new_moves) {
		return ChildIndexJ(
			child_from_move(prev_moves.p[0], new_moves.p[0]),
			child_from_move(prev_moves.p[1], new_moves.p[1])).index();
	}

	u32& get(ChildIndexJ index) {
		return this->children[index.index()];
	}
};

struct GameNodeJExp3 : BaseGameNodeJ {
	typedef CounterExp3::SelectInfo SelectInfo;
	typedef CounterExp3 Counter;

	PlayerCounterExp3 players[2];

	bool has_stats(ChildIndexJ index) const {
		return players[0].move[index.x()].lw != 0.f
			&& players[1].move[index.y()].lw != 0.f;
	}

	u32 pick_best(LcgPair& rng, u32 player_id) {
		return players[player_id].pick_best(rng);
	}

	template<int Variant>
	ChildIndexJ select(LcgPair& rng, SelectInfo& select_info) {

		int unvisited_count0 = (this->players[0].move[0].lw == 0.0) + (this->players[0].move[1].lw == 0.0) + (this->players[0].move[2].lw == 0.0);

		u32 p0, p1;
		if (unvisited_count0 > 0) {
			int pick = rng.get_u32(unvisited_count0);
			for (u32 c = 0; c < 3; ++c) {
				if (this->players[0].move[c].lw == 0.0) {
					if (pick == 0) {
						p0 = c;
						select_info.dist[0] = 1.0 / 3.0;
						break;
					}
					--pick;
				}
			}
		} else {
			auto& p = players[0];
			p0 = p.pick_move<Variant>(rng, select_info.dist[0]);
		}

		int unvisited_count1 = (this->players[1].move[0].lw == 0.0) + (this->players[1].move[1].lw == 0.0) + (this->players[1].move[2].lw == 0.0);
		if (unvisited_count1 > 0) {
			int pick = rng.get_u32(unvisited_count1);
			for (u32 c = 0; c < 3; ++c) {
				if (this->players[1].move[c].lw == 0.0) {
					if (pick == 0) {
						p1 = c;
						select_info.dist[1] = 1.0 / 3.0;
						break;
					}
					--pick;
				}
			}
		} else {
			auto& p = players[1];
			p1 = p.pick_move<Variant>(rng, select_info.dist[1]);
		}

		ChildIndexJ idx(p0, p1);

		return idx;
	}

	GameNodeJExp3() {
		memset(children, 0, sizeof(children));
	}
};

static f32 const c1 = f32(sqrt(2.0));
static f32 const c2 = f32(1.2);

extern f32 ucb_tab[1024];

void init_ucb_tab();

struct PlayerCounter {
	typedef CounterF Counter;

	Counter move[3];
	MinMax<> life[3];

	template<int Variant>
	u32 pick_move(LcgPair& rng) {
		u32 total_total = move[0].total + move[1].total + move[2].total;
		f32 m0 = move[0].ucb<Variant>(c2, total_total);
		f32 m1 = move[1].ucb<Variant>(c2, total_total);
		f32 m2 = move[2].ucb<Variant>(c2, total_total);
		
		if (m1 > m0) {
			// {1, ...}
			if (m2 > m1) {
				// {2}
				return 2;
			} else if (m2 < m1) {
				// {1}
				return 1;
			} else {
				// {1, 2}
				return 1 + (rng.next() >> 31);
			}
		} else if (m1 < m0) {
			// {0, ...}
			if (m2 > m0) {
				// {2}
				return 2;
			} else if (m2 < m0) {
				// {0}
				return 0;
			} else {
				// {0, 2}
				return (rng.next() >> 31) << 1;
			}
		} else {
			// {0, 1, ...}
			if (m2 > m1) {
				// {2}
				return 2;
			} else if (m2 < m1) {
				// {0, 1}
				return (rng.next() >> 31);
			} else {
				// {0, 1, 2}
				return rng.get_u32(3);
			}
		}
	}

	u32 pick_best(LcgPair& rng) {
		CounterF best = move[0];
		u32 best_children[3];
		u32 best_children_count = 1;
		best_children[0] = 0;

		for (u32 m = 1; m < 3; ++m) {

			auto& cand_best = move[m];

			if (best.less_visited(cand_best)) {
				best = cand_best;
				best_children_count = 1;
				best_children[0] = m;
			} else if (best == cand_best) {
				best_children[best_children_count++] = m;
			}
		}

		return best_children[rng.get_u32(best_children_count)];
	}
};

struct GameNodeJ : BaseGameNodeJ {
	typedef CounterF::SelectInfo SelectInfo;
	typedef PlayerCounter::Counter Counter;

	PlayerCounter players[2];
	MinMax<> life[ChildCount][2];
	
	bool has_stats(ChildIndexJ index) const {
		return players[0].move[index.x()].total != 0
			&& players[1].move[index.y()].total != 0;
	}

	u32 pick_best(LcgPair& rng, u32 player_id) {
		return players[player_id].pick_best(rng);
	}

	template<int Variant>
	ChildIndexJ select(LcgPair& rng, SelectInfo& select_info) {
		u32 p0, p1;
		int unvisited_count0 = (this->players[0].move[0].total == 0) + (this->players[0].move[1].total == 0) + (this->players[0].move[2].total == 0);

		if (unvisited_count0 > 0) {
			int pick = rng.get_u32(unvisited_count0);
			for (u32 c = 0; c < 3; ++c) {
				if (this->players[0].move[c].total == 0) {
					if (pick == 0) {
						p0 = c;
						break;
					}
					--pick;
				}
			}
		} else {
			auto& p = players[0];

			p0 = p.pick_move<Variant>(rng);
		}

		int unvisited_count1 = (this->players[1].move[0].total == 0) + (this->players[1].move[1].total == 0) + (this->players[1].move[2].total == 0);
		if (unvisited_count1 > 0) {
			int pick = rng.get_u32(unvisited_count1);
			for (u32 c = 0; c < 3; ++c) {
				if (this->players[1].move[c].total == 0) {
					if (pick == 0) {
						p1 = c;
						break;
					}
					--pick;
				}
			}
		} else {
			auto& p = players[1];
			p1 = p.pick_move<Variant>(rng);
		}

		ChildIndexJ idx(p0, p1);

		return idx;
	}

	GameNodeJ() {
		//memset(players, 0, sizeof(players));
	}
};


struct Writer {
	FILE* f;
	usize bytes;

	Writer(char const* path) {
		f = fopen(path, "a+");
		bytes = 0;
	}

	void write(void const* v, usize len) {
		fwrite(v, 1, len, f);
		bytes += len;
		if (bytes > (1 << 10)) {
			bytes &= (1 << 10) - 1;
		}
		fflush(f);
	}

	~Writer() {
		fclose(f);
	}
};

template<typename GN>
struct VisitedEdge : GN::SelectInfo {
	NodeId node;
	ChildIndexJ edge;
	u8 to_pos[2];

	VisitedEdge() {}

	VisitedEdge(NodeId node, ChildIndexJ edge, u8 to_pos0, u8 to_pos1, typename GN::SelectInfo const& select_info)
		: node(node), edge(edge), SelectInfo(select_info) {
		to_pos[0] = to_pos0;
		to_pos[1] = to_pos1;
	}

	bool operator==(VisitedEdge const& other) const {
		return node == other.node && edge == other.edge;
	}
};

template<int Variant, typename GN = GameNodeJ>
struct McBotJ : Bot {
	BasicGameTree<GN> tree;
	u32 sims, depth;
	Writer* writer;
	bool has_written_board;
	u32 updates;
	u32 cells_left;
	u32 pruned;

	CounterF amaf[2][256];
	vector<VisitedEdge<GN>> visited;
	usize total_updates, same_count;

	McBotJ(LcgPair& rng, Writer* writer = 0)
		: sims(0), depth(0), Bot(rng), writer(writer),
		  has_written_board(false),
		  total_updates(0),
		  same_count(0),
		  updates(0),
		  cells_left(0),
		  pruned(0) {

		init_ucb_tab();

		memset(amaf, 0, sizeof(amaf));
	}

	virtual BoardMoves move(int time) {
		if (deterministic) {
			for (u32 i = 0; i < iterations; ++i) {
				run();
			}
		} else {
			u32 start = ticks();
			time = time * 95 / 100; // Margin of error
			while (true) {
				for (u32 i = 0; i < 1000; ++i) {
					run();
				}

				if (ticks() - start > time) {
					break;
				}
			}
		}

		if (!has_written_board && writer) {
			auto p = win_prob();
			if (p > 0.97 || p < 0.03) {
#if 0
				u32 start = ticks();
				while (true) {
					for (u32 i = 0; i < 1000; ++i) {
						run();
					}

					if (ticks() - start > 5000) {
						break;
					}
				}

				if (p > 0.97 || p < 0.03) {
#else
				{
#endif
					char win = p > 0.5 ? 'y' : 'n';
					writer->write(&win, 1);
					u32 dest[32], opp[32];
					this->state.board.compact(state.board.headings[player_id].pos, state.board.headings[player_id ^ 1].pos, dest, opp);
					writer->write(dest, sizeof(dest));
					writer->write(opp, sizeof(opp));
					has_written_board = true;
				}
			}
		}

		return pick_next_move();
	}

	void advance(BoardMovesJ prev_moves, BoardMovesJ new_moves) {
		if (tree.root != 0) {
			auto& r = tree.get(tree.root);
			auto next_child = GN::child_index(prev_moves, new_moves);

			auto next_child_id = r.children[next_child];
			for (u32 i = 0; i < GN::ChildCount; ++i) {
				if (i != next_child && r.children[i]) {
					tree.free_subtree(r.children[i]);
				}
			}
			tree.free_node(tree.root);
			tree.root = next_child_id;
			--depth;
		}
	}

	virtual void update(Board board) {
		if (updates == 0) {
			Bot::update(board);
		} else {
			auto prev_moves = BoardMovesJ(state.board.headings[0].prev_move, state.board.headings[1].prev_move);
			Bot::update(board);
			advance(prev_moves, BoardMovesJ(state.board.headings[0].prev_move, state.board.headings[1].prev_move));
		}

		this->cells_left = 0;
		for (u32 y = 0; y < 16; ++y) {
			this->cells_left += popcount16(~board.get_row(y) & 0xffff);
		}
		++updates;
		visited.clear();
	}

	BoardMoves pick_next_move() {
		auto& r = tree.get(tree.root);

		u32 m = r.pick_best(rng, this->player_id);

		return move_from_child(this->state.board.headings[player_id].prev_move, m);
	}

	ChildIndexJ select(LcgPair& rng, NodeId node, u32 updates, typename GN::SelectInfo& select_info) {
		auto& n = tree.get(node);

		if (force_first_moves && updates < 3) {
			auto v = n.select<Variant>(rng, select_info);
			v.v[player_id] = 1;
			return v;
		}

		return n.select<Variant>(rng, select_info);
	}

	double win_prob() {
		if (tree.root != 0) {
			auto& r = tree.get(tree.root);
			auto m0 = r.players[this->player_id].move[0].ratio();
			auto m1 = r.players[this->player_id].move[1].ratio();
			auto m2 = r.players[this->player_id].move[2].ratio();
			return std::max(std::max(m0, m1), m2);
		} else {
			return 0.5;
		}
	}

	virtual void print_diag(double(&diag)[16][16]) {
		if (tree.root != 0) {
			auto& r = tree.get(tree.root);

			auto m0 = r.players[this->player_id].move[0].ratio();
			auto m1 = r.players[this->player_id].move[1].ratio();
			auto m2 = r.players[this->player_id].move[2].ratio();
			printf("[%2.2f] %2.2f %2.2f %2.2f %d %d (%2.2f MB)\n",
				std::max(std::max(m0, m1), m2),
				m0,
				m1,
				m2,
				sims,
				depth,
				tree.nodes_size() / (1024.0 * 1024.0));
			//printf("%f %d\n", (f64)this->tab_selections / this->selections, this->selections);
		}

		for (u32 y = 0; y < 16; ++y)
		for (u32 x = 0; x < 16; ++x) {
			diag[y][x] = this->amaf[this->player_id][y * 16 + x].ratio();
			this->amaf[this->player_id][y * 16 + x] = CounterF::empty();
		}
	}

	inline u32 expand_node(Board& board, MinMax<> life0, MinMax<> life1);
	void run();
};

template<int Variant, typename GN>
inline void warmup(McBotJ<Variant, GN>& self, GN& ch, Board& board, MinMax<> life0, MinMax<> life1) {
	for (u32 m = 0; m < 3; ++m) {
		auto next_pos0 = PlayerHeading::next_pos(board.headings[0].pos, board.headings[0].get_move(m));
		auto next_pos1 = PlayerHeading::next_pos(board.headings[1].pos, board.headings[1].get_move(m));
		bool wall_hit0 = board.is_wall(next_pos0);
		bool wall_hit1 = board.is_wall(next_pos1);

		ch.players[0].move[m] = wall_hit0 ? GN::Counter::unlikely() : GN::Counter::even_zero();
		ch.players[1].move[m] = wall_hit1 ? GN::Counter::unlikely() : GN::Counter::even_zero();

		if (Variant) {
			auto move_life0 = wall_hit0 ? MinMax<>::zero() : MinMax<>(std::max(life0.min, u8(1)), std::max(life0.max, u8(1)));
			auto move_life1 = wall_hit1 ? MinMax<>::zero() : MinMax<>(std::max(life1.min, u8(1)), std::max(life1.max, u8(1)));
			ch.players[0].life[m] = move_life0;
			ch.players[1].life[m] = move_life1;

			ch.life[m + 0 * 3][0] = move_life0;
			ch.life[m + 1 * 3][0] = move_life0;
			ch.life[m + 2 * 3][0] = move_life0;
			ch.life[0 + m * 3][1] = move_life1;
			ch.life[1 + m * 3][1] = move_life1;
			ch.life[2 + m * 3][1] = move_life1;
		}
	}
}

template<int Variant, typename GN>
inline u32 McBotJ<Variant, GN>::expand_node(Board& board, MinMax<> life0, MinMax<> life1) {
	u32 node = tree.alloc_node();
	auto& r = tree.get(node);
	warmup(*this, r, board, life0, life1);
	return node;
}

template<int Variant, typename GN>
void McBotJ<Variant, GN>::run() {
	GameState state(this->state);
	
	//u32 cur_depth = 0;
	u32 update_count = this->updates;
	u32 cells_left_count = this->cells_left;
	NodeId node = tree.root;
	MinMax<> result[2] = { MinMax<>(0, cells_left_count), MinMax<>(0, cells_left_count) };
	if (node == 0) {
		tree.root = node = expand_node(state.board, result[0], result[1]);
	}

	visited.clear();

	f32 p0;

	while (true) {
		typename GN::SelectInfo select_info;
		
		ChildIndexJ child = select(rng, node, update_count, select_info);
		++update_count;

		auto& n = tree.get(node);

		visited.push_back(VisitedEdge<GN>(node, child, state.board.headings[0].pos.compact(), state.board.headings[1].pos.compact(), select_info));
		u32 end_state = state.update(state.board.get_moves(child));
		++sims;

		if (end_state) {
			p0 = end_state == 1 ? 0.f : (end_state == 2 ? 1.f : 0.5f);
			if (Variant) {
				if (end_state & 1)
					result[0] = MinMax<>::zero();
				else {
					assert(result[0].max > 0);
					result[0].min = std::max(result[0].min, u8(1));
				}

				if (end_state & 2)
					result[1] = MinMax<>::zero();
				else {
					assert(result[1].max > 0);
					result[1].min = std::max(result[1].min, u8(1));
				}
			}
			break;
		}

		if (!n.has_stats(child)) {

			// TODO: Create "survivor"-subnode if there is no intersection after last move

			// Life is measured from the last node recorded in visited,
			// so we need to add 1 to the computed result.
			i32 score;
			if (Variant) {
				MinMax<> new_result[2];
				score_board4(state.board, false, this->rng, new_result);
				result[0] = result[0] & (new_result[0] + 1);
				result[1] = result[1] & (new_result[1] + 1);
				score = (i32)new_result[0].max - (i32)new_result[1].max;
			} else {
				u32 new_result[2];
				score_board_just_max(state.board, false, this->rng, new_result);
				score = (i32)new_result[0] - (i32)new_result[1];
			}
			p0 = score < 0 ? 0.f : (score > 0 ? 1.f : 0.5f);
			break;
		}

		auto child_index = child.index();
		auto child_id = n.children[child_index];

		if (Variant) {

			// TODO: Will this ever improve bounds?
			//result[0] = result[0] & n.players[0].life[child.v[0]];
			//result[1] = result[1] & n.players[1].life[child.v[1]];
			result[0] = result[0] & n.life[child_index][0];
			result[1] = result[1] & n.life[child_index][1];

			if (result[0].min > 0) {
				--result[0].min;
			}
			result[0].max -= 1;
			assert(result[0].is_valid());

			if (result[1].min > 0) {
				--result[1].min;
			}
			result[1].max -= 1;
			assert(result[1].is_valid());
		}

		if (child_id == 0) {
			tree.get(node).children[child_index] = child_id = expand_node(state.board, result[0], result[1]);

			// TEMP
			auto& ch = tree.get(child_id);
			/*
			if (ch.players[1].life[0].max == 0 && ch.players[1].life[1].max == 0 && ch.players[1].life[2].max == 0) {
				printf("x");
			}
			*/
		}

		node = child_id;
	}

	this->depth = std::max(this->depth, (u32)visited.size());

	f32 scores[2] = { p0, 1.f - p0 };
	bool prop[2] = { true, true };

	for (u32 i = (u32)visited.size(); i-- > 0; ) {
		auto& v = visited[i];

		auto& n = tree.get(v.node);

		n.players[0].move[v.edge.v[0]].add<Variant>(scores[0], 0, v);
		n.players[1].move[v.edge.v[1]].add<Variant>(scores[1], 1, v);
		amaf[0][v.to_pos[0]].add<Variant>(scores[0], 0, CounterF::SelectInfo());
		amaf[1][v.to_pos[1]].add<Variant>(scores[1], 1, CounterF::SelectInfo());

		if (Variant) {
			auto& l = n.life[v.edge.index()];
		
			for (u32 p = 0; p < 2; ++p) {
				if (prop[p]) {
					bool new_prop = l[p].and_(result[p]);
					prop[p] = new_prop;

					if (new_prop) {
						auto m = v.edge.v[p];
						// Player p has no control over which opponent move is chosen, so we use | instead of ^
						if (p == 0) {
							n.players[p].life[m] = n.life[0 + m * 3][p] | n.life[1 + m * 3][p] | n.life[2 + m * 3][p];
						} else {
							n.players[p].life[m] = n.life[m + 0 * 3][p] | n.life[m + 1 * 3][p] | n.life[m + 2 * 3][p];
						}
					}
				}
			}

			/*
			Prune options for a move:
			  A move is always worse than all unpruned moves of the opponent:
				p.life[m] < life[m, 0..<3] filter { # would_make } fold(or)

			  A move is always worse than all the other possible unpruned moves:
				p.life[m] < p.life[0..<3 filter {# != m}] fold(or)
			*/

			MinMax<> combined[2];

			if (prop[0] || prop[1]) {
				for (u32 p = 0; p < 2; ++p) {
					combined[p] = n.players[p].life[0] ^ n.players[p].life[1] ^ n.players[p].life[2];
				}
			}

			for (u32 p = 0; p < 2; ++p) {
				if (prop[p]) {
					for (u32 m = 0; m < 3; ++m) {
						if (!n.players[p].move[m].is_dead()
						  && (n.players[p].life[m] < combined[0] || n.players[p].life[m] < combined[1])) {
							// Prune because this move is worse than all other moves we can make
							// or all moves the opponent can make.
							n.players[p].move[m] = GN::Counter::unlikely();
							//printf("Pruned %d %d\n", p, m);
							++this->pruned;

							/* TODO: Free pruned subtrees
							if (p == 0) {
								for (u32 m2 = 0; m2 < 3; ++m2) {
									tree.free_subtree();
								}
							}
							*/
						}
					}
			
					result[p] = combined[p] + 1;
				}
			}
		}
	}

	this->total_updates += visited.size();
}

#endif // MCBOTJ_HPP
