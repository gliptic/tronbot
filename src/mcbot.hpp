#ifndef MCBOT_HPP
#define MCBOT_HPP 1

#include "bot.hpp"
#include "gametree.hpp"
#include "rand.hpp"

struct McBot : Bot {

	GameTree tree;
	int variant;
	u32 sims, depth;

	McBot(int variant, LcgPair& rng) : variant(variant), sims(0), depth(0), Bot(rng) {
	}

	virtual BoardMoves move(int time) {
		if (deterministic) {
			for (u32 i = 0; i < iterations; ++i) {
				run(this->state);
			}
		} else {
			u32 start = ticks();
			u32 utime = (u32)time * 95 / 100; // Margin of error
			while (true) {
				for (u32 i = 0; i < 1000; ++i) {
					run(this->state);
				}

				if (ticks() - start > time) {
					break;
				}
			}
		}

		return pick_next_move(rng);
	}

	void advance(BoardMoves prev_move, BoardMoves new_move) {
		if (tree.root != 0) {
			auto& r = tree.get(tree.root);
			auto next_child = child_from_move(prev_move, new_move);
			auto next_child_id = r.children[next_child];
			for (u32 i = 0; i < 3; ++i) {
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
		auto prev_move0 = state.prev_moves[player_id];
		auto prev_move1 = state.prev_moves[player_id ^ 1];
		Bot::update(board);
		advance(prev_move0, state.prev_moves[player_id]);
		advance(prev_move1, state.prev_moves[player_id ^ 1]);
	}

	f64 ucb(NodeId node, double c) {
		auto& n = tree.get(node);
		// TODO: Table of c * sqrt(log(n.total + 1) / n.total)
		return (double)n.wins / n.total + c * sqrt(log(n.total + 1) / n.total);
	}

	f64 ratio(NodeId node) {
		auto& n = tree.get(node);
		// TODO: Table of c * sqrt(log(n.total + 1) / n.total)
		return (double)n.wins / n.total;
	}

	BoardMoves pick_next_move(LcgPair& rng) {
		auto& r = tree.get(tree.root);

		if (variant == 0) {
			i32 best = -1, best_wins = -1;
			u32 best_children[3];
			u32 best_children_count = 0;

			for (u32 i = 0; i < 3; ++i) {
				auto child = r.children[i];

				auto& n = tree.get(child);
				i32 cand_best = n.total;
				if (cand_best > best || (cand_best == best && n.wins > best_wins)) {
					best = cand_best;
					best_wins = n.wins;
					best_children_count = 1;
					best_children[0] = i;
				} else if (cand_best == best && n.wins == best_wins) {
					best_children[best_children_count++] = i;
				}
			}

			return move_from_child(this->state.prev_moves[player_id], best_children[rng.get_u32(best_children_count)]);
		} else {
			f64 best = -1.0;
			u32 best_children[3];
			u32 best_children_count = 0;

			for (int i = 0; i < 3; ++i) {
				auto child = r.children[i];

				f64 cand_best = ratio(child);
				if (cand_best > best) {
					best = cand_best;
					best_children_count = 1;
					best_children[0] = i;
				} else if (cand_best == best) {
					best_children[best_children_count++] = i;
				}
			}

			return move_from_child(this->state.prev_moves[player_id], best_children[rng.get_u32(best_children_count)]);
		}
	}

	u32 select(LcgPair& rng, NodeId node) {
		auto& n = tree.get(node);

		int unvisited_count = (n.children[0] == 0) + (n.children[1] == 0) + (n.children[2] == 0);
		if (unvisited_count > 0) {
			int pick = rng.get_u32(unvisited_count);
			for (int i = 0; i < 3; ++i) {
				if (n.children[i] == 0) {
					if (pick == 0) {
						return i;
					}
					--pick;
				}
			}
			abort();
		} else {
			f64 best = -1.0;
			u32 best_children[3];
			u32 best_children_count = 0;

			for (int i = 0; i < 3; ++i) {
				auto child = n.children[i];

				f64 cand_best = ucb(child, sqrt(2.0));
				if (cand_best > best) {
					best = cand_best;
					best_children_count = 1;
					best_children[0] = i;
				} else if (cand_best == best) {
					best_children[best_children_count++] = i;
				}
			}

			return best_children[rng.get_u32(best_children_count)];
		}
	}

	/*
	virtual void print_diag() {
	printf("%d M\n", tree.nodes.size() * sizeof(GameNode));
	}
	*/

	virtual void print_diag() {
		//printf("%d M\n", tree.nodes.size() * sizeof(GameNodeJ));
		if (tree.root != 0) {
			auto& r = tree.get(tree.root);

			if (tree.get(r.children[0]).total == 0
				|| tree.get(r.children[1]).total == 0
				|| tree.get(r.children[2]).total == 0) {
				printf("??\n");
			}

			auto m0 = ratio(r.children[0]);
			auto m1 = ratio(r.children[1]);
			auto m2 = ratio(r.children[2]);
			printf("[%2.2f] %2.2f %2.2f %2.2f %d %d (%2.2f MB)\n",
				std::max(std::max(m0, m1), m2),
				m0,
				m1,
				m2,
				sims / 2,
				depth / 2,
				tree.nodes.size() * sizeof(GameNode) / (1024.0 * 1024.0));
		}
	}

	void run(GameState const& cur_state) {
		vector<NodeId> visited;

		NodeId node = tree.root;
		if (node == 0) {
			// Unvisited
			node = tree.alloc_node();
			tree.root = node;
		}

		visited.push_back(node);

		GameStateTurnBased state(this->player_id, cur_state);

		u32 end_state = 0;
		u32 cur_depth = 0;

		while (end_state == 0) {
			u32 child = select(rng, node);

			auto child_id = tree.get(node).children[child];
			if (child_id == 0) {
				// Unvisited
				child_id = tree.alloc_node();
				tree.get(node).children[child] = child_id;
				node = child_id;
				visited.push_back(node);
				auto next_move = state.get_move(child);
				end_state = state.do_move(next_move);
				++sims;
				++cur_depth;
				break;
			}

			node = child_id;
			visited.push_back(node);
			auto next_move = state.get_move(child);
			end_state = state.do_move(next_move);
			++sims;
			++cur_depth;
			//assert(end_state == 0);
		}

		this->depth = std::max(this->depth, cur_depth);

		// Play-out
		while (end_state == 0) {
			auto moves = state.board.legal_moves((Player)state.turn);
			BoardMoves next_move;

			if (moves == NONE) {
				next_move = state.prev_moves[state.turn];
			} else {
				while (true) {
					int move = rng.next() & 3;
					if (moves & (1 << move)) {
						next_move = BoardMoves(1 << move);
						break;
					}
				}
			}

			end_state = state.do_move(next_move);
			++sims;
		}

		for (u32 i = 0; i < visited.size(); ++i) {
			NodeId node = visited[i];
			u32 player = (i + (u32)this->player_id + 1) % 2;

			auto& n = tree.get(node);
			++n.total;
			if (((end_state >> player) & 1) == 0) {
				++n.wins;
			}
		}
	}
};

#endif // MCBOT_HPP
