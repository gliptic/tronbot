#ifndef GAMETREE_HPP
#define GAMETREE_HPP 1

#include "game_enums.hpp"
#include "board.hpp"
#include <vector>
#include <algorithm>

typedef u32 NodeId;
using std::vector;

static u32 const iterations = 20000; // 00;
static bool const deterministic = false;

struct GameNode {
	u32 wins, total;
	NodeId children[3];

	static u32 const ChildCount = 3;

	GameNode()
		: wins(0), total(0) {

		children[0] = 0;
		children[1] = 0;
		children[2] = 0;
	}
};

template<typename NodeT>
struct BasicGameTree1 {
	vector<NodeT> nodes;
	NodeId root;
	NodeId first_free;

	BasicGameTree1() : first_free(0), root(0) {
		// Dummy node
		nodes.push_back(NodeT());
	}

	NodeId alloc_node() {
		NodeId index = this->first_free;
		if (index == 0) {
			index = NodeId(nodes.size());
			nodes.push_back(NodeT());
		} else {
			this->first_free = nodes[index].children[0];
			nodes[index] = NodeT();
		}

		return index;
	}

	NodeT& get(NodeId index) {
		return nodes[index];
	}

	usize nodes_size() const {
		return nodes.size() * sizeof(NodeT);
	}

	void free_node(NodeId index) {
		nodes[index].children[0] = first_free;
		first_free = index;
	}

	void free_subtree(NodeId index) {
		auto& n = get(index);
		for (u32 i = 0; i < NodeT::ChildCount; ++i) {
			if (n.children[i]) free_subtree(n.children[i]);
		}
		free_node(index);
	}
};

template<typename NodeT>
struct BasicGameTree {
	u8* nodes;
	u32 nodes_next;
	u32 nodes_cap;
	NodeId root;
	NodeId first_free;

	BasicGameTree() : first_free(0), root(0) {
		// Dummy node
		//nodes.push_back(NodeT());
		u32 cap = 1024 * sizeof(NodeT);
		nodes = (u8 *)malloc(cap);
		nodes_next = sizeof(NodeT);
		nodes_cap = cap;
	}

	~BasicGameTree() {
		free(nodes);
	}

	usize nodes_size() const {
		return nodes_next;
	}

	NodeId alloc_node() {
		NodeId index = this->first_free;
		if (index == 0) {
			index = nodes_next;
			get(index) = NodeT();
			nodes_next += sizeof(NodeT);
			if (nodes_next == nodes_cap) {
				nodes_cap *= 2;
				nodes = (u8 *)realloc(nodes, nodes_cap);
			}
		} else {
			NodeT& next = get(index);
			this->first_free = next.children[0];
			next = NodeT();
		}

		return index;
	}

	NodeT& get(NodeId index) {
		return *(NodeT *)(nodes + index);
	}

	void free_node(NodeId index) {
		get(index).children[0] = first_free;
		first_free = index;
	}

	void free_subtree(NodeId index) {
		auto& n = get(index);
		for (u32 i = 0; i < NodeT::ChildCount; ++i) {
			if (n.children[i]) free_subtree(n.children[i]);
		}
		free_node(index);
	}
};

#endif // GAMETREE_HPP
