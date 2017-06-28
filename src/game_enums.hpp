#ifndef GAMEENUMS_h
#define GAMEENUMS_h

#include <cstdint>

typedef uint8_t u8;
typedef uint64_t u64;
typedef uint32_t u32;
typedef int32_t i32;
typedef uint16_t u16;
typedef double f64;
typedef float f32;
typedef size_t usize;

struct Vec {
	Vec() : Vec(0, 0) {
	}

	Vec(u32 x, u32 y) {
		v[0] = x;
		v[1] = y;
	}

	bool operator==(Vec const& other) const {
		return this->v[0] == other.v[0] && this->v[1] == other.v[1];
	}

	u8 compact() const {
		return u8((v[0] + (v[1] << 4)) & 0xff);
	}

	u32 x() const { return v[0]; }
	u32 y() const { return v[1]; }

	u32 v[2];
};

enum Player {
	Pl1 = 0, Pl2 = 1
};

enum BoardMoves {
	NONE = 0,
	RIGHT = 1 << 0,
	DOWN = 1 << 1,
	LEFT = 1 << 2,
	UP = 1 << 3,
};

struct BoardMovesJ {
	BoardMoves p[2];

	BoardMovesJ(BoardMoves p0, BoardMoves p1) {
		p[0] = p0;
		p[1] = p1;
	}
};

inline BoardMoves opposite_move(BoardMoves move) {
	switch (move) {
	case UP: return DOWN;
	case LEFT: return RIGHT;
	case DOWN: return UP;
	case RIGHT: return LEFT;
	}

	return NONE;
}

inline BoardMoves operator|(BoardMoves a, BoardMoves b) {
	return (BoardMoves)((u32)a | (u32)b);
}

inline BoardMoves move_from_child(BoardMoves prev_move, u32 child) {

	u32 prev = prev_move;
	u32 rot = (child - 1) & 3;
	return BoardMoves(((prev >> rot) | (prev << (4 - rot))) & 15);
}

inline u32 child_from_move(BoardMoves prev_move, BoardMoves new_move) {
	u32 prev = prev_move;
	if (BoardMoves(((prev << 1) | (prev >> (4 - 1))) & 15) == new_move) {
		return 0;
	} else if (BoardMoves(((prev >> 1) | (prev << (4 - 1))) & 15) == new_move) {
		return 2;
	} else {
		return 1;
	}
}

#endif
