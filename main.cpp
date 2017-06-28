#include <cstdint>
#include "src/parser.hpp"
#include "src/bot.hpp"

#include "src/silly_bot.hpp"
#include "src/mcbotj.hpp"

#include "driver.hpp"
#include "src/thread.hpp"
#include "gfx.hpp"

template<typename T>
inline void asserteq(T x, T y) {
	if (x != y) {
		printf("Fail\n"); abort();
	}
}

void test() {
	asserteq(move_from_child(DOWN, 0), LEFT);
	asserteq(move_from_child(LEFT, 0), UP);
	asserteq(move_from_child(UP, 0), RIGHT);
	asserteq(move_from_child(RIGHT, 0), DOWN);

	asserteq(move_from_child(DOWN, 1), DOWN);
	asserteq(move_from_child(LEFT, 1), LEFT);
	asserteq(move_from_child(UP, 1), UP);
	asserteq(move_from_child(RIGHT, 1), RIGHT);

	asserteq(move_from_child(DOWN, 2), RIGHT);
	asserteq(move_from_child(LEFT, 2), DOWN);
	asserteq(move_from_child(UP, 2), LEFT);
	asserteq(move_from_child(RIGHT, 2), UP);

	// Moves
	asserteq(child_from_move(DOWN, DOWN), (u32)1);
	asserteq(child_from_move(LEFT, LEFT), (u32)1);
	asserteq(child_from_move(UP, UP), (u32)1);
	asserteq(child_from_move(RIGHT, RIGHT), (u32)1);

	asserteq(child_from_move(DOWN, LEFT), (u32)0);
	asserteq(child_from_move(LEFT, UP), (u32)0);
	asserteq(child_from_move(UP, RIGHT), (u32)0);
	asserteq(child_from_move(RIGHT, DOWN), (u32)0);

	asserteq(child_from_move(DOWN, RIGHT), (u32)2);
	asserteq(child_from_move(LEFT, DOWN), (u32)2);
	asserteq(child_from_move(UP, LEFT), (u32)2);
	asserteq(child_from_move(RIGHT, UP), (u32)2);
}

#define SLEE 1
#define RENDER 1

template<int Variant, typename GN = GameNodeJ>
struct BotProcess : BaseBot {
	Thread thread;
	Mutex mutex;
	LcgPair rng;
	McBotJ<Variant, GN> bot;
	bool running;

	TL_NON_COPYABLE(BotProcess);
	TL_NON_MOVABLE(BotProcess);

	BotProcess(LcgPair& rng)
		: bot(this->rng, 0)
		, thread(Thread::create(true, run, this))
		, rng(rng)
		, running(false) {
	}

	~BotProcess() {
		{
			auto lock(mutex.lock());
			running = false;
		}
		thread.join();
	}

	static void run(void* p) {
		BotProcess& self = *(BotProcess *)p;
		for (;;) {
			auto lock(self.mutex.lock());

			if (!self.running)
				return;

			for (u32 i = 0; i < 1000; ++i) {
				self.bot.run();
			}
		}
	}

	virtual BoardMoves move(int time) {
		auto lock(mutex.lock());
		return bot.move(time);
	}

	// Update
	virtual void round(int time) {
	}

	virtual void update(Board board) {
		auto lock(mutex.lock());
		if (!running) {
			running = true;
			thread.resume();
		}
		return bot.update(board);
	}

	// Settings
	virtual void timebank(int time) {}
	virtual void time_per_move(int time) {}
	virtual void your_bot(std::string name) {}

	virtual void your_bot_id(Player player_id) {
		auto lock(mutex.lock());
		return bot.your_bot_id(player_id);
	}

	virtual void field_width(int width) {
		auto lock(mutex.lock());
		return bot.field_width(width);
	}

	virtual void player_names(std::string player1, std::string player2) {
	}

	virtual void field_height(int height) {
		auto lock(mutex.lock());
		return bot.field_height(height);
	}

	virtual void print_diag(double(&diag)[16][16]) {
		auto lock(mutex.lock());
		return bot.print_diag(diag);
	}
};

#pragma comment(lib, "d3d9.lib")

int test_skip() {
	f32 sqrt2 = sqrtf(2.f);
	bool in_seq = false;
	f32 prev = 0.f;
	for (f32 r = 0.f; r < 0.99f; r += 0.01f) {
		u32 count = 0;
		u32 max = 10000;
		for (u32 total = 1; total < max; ++total) {
			CounterF c(r, total);
			CounterF c2(c);
			c2.add<0>(1.f, 0, CounterF::SelectInfo());
			if (c.ucb(sqrt2, 0) <= c2.ucb(sqrt2, 0)) {
				++count;
			}
		}
		if (count >= max - 1) {
			if (!in_seq) {
				printf(">= %2.2f\n", r);
				in_seq = true;
			}
			prev = r;
		} else {
			if (in_seq) {
				printf("<= %2.2f\n", prev);
				in_seq = false;
			}
		}
	}

	return 0;
}

int bench() {
	// For debug:
	scr = tigrWindow(1024, 768, "mcbot", 0);
	LcgPair rng(2, 0);
	McBotJ<1, GameNodeJ> bot(rng, 0); // Exp3
	bot.field_width(16);
	bot.field_height(16);
	bot.your_bot_id(Pl1);
	GameState state(GameState::standard());
	bot.update(state.board);
	auto before = ticks();
	for (int i = 0; i < 400000; ++i) {
		if (i == 0) debug_board = true;
		bot.run();
	}
	auto after = ticks();
	printf("%d ms\n", after - before);
	return 0;
}


int tests() {
#if 0
	Bot bot;
	Parser parser(bot);
	parser.Parse();
#else
	u32 wins[3] = {0, 0, 0};

	LcgPair rng(2, 0);

	test();

	int f = 0;
#if SLEE
	int const skip = 1;
#else
	int const skip = 100;
	//u32 prev_ticks = ticks();
#endif

	Writer w("boards.dat");
	Writer* pw = 0;

#if RENDER
	scr = tigrWindow(1024, 768, "mcbot", 0);
	auto checkClosed = [&] { return tigrClosed(scr); };
#else
	auto checkClosed = [&] { return false; };
#endif

	for (int c = 0; c < 1000 && !checkClosed(); ++c) {

		//Game game(rng, std::make_unique<McBotJ<1>>(rng), std::make_unique<McBotJ<0>>(rng));
		//Game game(rng, std::make_unique<BotProcess>(1, rng), std::make_unique<BotProcess>(1, rng));
		Game game(rng, std::make_unique<BotProcess<1, GameNodeJ>>(rng), std::make_unique<BotProcess<0>>(rng));
		//Game game(rng, std::make_unique<BotProcess<1>>(rng), std::make_unique<SillyBot>(rng, 0));

		rng.next();
		
		u32 status;
		
		do {
#if SLEE
			//_sleep(50);
#endif
			auto moves = game.update_bots();

#if RENDER
			if ((f % skip) == 0) {
				double diag[16][16];
				game.bot[0]->print_diag(diag);
				BitBoard* boards[1] = { &game.board };
				TPixel colors[1] = { tigrRGB(100, 100, 100) };
				render_board(1, boards, colors, game.board.headings, &diag);
				game.bot[1]->print_diag(diag);

				printf("[%d %d %d]\n", wins[0], wins[1], wins[2]);
			}
#endif
			status = game.update_(moves);
			++f;
		} while (status == 0);

#if SLEE
		_sleep(2000);
#endif

		if (status == 3) ++wins[2];
		if (status == 2) ++wins[0];
		if (status == 1) ++wins[1];
	}

	system("cls");
	printf("[%d %d %d]\n", wins[0], wins[1], wins[2]);

#if RENDER
	tigrFree(scr);
#endif

#endif
	return 0;
}

int main() {
	//return test_skip();
	return bench();
	//return tests();
}

