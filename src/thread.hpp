#ifndef THREAD_HPP
#define THREAD_HPP 1

#include <utility>
using std::move;

#ifdef _WIN32
#include <windows.h>

#define TL_NON_COPYABLE(name) \
	name(name const&) = delete; \
	name& operator=(name const&) = delete
#define TL_MOVABLE(name) \
	name(name&&) = default; \
	name& operator=(name&&) = default
#define TL_NON_MOVABLE(name) \
	name(name&&) = delete; \
	name& operator=(name&&) = delete

struct ThreadRef {
	HANDLE h;
	DWORD thread_id;

	ThreadRef()
		: h(INVALID_HANDLE_VALUE) {
	}

	ThreadRef(HANDLE h, DWORD thread_id)
		: h(h), thread_id(thread_id) {
	}

	void suspend() {
		SuspendThread(this->h);
	}

	void resume() {
		ResumeThread(this->h);
	}

	bool join() {
		return WAIT_OBJECT_0 == WaitForSingleObject(this->h, INFINITE);
	}

	void close() { CloseHandle(this->h); }

	static ThreadRef create(bool suspended, void(*f)(void* param), void* param) {
		DWORD thread_id;
		auto handle = CreateThread(NULL, 0, (LPTHREAD_START_ROUTINE)f, param, suspended ? CREATE_SUSPENDED : 0, &thread_id);

		return ThreadRef(handle, thread_id);
	}
};

struct Thread : ThreadRef {
	TL_NON_COPYABLE(Thread);
	//TL_MOVABLE(Thread);

	Thread(Thread&& other)
		: ThreadRef(other.h, other.thread_id) {
		other.h = INVALID_HANDLE_VALUE;
	}

	Thread& operator=(Thread&& other) = delete;

	explicit Thread(ThreadRef ref)
		: ThreadRef(ref) {
	}

	static Thread create(bool suspended, void(*f)(void* param), void* param) {
		return Thread(ThreadRef::create(suspended, f, param));
	}

	~Thread() {
		this->close();
	}
};

struct ThreadScoped : Thread {
	using Thread::Thread;

	ThreadScoped(Thread&& th)
		: Thread(move(th)) {
	}

	TL_NON_COPYABLE(ThreadScoped);
	TL_MOVABLE(ThreadScoped);

	template<typename F>
	static ThreadScoped create(bool suspended, F&& f) {
		F* p = new F(move(f));
		return Thread::create(suspended, [](void* p) {
			(*(F *)p)();
			delete p;
		}, p);
	}

	~ThreadScoped() {
		if (this->h != INVALID_HANDLE_VALUE) {
			this->join();
		}
	}
};

struct Mutex {
	struct Lock {
		TL_NON_COPYABLE(Lock);
		TL_MOVABLE(Lock);

		HANDLE mutex;

		Lock(HANDLE mutex) : mutex(mutex) {
		}

		~Lock() {
			ReleaseMutex(mutex);
		}
	};

	HANDLE mutex;
	
	Mutex()
		: mutex(CreateMutex(NULL, FALSE, NULL)) {
	}

	Lock lock() {
		WaitForSingleObject(mutex, INFINITE);
		return Lock(mutex);
	}
};
#endif

#endif
