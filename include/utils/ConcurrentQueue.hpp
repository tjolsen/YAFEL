#ifndef _YAFEL_CONCURRENTQUEUE_HPP
#define _YAFEL_CONCURRENTQUEUE_HPP

#include <list>
#include <mutex>

template<typename T>
class ConcurrentQueue {

public:
    
    void enqueue(T x) {
	list<T> tmp;
	tmp.push_back(move(x));
	{
	    std::lock_guard<std::mutex>> lock(mtx);
	    q.splice(end(q), tmp);
	}
    }

    T dequeue() {
	std::lock_guard<std::mutex> lock(mtx);
	T retval = q.front();
	q.pop_front();
	return retval;
    }

private:
    std::mutex mtx;
    std::list<T> q;
    
};

#endif
