#ifndef YAFEL_TASK_SCHEDULER
#define YAFEL_TASK_SCHEDULER

/**
 * \file
 */

#include "yafel_globals.hpp"

#include <thread>
#include <mutex>
#include <deque>
#include <future>
#include <functional>
#include <type_traits>
#include <atomic>
#include <vector>
#include <tuple>

//Platform-specific
#include <pthread.h>
#include <sched.h>

YAFEL_NAMESPACE_OPEN


namespace worker_global {

thread_local int worker_id{-1};

}

/**
 * \class TaskScheduler
 *
 * Task scheduling system written in pure C++14 with no external
 * library dependencies. The class is modeled on one presented in
 * a talk given by Sean Parent called "Better Code: Concurrency".
 * Search for talk on Youtube.
 * Implements a "task stealing" scheduling system.
 * Performance over the more naive "ThreadPool" is gained by spinning
 * over N task queues to avoid contention, both on the enqueue and 
 * dequeue side.
 * Combines ideas from the Sean Parent implementation and the ThreadPool
 * implementation (see ThreadPool in this project) from github by returning
 * a std::future<T()> from the TaskScheduler::enqueue() method.
 */
class TaskScheduler
{
public:
    using lock_t = std::unique_lock<std::mutex>;

    constexpr static int spin_count{32};
    const int count;

    std::vector<std::thread> workers;
    std::vector<std::deque<std::function<void()>>> tq;
    std::vector<std::mutex> mtxs;
    std::vector<std::condition_variable> cvs;
    std::vector<bool> done;
    std::atomic<int> _index;

    inline TaskScheduler(int nthreads = std::thread::hardware_concurrency())
            : count(nthreads),
              workers(nthreads),
              tq(nthreads),
              mtxs(nthreads),
              cvs(nthreads),
              done(nthreads, false),
              _index(0)
    {
        cpu_set_t cpus;

        for (int i = 0; i < nthreads; ++i) {
            workers[i] = std::thread([i, this]() {
                // Set thread-local "global" worker id.
                // Supposed to be accessible by functions run by workers.
                // Useful for "reduction" operations that accumulate into
                // multiple buffers (eg: summing vector entries into a vector of partial sums,
                // which is subsequently reduced later).
                worker_global::worker_id = i;
                run(i);
            });

            //Set thread affinity. Modulo just in case some idiot decides to spawn more threads than cores.
            // Wraparound to zero seems the most sensible way to handle this without yelling at the user
            // too much.
            CPU_ZERO(&cpus);
            CPU_SET(i % std::thread::hardware_concurrency(),&cpus);
            auto ret = pthread_setaffinity_np(workers[i].native_handle(), sizeof(cpu_set_t), &cpus);
            if(ret != 0) {
                throw std::runtime_error("Failed to set thread affinity");
            }
        }
    }

    inline ~TaskScheduler()
    {
        for (int i = 0; i < count; ++i) {
            {
                std::unique_lock<std::mutex> lock{mtxs[i]};
                done[i] = true;
            }
            cvs[i].notify_all();
        }

        for (auto &w : workers)
            w.join();
    }

    template<typename F, typename ...Args>
    auto enqueue(F &&f, Args &&...args)
    -> std::future<typename std::result_of<F(Args...)>::type>
    {

        using ret_type = typename std::result_of<F(Args...)>::type;

        auto task = std::make_shared<std::packaged_task<ret_type()> >
                (
                        [f=std::forward<F>(f),
                         args = std::forward_as_tuple(args...)/*std::make_tuple(std::forward<Args>(args)...)*/]() {
                            return std::apply(f,args);
                        }
                );
                //(std::bind(std::forward<F>(f), std::forward<Args>(args)...));

        auto ret = task->get_future();
        auto i = _index++;

        //try to push without blocking
        for (int n = 0; n != count * spin_count; ++n) {
            auto idx = (i + n) % count;
            {
                lock_t lock{mtxs[idx], std::try_to_lock};
                if (!lock) continue;

                tq[idx].emplace_back([task]() { (*task)(); });
            }
            cvs[idx].notify_one();
            return ret;
        }
        {
            int idx = i%count;
            lock_t lock{mtxs[idx]};
            tq[idx].emplace_back([task]() { (*task)(); });
        }
        cvs[i%count].notify_one();
        return ret;
    }

private:

    void run(unsigned i)
    {

        while (true) {
            std::function<void()> f;
            bool try_failure = true;
            for (int n = 0; n != count * spin_count; ++n) {
                auto idx = (i + n) % count;
                lock_t lock{mtxs[idx], std::try_to_lock};

                if (!lock || tq[idx].empty()) continue;

                f = std::move(tq[idx].front());
                tq[idx].pop_front();
                try_failure = false;
                break;
            }
            if (try_failure) {
                lock_t lock{mtxs[i]};
                while (tq[i].empty() && !done[i])
                    cvs[i].wait(lock);
                if (tq[i].empty())
                    break;
                f = std::move(tq[i].front());
                tq[i].pop_front();
            }
            if (!f) break;

            f();
        } //end while

    }//end run


};

YAFEL_NAMESPACE_CLOSE

#endif
