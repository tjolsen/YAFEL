//
// Created by tyler on 9/15/17.
//

#ifndef YAFEL_TASK_HPP
#define YAFEL_TASK_HPP

#include "yafel_globals.hpp"
#include "utils/SmallVector.hpp"
#include <future>
#include <memory>
#include <type_traits>
#include <iostream>

YAFEL_NAMESPACE_OPEN

class TaskScheduler;

template<typename TASKSCHEDULER = TaskScheduler>
class Task
{
public:
    Task(TASKSCHEDULER &ts) : TS(ts) {}

    template<typename F, typename ...Args>
    inline auto addChild(F &&f, Args &&...args)
    {
        auto [task,fut] = TS.createTask(std::forward<F>(f), std::forward<Args>(args)...);

        children.push_back(task);

        return std::make_pair(std::move(task),std::move(fut));
    };

    inline void operator()() {
        (*internal_task)();
    }

private:
    // TaskScheduler is the one that actually handles task creation.
    // The create method is private and will be called by the scheduler
    // via its friend access.
    friend class TaskScheduler;



    template<typename F, typename ...Args>
    auto create(F &&f, Args &&...args) {
        using ret_type = std::result_of_t<F(Args...)>;

        std::promise<ret_type> task_promise;
        auto fut = task_promise.get_future();

        auto task_func =
                [this, task_promise=std::move(task_promise),
                        f = std::forward<F>(f),
                        args = std::forward_as_tuple(args...)]() mutable {

                    auto ret = std::apply(f,args);
                    task_promise.set_value(std::move(ret));

                    //Enqueue any children
                    for(auto &child : this->children) {
                        std::cout << "enqueueing child" << "\n";
                        this->TS.enqueue(child);
                    }
                };


        internal_task = std::make_unique<TaskImpl<decltype(task_func)> >(std::move(task_func));





        return fut;
    }




private:
    static constexpr int ReservedChildren = 2;
    SmallVector<std::shared_ptr<Task>, ReservedChildren> children;
    TASKSCHEDULER& TS;




    //----------------------------------------------------
    // Internal structures for a type-erased callable type that acts like a void() type.
    // Replacing the std::function, since I can't quite get into its guts the way I'd like to.
    struct TaskBase
    {
        virtual void operator()() = 0;
    };

    template<typename F>
    struct TaskImpl : public TaskBase
    {

        TaskImpl(F &&f) : callable(std::forward<F>(f)) {}

        void operator()() {
            callable();
        }

        F callable;
    };


    std::unique_ptr<TaskBase> internal_task;
};


YAFEL_NAMESPACE_CLOSE


#endif //YAFEL_TASK_HPP
