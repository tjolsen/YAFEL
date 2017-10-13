//
// Created by tyler on 7/16/17.
//

#ifndef YAFEL_SMALLVECTOR_HPP
#define YAFEL_SMALLVECTOR_HPP

#include "yafel_globals.hpp"
#include <initializer_list>
#include <algorithm>
#include <type_traits>
#include <stdexcept>

YAFEL_NAMESPACE_OPEN

/**
 * \class SmallVectorImpl
 * \brief Implementation of an STL-like vector class for use with SmallVector
 * @tparam T
 */
template<typename T>
class SmallVectorImpl
{
public:
    using iterator = T *;
    using const_iterator = T const *;
    using reference = T &;
    using const_reference = T const &;



protected:
    T *p_begin, *p_end;
    std::size_t _capacity{0};
    bool _isSmall{false};
    using storage_type = std::aligned_storage_t<sizeof(T), alignof(T)>;


    SmallVectorImpl() : p_begin(nullptr), p_end(nullptr) {}

    SmallVectorImpl(T *begin, T *end, std::size_t capacity, bool small)
            : p_begin(begin), p_end(end), _capacity(capacity), _isSmall(small) {}


    void freeBuffer()
    {
        if (p_begin) {
            delete[] reinterpret_cast<storage_type *>(p_begin);
        }
    }

    inline static void destroyRange(iterator it_start, iterator it_end)
    {
        if constexpr (!std::is_trivially_destructible<T>::value)
        {
            while (it_end != it_start) {
                --it_end;
                it_end->~T();
            }
        }
    }

public:
    bool isSmall() const { return _isSmall; }

    T &operator[](std::size_t idx) { return p_begin[idx]; }

    T const &operator[](std::size_t idx) const { return p_begin[idx]; }

    inline iterator begin() { return p_begin; }

    inline iterator end() { return p_end; }

    inline const_iterator cbegin() const { return p_begin; }

    inline const_iterator cend() const { return p_end; }

    inline iterator &data() { return p_begin; }

    inline iterator &dataEnd() { return p_end; }

    void reserve(std::size_t newCapacity)
    {
        if (newCapacity > _capacity) {

            auto ptr = reinterpret_cast<T *>(new storage_type[newCapacity]);
            if (ptr) {
                std::move(p_begin, p_end, ptr);
                auto oldSize = size();
                auto oldBuffer = p_begin;
                p_begin = ptr;
                p_end = ptr + oldSize;
                _capacity = newCapacity;

                if (isSmall()) {
                    _isSmall = false;
                } else {
                    delete[] reinterpret_cast<storage_type *>(oldBuffer);
                }
            }
        }
    }

    void resize(std::size_t newSize)
    {
        auto old_size = this->size();
        reserve(newSize);
        p_end = p_begin + newSize;
        if (newSize < old_size) {
            destroyRange(p_end, p_begin + old_size);
        }
    }

    void resize(std::size_t newSize, const T &fillElement)
    {
        reserve(newSize);
        auto new_p_end = p_begin + newSize;
        for (auto ptr = p_end; ptr < new_p_end; ++ptr) {
            *ptr = fillElement;
        }
        p_end = new_p_end;
    }

    inline reference back() { return *(p_end - 1); }

    inline const_reference back() const { return *(p_end - 1); }

    inline reference front() { return *p_begin; }

    inline const_reference front() const { return *p_begin; }

    void push_back(const T &element)
    {
        if (size() >= capacity()) {
            reserve(2 * _capacity);
        }
        new(p_end) T(element);
        ++p_end;
    }

    template<typename ...Args>
    void emplace_back(Args ...args)
    {
        if (size() >= capacity()) {
            reserve(2 * _capacity);
        }
        new(p_end) T(args...);
        ++p_end;
    }

    void pop_back()
    {
        if (size() > 0) {
            if constexpr (!std::is_trivially_destructible<T>::value)
            {
                --p_end->~T();
            } else {
                --p_end;
            }
        }
        throw std::runtime_error("Bad pop on empty SmallVector");
    }

    std::size_t size() const { return p_end - p_begin; }

    std::size_t &capacity() { return _capacity; }

    std::size_t const &capacity() const { return _capacity; }
};


/**
 * \class SmallVector
 * \brief Small-size-optimized vector class inspired by the LLVM SmallVector
 * (but without the need to extricate that beast from their source tree).
 * 
 * The SmallVector (together with SmallVectorImpl) provide a relatively
 * feature-complete container class comparable in interface to the std::vector.
 * The only difference is that it follows different iterator invalidation
 * rules, so a small-size optimization can be layered on top of an otherwise
 * efficient vector container.
 *
 * @tparam T 
 * @tparam N
 */
template<typename T, std::size_t N>
class SmallVector : public SmallVectorImpl<T>
{
    using storage_type = typename SmallVectorImpl<T>::storage_type;
public:

    SmallVector()
            : SmallVectorImpl<T>(reinterpret_cast<T *>(&buffer[0]),
                                 reinterpret_cast<T *>(&buffer[0]), N, true) {}

    SmallVector(std::size_t n)
            : SmallVectorImpl<T>(reinterpret_cast<T *>(&buffer[0]),
                                 reinterpret_cast<T *>(&buffer[0]), N, true)
    {
        this->resize(n);
    }

    SmallVector(std::size_t n, const T &val) : SmallVector(n)
    {
        for (auto it = this->begin(); it < this->end(); ++it) {
            new(it) T(val);
        }
    }

    SmallVector(std::initializer_list<T> L) : SmallVector()
    {
        this->reserve(L.size());
        for (auto const &l : L) {
            this->push_back(l);
        }
    }

    SmallVector(const SmallVector<T, N> &rhs) :
            SmallVectorImpl<T>(reinterpret_cast<T *>(&buffer[0]),
                               reinterpret_cast<T *>(&buffer[0]), N, true)
    {
        this->resize(rhs.size());
        std::copy(rhs.cbegin(), rhs.cend(), this->begin());
    }

    template<std::size_t N2>
    SmallVector(const SmallVector<T, N2> &rhs) :
            SmallVectorImpl<T>(reinterpret_cast<T *>(&buffer[0]),
                               reinterpret_cast<T *>(&buffer[0]), N, true)
    {
        this->resize(rhs.size());
        std::copy(rhs.cbegin(), rhs.cend(), this->begin());
    }

    SmallVector(SmallVector<T, N> &&rhs)
    {
        if (rhs.isSmall()) {
            this->p_begin = reinterpret_cast<T *>(&buffer[0]);
            std::copy(rhs.cbegin(), rhs.cend(), this->begin());
            this->p_end = this->p_begin + rhs.size();
            this->_capacity = N;
        } else {
            this->p_begin = rhs.begin();
            this->p_end = rhs.end();
            this->_capacity = rhs.capacity();
            this->_isSmall = false;

            rhs.data() = nullptr;
            rhs.dataEnd() = nullptr;
            rhs.capacity() = 0;
        }
    }

    template<std::size_t N2>
    SmallVector(SmallVector<T, N2> &&rhs)
    {
        if (rhs.isSmall()) {
            //must copy from rhs
            if (N >= N2 || rhs.size() <= N) {
                //rhs fits into small-size buffer
                this->p_begin = reinterpret_cast<T *>(&buffer[0]);
                this->_capacity = N;
                this->_isSmall = true;
                std::move(rhs.cbegin(), rhs.cend(), this->begin());
            } else {
                this->resize(rhs.size());
                std::move(rhs.cbegin(), rhs.cend(), this->begin());
            }
        } else {
            //move from rhs vector
            this->p_begin = rhs.begin();
            this->p_end = rhs.end();
            this->_capacity = rhs.capacity();
            this->_isSmall = false;

            rhs.data() = nullptr;
            rhs.dataEnd() = nullptr;
            rhs.capacity() = 0;
        }
    }

    //Copy assignment
    template<std::size_t N2>
    SmallVector<T, N> &operator=(const SmallVector<T, N2> &rhs)
    {
        this->copy(rhs);
        return *this;
    }

    SmallVector<T, N> &operator=(const SmallVector<T, N> &rhs)
    {
        this->copy(rhs);
        return *this;
    }

    //Move Assignment
    template<std::size_t N2>
    SmallVector<T, N> &operator=(SmallVector<T, N2> &&rhs)
    {
        this->move(std::forward<SmallVector<T, N2>>(rhs));
        return *this;
    }

    ~SmallVector()
    {
        this->destroyRange(this->begin(), this->end());
        if (&((*this)[0]) != reinterpret_cast<T *>(&buffer[0])) {
            this->freeBuffer();
        }
    }

protected:

    template<std::size_t N2>
    void copy(const SmallVector<T, N2> &rhs)
    {
        if (this->capacity() >= rhs.size()) {
            //current storage (small or not) has space for rhs
            std::copy(rhs.cbegin(), rhs.cend(), this->begin());
            this->dataEnd() = this->data() + rhs.size();
        } else {
            this->resize(rhs.size());
            std::copy(rhs.cbegin(), rhs.cend(), this->begin());
        }
    }

    template<std::size_t N2>
    void move(SmallVector<T, N2> &&rhs)
    {
        if (rhs.isSmall()) {
            copy(rhs);
        } else {
            //clean up current buffer and free if (!this->isSmall())
            this->destroyRange(this->begin(), this->end());
            if (!this->isSmall()) {
                this->freeBuffer();
            }

            //move from heap-allocated rhs buffer
            this->_isSmall = false;
            this->data() = rhs.data();
            this->dataEnd() = rhs.dataEnd();
            this->capacity() = rhs.capacity();

            rhs.data() = nullptr;
            rhs.dataEnd() = nullptr;
            rhs.capacity() = 0;
        }
    }

    storage_type buffer[N];
};


YAFEL_NAMESPACE_CLOSE


#endif //PARTICLEDYN_SMALLVECTOR_HPP
