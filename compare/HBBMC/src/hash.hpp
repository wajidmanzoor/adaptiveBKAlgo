#ifndef CUCKOO_HASH_H
#define CUCKOO_HASH_H

#include <cstdint>
#include <immintrin.h>
#include <cstdint>
#include <cstring>
#include <cstdio>
#include <utility>
#include <cassert>

#define EMPTY 0xffffffff

class CuckooHash {

private:
    int32_t capacity;
    int32_t mask;
    int32_t size;
    const int32_t buff_size = sizeof(int32_t);
    int32_t *hashtable;

    void rehash(int32_t **_table) {
        int32_t old_capacity = capacity;
        mask = mask == 0 ? 1 : ((mask << 1) | 1);
        capacity = (mask + 1) * buff_size;
        int32_t *new_hash = new int32_t[capacity];
        memset((new_hash), EMPTY, sizeof(int32_t) * capacity);
        for (int32_t i = 0; i < old_capacity; ++i) {
            if ((*_table)[i] != EMPTY) insert((*_table)[i], &new_hash);
        }
        std::swap((*_table), new_hash);
        delete[] new_hash;
    }

    void insert(const int32_t &_u, int32_t **_table) {

        int32_t hs = hash1(_u);
        for (int32_t i = 0; i < buff_size; ++i) {
            if ((*_table)[hs * buff_size + i] == EMPTY) {
                (*_table)[hs * buff_size + i] = _u;
                return;
            }
        }

        hs = hash2(_u);
        for (int32_t i = 0; i < buff_size; ++i) {
            if ((*_table)[hs * buff_size + i] == EMPTY) {
                (*_table)[hs * buff_size + i] = _u;
                return;
            }
        }

        bool use_hash1 = true;
        int32_t u = _u;
        for (int32_t i = 0; i < mask; ++i) {
            int32_t replaced;
            if (use_hash1) hs = hash1(u);
            else hs = hash2(u);
            int32_t j = 0;
            for (; j < buff_size; ++j) {
                if ((*_table)[hs * buff_size + j] == EMPTY) break;
            }
            if (buff_size == j) {
                replaced = std::move((*_table)[hs * buff_size]);
                j = 1;
                for (; j < buff_size; j++) {
                    (*_table)[hs * buff_size + j - 1] =
                            std::move((*_table)[hs * buff_size + j]);
                }
                (*_table)[hs * buff_size + j - 1] = u;
            } else {
                replaced = std::move((*_table)[hs * buff_size + j]);
                (*_table)[hs * buff_size + j] = u;
            }
            use_hash1 = hs == hash2(replaced);
            u = std::move(replaced);
            if (u == EMPTY) return;
        }
        rehash(_table);
        insert(u, _table);
    }

    int32_t hash1(const int32_t &x) { return x & mask; }

    int32_t hash2(const int32_t &x) { return ~x & mask; }

public:
    CuckooHash() {
        capacity = 0;
        hashtable = NULL;
        mask = 0;
        size = 0;
    }

    ~CuckooHash() {
        if (hashtable) delete[] hashtable;
    }

    void reserve(int32_t _size) {
        if (capacity >= _size) return;
        mask = mask == 0 ? 1 : ((mask << 1) | 1);
        while (_size >= mask * buff_size) mask = (mask << 1) | 1;
        capacity = (mask + 1) * buff_size;
        if (hashtable) delete[] hashtable;
        hashtable = new int32_t[capacity];
        memset(hashtable, EMPTY, sizeof(int32_t) * capacity);
    }

    void insert(const int32_t &_u) {
        if (find(_u)) return;
        insert(_u, &hashtable);
        size++;
    }

    bool find(const int32_t &_u) {
        int32_t hs1 = hash1(_u);
        int32_t hs2 = hash2(_u);

        __m128i cmp = _mm_set1_epi32(_u);
        __m128i b1 = _mm_load_si128((__m128i * ) & hashtable[buff_size * hs1]);
        __m128i b2 = _mm_load_si128((__m128i * ) & hashtable[buff_size * hs2]);
        __m128i flag = _mm_or_si128(_mm_cmpeq_epi32(cmp, b1), _mm_cmpeq_epi32(cmp, b2));

        return _mm_movemask_epi8(flag) != 0;
    }

    int32_t get_capacity() {
        return capacity;
    }

    int32_t get_size() {
        return size;
    }

    int32_t get_mask() {
        return mask;
    }

    int32_t *get_hashtable() {
        return hashtable;
    }
};

#endif