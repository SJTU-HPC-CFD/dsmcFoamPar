/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2021 hyStrath
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of hyStrath, a derivative work of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::dsmcFastRng

Description
    Lightweight per-thread RNG for DSMC OpenMP hot paths.

\*---------------------------------------------------------------------------*/

#ifndef dsmcFastRng_H
#define dsmcFastRng_H

#include "label.H"
#include "scalar.H"
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <new>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class dsmcFastRng
{
    uint64_t s0_;
    uint64_t s1_;

    static inline uint64_t splitMix64(uint64_t z)
    {
        uint64_t v = z;
        v = (v ^ (v >> 30))*0xbf58476d1ce4e5b9ULL;
        v = (v ^ (v >> 27))*0x94d049bb133111ebULL;
        return v ^ (v >> 31);
    }

public:

    explicit dsmcFastRng(uint64_t seed = 1)
    {
        reset(seed);
    }

    inline void reset(uint64_t seed)
    {
        s0_ = splitMix64(seed + 0x9e3779b97f4a7c15ULL);
        s1_ = splitMix64(seed + 0x9e3779b97f4a7c15ULL + 1ULL);
    }

    inline uint64_t next()
    {
        const uint64_t oldS0 = s0_;
        uint64_t oldS1 = s1_;
        oldS1 ^= oldS1 << 23;
        s1_ = oldS1 ^ oldS0 ^ (oldS1 >> 18) ^ (oldS0 >> 5);
        s0_ = oldS1;
        return s1_ + oldS0;
    }

    inline scalar sample01()
    {
        return scalar((next() >> 11)*0x1.0p-53);
    }

    inline label position(const label n)
    {
        if (n <= 1)
        {
            return 0;
        }

        const uint64_t bound = uint64_t(n);

        #if defined(__SIZEOF_INT128__)
        uint64_t x = next();
        __uint128_t m = __uint128_t(x)*__uint128_t(bound);
        uint64_t low = uint64_t(m);

        if (low < bound)
        {
            const uint64_t threshold = uint64_t(-bound)%bound;
            while (low < threshold)
            {
                x = next();
                m = __uint128_t(x)*__uint128_t(bound);
                low = uint64_t(m);
            }
        }

        return label(m >> 64);
        #else
        return label(next()%bound);
        #endif
    }

    static inline uint64_t collisionSeed
    (
        const label procI,
        const label threadI,
        const label timeIndex
    )
    {
        return uint64_t(procI)*1000000ULL
             + uint64_t(threadI)*10000ULL
             + uint64_t(timeIndex + 1);
    }

    static inline uint64_t moveSeed
    (
        const label procI,
        const label threadI,
        const label timeIndex
    )
    {
        uint64_t seed = 0x6d6f7665524e47ULL;
        seed ^= uint64_t(procI + 1)*0x9e3779b97f4a7c15ULL;
        seed ^= uint64_t(threadI + 1)*0xbf58476d1ce4e5b9ULL;
        seed ^= uint64_t(timeIndex + 1)*0x94d049bb133111ebULL;
        return seed;
    }
};


class alignas(64) dsmcPaddedFastRng
{
    dsmcFastRng rng_;

    static inline void* alignedAllocate(const std::size_t size)
    {
        void* ptr = nullptr;
        if (posix_memalign(&ptr, 64, size) != 0)
        {
            throw std::bad_alloc();
        }
        return ptr;
    }

public:

    explicit dsmcPaddedFastRng(uint64_t seed = 1)
    :
        rng_(seed)
    {}

    static void* operator new(const std::size_t size)
    {
        return alignedAllocate(size);
    }

    static void operator delete(void* ptr) noexcept
    {
        std::free(ptr);
    }

    static void* operator new[](const std::size_t) = delete;
    static void operator delete[](void*) = delete;

    inline void reset(uint64_t seed)
    {
        rng_.reset(seed);
    }

    inline dsmcFastRng& rng()
    {
        return rng_;
    }

    inline const dsmcFastRng& rng() const
    {
        return rng_;
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
