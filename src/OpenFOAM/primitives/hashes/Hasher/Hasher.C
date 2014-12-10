/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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

Description
    Hashing functions, mostly from Bob Jenkins
\*---------------------------------------------------------------------------*/

#include "Hasher.H"
#include "HasherInt.H"

#if defined (__GLIBC__)
#  include <endian.h>
#endif

// Left-rotate a 32-bit value and carry by nBits
#define bitRotateLeft(x, nBits)  (((x) << (nBits)) | ((x) >> (32 - (nBits))))


// ----------------------------------------------------------------------------
// lookup3.c, by Bob Jenkins, May 2006, Public Domain.
//
// These are functions for producing 32-bit hashes for hash table lookup.
// hashword(), hashlittle(), hashlittle2(), hashbig(), mix(), and final()
// are externally useful functions.  Routines to test the hash are included
// if SELF_TEST is defined.  You can use this free for any purpose.  It's in
// the public domain.  It has no warranty.
//
// You probably want to use hashlittle().  hashlittle() and hashbig()
// hash byte arrays.  hashlittle() is is faster than hashbig() on
// little-endian machines.  Intel and AMD are little-endian machines.
// On second thought, you probably want hashlittle2(), which is identical to
// hashlittle() except it returns two 32-bit hashes for the price of one.
// You could implement hashbig2() if you wanted but I haven't bothered here.
//
// If you want to find a hash of, say, exactly 7 integers, do
//   a = i1;  b = i2;  c = i3;
//   mix(a,b,c);
//   a += i4; b += i5; c += i6;
//   mix(a,b,c);
//   a += i7;
//   final(a,b,c);
// then use c as the hash value.  If you have a variable length array of
// 4-byte integers to hash, use hashword().  If you have a byte array (like
// a character string), use hashlittle().  If you have several byte arrays, or
// a mix of things, see the comments above hashlittle().
//
// Why is this so big?  I read 12 bytes at a time into 3 4-byte integers,
// then mix those integers.  This is fast (you can do a lot more thorough
// mixing with 12*3 instructions on 3 integers than you can with 3 instructions
// on 1 byte), but shoehorning those bytes into integers efficiently is messy.
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// mix -- mix 3 32-bit values reversibly.
//
// This is reversible, so any information in (a,b,c) before mix_hash() is
// still in (a,b,c) after mix_hash().
//
// If four pairs of (a,b,c) inputs are run through mix_hash(), or through
// mix_hash() in reverse, there are at least 32 bits of the output that
// are sometimes the same for one pair and different for another pair.
// This was tested for:
// * pairs that differed by one bit, by two bits, in any combination
//   of top bits of (a,b,c), or in any combination of bottom bits of
//   (a,b,c).
// * "differ" is defined as +, -, ^, or ~^.  For + and -, I transformed
//   the output delta to a Gray code (a^(a>>1)) so a string of 1's (as
//   is commonly produced by subtraction) look like a single 1-bit
//   difference.
// * the base values were pseudorandom, all zero but one bit set, or
//   all zero plus a counter that starts at zero.
//
// Some k values for my "a-=c; a^=rot(c,k); c+=b;" arrangement that
// satisfy this are
//     4  6  8 16 19  4
//     9 15  3 18 27 15
//    14  9  3  7 17  3
// Well, "9 15 3 18 27 15" didn't quite get 32 bits diffing
// for "differ" defined as + with a one-bit base and a two-bit delta.  I
// used http://burtleburtle.net/bob/hash/avalanche.html to choose
// the operations, constants, and arrangements of the variables.
//
// This does not achieve avalanche.  There are input bits of (a,b,c)
// that fail to affect some output bits of (a,b,c), especially of a.  The
// most thoroughly mixed value is c, but it doesn't really even achieve
// avalanche in c.
//
// This allows some parallelism.  Read-after-writes are good at doubling
// the number of bits affected, so the goal of mixing pulls in the opposite
// direction as the goal of parallelism.  I did what I could.  Rotates
// seem to cost as much as shifts on every machine I could lay my hands
// on, and rotates are much kinder to the top and bottom bits, so I used
// rotates.
// ----------------------------------------------------------------------------

#define bitMixer(a, b, c)                                                     \
    {                                                                         \
        a -= c; a ^= bitRotateLeft(c, 4); c += b;                             \
        b -= a; b ^= bitRotateLeft(a, 6); a += c;                             \
        c -= b; c ^= bitRotateLeft(b, 8); b += a;                             \
        a -= c; a ^= bitRotateLeft(c,16); c += b;                             \
        b -= a; b ^= bitRotateLeft(a,19); a += c;                             \
        c -= b; c ^= bitRotateLeft(b, 4); b += a;                             \
    }


// ----------------------------------------------------------------------------
// final -- final mixing of 3 32-bit values (a,b,c) into c
//
// Pairs of (a,b,c) values differing in only a few bits will usually
// produce values of c that look totally different.  This was tested for
// * pairs that differed by one bit, by two bits, in any combination
//   of top bits of (a,b,c), or in any combination of bottom bits of
//   (a,b,c).
// * "differ" is defined as +, -, ^, or ~^.  For + and -, I transformed
//   the output delta to a Gray code (a^(a>>1)) so a string of 1's (as
//   is commonly produced by subtraction) look like a single 1-bit
//   difference.
// * the base values were pseudorandom, all zero but one bit set, or
//   all zero plus a counter that starts at zero.
//
// These constants passed:
//  14 11 25 16 4 14 24
//  12 14 25 16 4 14 24
// and these came close:
//   4  8 15 26 3 22 24
//  10  8 15 26 3 22 24
//  11  8 15 26 3 22 24
// ----------------------------------------------------------------------------

#define bitMixerFinal(a, b, c)                                                \
    {                                                                         \
        c ^= b; c -= bitRotateLeft(b, 14);                                    \
        a ^= c; a -= bitRotateLeft(c, 11);                                    \
        b ^= a; b -= bitRotateLeft(a, 25);                                    \
        c ^= b; c -= bitRotateLeft(b, 16);                                    \
        a ^= c; a -= bitRotateLeft(c, 4);                                     \
        b ^= a; b -= bitRotateLeft(a, 14);                                    \
        c ^= b; c -= bitRotateLeft(b, 24);                                    \
    }


// * * * * * * * * * * * * * * Static Functions  * * * * * * * * * * * * * * //

// ----------------------------------------------------------------------------
// hashlittle() -- hash a variable-length key into a 32-bit value
//   k       : the key (the unaligned variable-length array of bytes)
//   length  : the length of the key, counting by bytes
//   initval : can be any 4-byte value
// Returns a 32-bit value.  Every bit of the key affects every bit of
// the return value.  Two keys differing by one or two bits will have
// totally different hash values.
//
// The best hash table sizes are powers of 2.  There is no need to do
// mod a prime (mod is sooo slow!).  If you need less than 32 bits,
// use a bitmask.  For example, if you need only 10 bits, do
//   h = (h & hashmask(10));
// In which case, the hash table should have hashsize(10) elements.
//
// If you are hashing n strings (uint8_t **)k, do it like this:
//   for (i=0, h=0; i<n; ++i) h = hashlittle( k[i], len[i], h);
//
// By Bob Jenkins, 2006.  bob_jenkins@burtleburtle.net.  You may use this
// code any way you wish, private, educational, or commercial.  It's free.
//
// Use for hash table lookup, or anything where one collision in 2^^32 is
// acceptable.  Do NOT use for cryptographic purposes.
// ----------------------------------------------------------------------------

//- specialized little-endian code
#if !defined (__BYTE_ORDER) || (__BYTE_ORDER == __LITTLE_ENDIAN)
static unsigned jenkins_hashlittle
(
    const void *key,
    size_t length,
    unsigned initval
)
{
    uint32_t a, b, c;
    union { const void *ptr; size_t i; } u; // to cast key to (size_t) happily

    // Set up the internal state
    a = b = c = 0xdeadbeef + static_cast<uint32_t>(length) + initval;

    u.ptr = key;
    if ((u.i & 0x3) == 0)
    {
        // 32-bit chunks
        const uint32_t *k = reinterpret_cast<const uint32_t*>(key);

        // all but last block: aligned reads and affect 32 bits of (a,b,c)
        while (length > 12)
        {
            a += k[0];
            b += k[1];
            c += k[2];
            bitMixer(a,b,c);
            length -= 12;
            k += 3;
        }

        // handle the last (probably partial) block byte-wise
        const uint8_t *k8 = reinterpret_cast<const uint8_t*>(k);
        switch (length)
        {
            case 12: c += k[2]; b += k[1]; a += k[0]; break;
            case 11: c += static_cast<uint32_t>(k8[10]) << 16; // fall through
            case 10: c += static_cast<uint32_t>(k8[9])  << 8;  // fall through
            case 9 : c += k8[8];                               // fall through
            case 8 : b += k[1]; a += k[0]; break;
            case 7 : b += static_cast<uint32_t>(k8[6]) << 16;  // fall through
            case 6 : b += static_cast<uint32_t>(k8[5]) << 8;   // fall through
            case 5 : b += k8[4];                               // fall through
            case 4 : a += k[0]; break;
            case 3 : a += static_cast<uint32_t>(k8[2]) << 16;  // fall through
            case 2 : a += static_cast<uint32_t>(k8[1]) << 8;   // fall through
            case 1 : a += k8[0]; break;
            case 0 : return c;  // zero-length requires no mixing
        }
    }
    else if ((u.i & 0x1) == 0)
    {
        // 16-bit chunks
        const uint16_t *k = reinterpret_cast<const uint16_t*>(key);

        // all but last block: aligned reads and different mixing
        while (length > 12)
        {
            a += k[0] + (static_cast<uint32_t>(k[1]) << 16);
            b += k[2] + (static_cast<uint32_t>(k[3]) << 16);
            c += k[4] + (static_cast<uint32_t>(k[5]) << 16);
            bitMixer(a,b,c);
            length -= 12;
            k += 6;
        }

        // handle the last (probably partial) block
        const uint8_t *k8 = reinterpret_cast<const uint8_t*>(k);
        switch (length)
        {
            case 12:
                c += k[4] + (static_cast<uint32_t>(k[5]) << 16);
                b += k[2] + (static_cast<uint32_t>(k[3]) << 16);
                a += k[0] + (static_cast<uint32_t>(k[1]) << 16);
                break;
            case 11:
                c += static_cast<uint32_t>(k8[10]) << 16;
                // fall through
            case 10:
                c += k[4];
                b += k[2] + (static_cast<uint32_t>(k[3]) << 16);
                a += k[0] + (static_cast<uint32_t>(k[1]) << 16);
                break;
            case 9 :
                c += k8[8];
                // fall through
            case 8 :
                b += k[2] + (static_cast<uint32_t>(k[3]) << 16);
                a += k[0] + (static_cast<uint32_t>(k[1]) << 16);
                break;
            case 7 :
                b += static_cast<uint32_t>(k8[6]) << 16;
                // fall through
            case 6 :
                b += k[2];
                a += k[0] + (static_cast<uint32_t>(k[1]) << 16);
                break;
            case 5 :
                b += k8[4];
                // fall through
            case 4 :
                a += k[0] + (static_cast<uint32_t>(k[1]) << 16);
                break;
            case 3 :
                a += static_cast<uint32_t>(k8[2]) << 16;
                // fall through
            case 2 :
                a += k[0];
                break;
            case 1 :
                a += k8[0];
                break;
            case 0 : return c;     // zero-length requires no mixing
        }
    }
    else
    {
        const uint8_t *k = reinterpret_cast<const uint8_t*>(key);

        // all but the last block: affect some 32 bits of (a,b,c)
        while (length > 12)
        {
            a += k[0];
            a += static_cast<uint32_t>(k[1]) << 8;
            a += static_cast<uint32_t>(k[2]) << 16;
            a += static_cast<uint32_t>(k[3]) << 24;
            b += k[4];
            b += static_cast<uint32_t>(k[5]) << 8;
            b += static_cast<uint32_t>(k[6]) << 16;
            b += static_cast<uint32_t>(k[7]) << 24;
            c += k[8];
            c += static_cast<uint32_t>(k[9])  << 8;
            c += static_cast<uint32_t>(k[10]) << 16;
            c += static_cast<uint32_t>(k[11]) << 24;

            bitMixer(a,b,c);
            length -= 12;
            k += 12;
        }

        // last block: affect all 32 bits of (c)
        switch (length) // most case statements fall through
        {
            case 12: c += static_cast<uint32_t>(k[11]) << 24;
            case 11: c += static_cast<uint32_t>(k[10]) << 16;
            case 10: c += static_cast<uint32_t>(k[9]) << 8;
            case 9 : c += k[8];

            case 8 : b += static_cast<uint32_t>(k[7]) << 24;
            case 7 : b += static_cast<uint32_t>(k[6]) << 16;
            case 6 : b += static_cast<uint32_t>(k[5]) << 8;
            case 5 : b += k[4];

            case 4 : a += static_cast<uint32_t>(k[3]) << 24;
            case 3 : a += static_cast<uint32_t>(k[2]) << 16;
            case 2 : a += static_cast<uint32_t>(k[1]) << 8;
            case 1 : a += k[0];
                break;

            case 0 : return c;
        }
    }

    bitMixerFinal(a,b,c);
    return c;
}
#endif




// ----------------------------------------------------------------------------
// hashbig():
// This is the same as hashword() on big-endian machines.  It is different
// from hashlittle() on all machines.  hashbig() takes advantage of
// big-endian byte ordering.
// ----------------------------------------------------------------------------
// specialized big-endian code
#if !defined (__BYTE_ORDER) || (__BYTE_ORDER == __BIG_ENDIAN)
static unsigned jenkins_hashbig
(
    const void *key,
    size_t length,
    unsigned initval
)
{
    uint32_t a, b, c;
    union { const void *ptr; size_t i; } u; // to cast key to (size_t) happily

    // Set up the internal state
    a = b = c = 0xdeadbeef + static_cast<uint32_t>(length) + initval;

    u.ptr = key;
    if ((u.i & 0x3) == 0)
    {
        // 32-bit chunks
        const uint32_t *k = reinterpret_cast<const uint32_t*>(key);

        // all but last block: aligned reads and affect 32 bits of (a,b,c)
        while (length > 12)
        {
            a += k[0];
            b += k[1];
            c += k[2];
            bitMixer(a,b,c);
            length -= 12;
            k += 3;
        }

        // handle the last (probably partial) block byte-wise
        const uint8_t *k8 = reinterpret_cast<const uint8_t*>(k);

        switch (length) // most the case statements fall through
        {
            case 12: c += k[2]; b += k[1]; a += k[0]; break;
            case 11: c += static_cast<uint32_t>(k8[10]) << 8; // fall through
            case 10: c += static_cast<uint32_t>(k8[9]) << 16; // fall through
            case 9 : c += static_cast<uint32_t>(k8[8]) << 24; // fall through
            case 8 : b += k[1]; a += k[0]; break;
            case 7 : b += static_cast<uint32_t>(k8[6]) << 8;  // fall through
            case 6 : b += static_cast<uint32_t>(k8[5]) << 16; // fall through
            case 5 : b += static_cast<uint32_t>(k8[4]) << 24; // fall through
            case 4 : a += k[0]; break;
            case 3 : a += static_cast<uint32_t>(k8[2]) << 8;  // fall through
            case 2 : a += static_cast<uint32_t>(k8[1]) << 16; // fall through
            case 1 : a += static_cast<uint32_t>(k8[0]) << 24; break;
            case 0 : return c;
        }
    }
    else
    {
        // need to read the key one byte at a time
        const uint8_t *k = reinterpret_cast<const uint8_t*>(key);

        // all but the last block: affect some 32 bits of (a,b,c)
        while (length > 12)
        {
            a += static_cast<uint32_t>(k[0]) << 24;
            a += static_cast<uint32_t>(k[1]) << 16;
            a += static_cast<uint32_t>(k[2]) << 8;
            a += static_cast<uint32_t>(k[3]);
            b += static_cast<uint32_t>(k[4]) << 24;
            b += static_cast<uint32_t>(k[5]) << 16;
            b += static_cast<uint32_t>(k[6]) << 8;
            b += static_cast<uint32_t>(k[7]);
            c += static_cast<uint32_t>(k[8]) << 24;
            c += static_cast<uint32_t>(k[9]) << 16;
            c += static_cast<uint32_t>(k[10]) << 8;
            c += static_cast<uint32_t>(k[11]);

            bitMixer(a,b,c);
            length -= 12;
            k += 12;
        }

        // last block: affect all 32 bits of (c)
        switch (length) // the case statements fall through
        {
            case 12: c += k[11];
            case 11: c += static_cast<uint32_t>(k[10]) << 8;
            case 10: c += static_cast<uint32_t>(k[9]) << 16;
            case 9 : c += static_cast<uint32_t>(k[8]) << 24;
            case 8 : b += k[7];
            case 7 : b += static_cast<uint32_t>(k[6]) << 8;
            case 6 : b += static_cast<uint32_t>(k[5]) << 16;
            case 5 : b += static_cast<uint32_t>(k[4]) << 24;
            case 4 : a += k[3];
            case 3 : a += static_cast<uint32_t>(k[2]) << 8;
            case 2 : a += static_cast<uint32_t>(k[1]) << 16;
            case 1 : a += static_cast<uint32_t>(k[0]) << 24;
                break;
            case 0 : return c;
        }
    }

    bitMixerFinal(a,b,c);
    return c;
}
#endif


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //


unsigned Foam::Hasher
(
    const void *key,
    size_t length,
    unsigned initval
)
{
#ifdef __BYTE_ORDER
# if (__BYTE_ORDER == __BIG_ENDIAN)
    return jenkins_hashbig(key, length, initval);
# else
    return jenkins_hashlittle(key, length, initval);
# endif
#else
    // endian-ness not known at compile-time: runtime endian test
    const short endianTest = 0x0100;

    // yields 0x01 for big endian
    if (*(reinterpret_cast<const char*>(&endianTest)))
    {
        return jenkins_hashbig(key, length, initval);
    }
    else
    {
        return jenkins_hashlittle(key, length, initval);
    }
#endif
}


// ----------------------------------------------------------------------------
//  This works on all machines.  To be useful, it requires
//  -- that the key be an array of uint32_t's, and
//  -- that the length be the number of uint32_t's in the key
//
//  The function hashword() is identical to hashlittle() on little-endian
//  machines, and identical to hashbig() on big-endian machines,
//  except that the length has to be measured in uint32_ts rather than in
//  bytes.  hashlittle() is more complicated than hashword() only because
//  hashlittle() has to dance around fitting the key bytes into registers.
// ----------------------------------------------------------------------------
unsigned Foam::HasherInt
(
    const uint32_t *k,
    size_t length,
    unsigned seed
)
{
    uint32_t a, b, c;

    // Set up the internal state
    a = b = c = 0xdeadbeef + (static_cast<uint32_t>(length) << 2) + seed;

    // handle most of the key
    while (length > 3)
    {
        a += k[0];
        b += k[1];
        c += k[2];
        bitMixer(a,b,c);
        length -= 3;
        k += 3;
    }

    // handle the last 3 uint32_t's
    switch (length)  // all case statements fall through
    {
        case 3 : c += k[2];
        case 2 : b += k[1];
        case 1 : a += k[0];
            bitMixerFinal(a,b,c);
        case 0 :  // case 0: nothing left to add
            break;
    }

    return c;
}


// ----------------------------------------------------------------------------
// hashword2() -- same as hashword(), but take two seeds and return two
// 32-bit values.  pc and pb must both be non-null, and *pc and *pb must
// both be initialized with seeds.  If you pass in (*pb)==0, the output
// (*pc) will be the same as the return value from hashword().
// ----------------------------------------------------------------------------
unsigned Foam::HasherDual
(
    const uint32_t *k,
    size_t length,
    unsigned& hash1,  // IN: seed OUT: primary hash value
    unsigned& hash2   // IN: more seed OUT: secondary hash value
)
{
    uint32_t a, b, c;

    // Set up the internal state
    a = b = c = 0xdeadbeef + (static_cast<uint32_t>(length) << 2) + hash1;
    c += hash2;

    // handle most of the key
    while (length > 3)
    {
        a += k[0];
        b += k[1];
        c += k[2];
        bitMixer(a,b,c);
        length -= 3;
        k += 3;
    }

    // handle the last 3 uint32_t's
    switch (length)  // all case statements fall through
    {
        case 3 : c += k[2];
        case 2 : b += k[1];
        case 1 : a += k[0];
            bitMixerFinal(a,b,c);
        case 0 :  // case 0: nothing left to add
            break;
    }

    // report the result
    hash1 = c;
    hash2 = b;

    // return primary hash value
    return c;
}


// ************************************************************************* //
