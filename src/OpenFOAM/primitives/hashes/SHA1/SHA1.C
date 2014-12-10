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
    Functions to compute SHA1 message digest of files or memory blocks
    according to the NIST specification FIPS-180-1.

    Adapted from the gnulib implementation written by Scott G. Miller with
    credits to Robert Klep <robert@ilse.nl> -- Expansion function fix

    Copyright (C) 2000, 2001, 2003, 2004, 2005, 2006, 2008 Free Software
    Foundation, Inc.

\*---------------------------------------------------------------------------*/

#include "SHA1.H"
#include "IOstreams.H"

#include <cstring>

#if defined (__GLIBC__)
#  include <endian.h>
#endif


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

//! \cond fileScope
//  The bytes used to pad buffer to the next 64-byte boundary.
//  (RFC 1321, 3.1: Step 1)
static const unsigned char fillbuf[64] = { 0x80, 0 /* , 0, 0, ...  */ };
//! \endcond


// * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * * //

inline uint32_t Foam::SHA1::swapBytes(uint32_t n)
{
#ifdef __BYTE_ORDER
# if (__BYTE_ORDER == __BIG_ENDIAN)
    return n;
# else
    return
    (
        ((n) << 24)
      | (((n) & 0xff00) << 8)
      | (((n) >> 8) & 0xff00)
      | ((n) >> 24)
    );
# endif

#else

    const short x = 0x0100;

    // yields 0x01 for big endian
    if (*(reinterpret_cast<const char*>(&x)))
    {
        return n;
    }
    else
    {
        return
        (
            ((n) << 24)
          | (((n) & 0xff00) << 8)
          | (((n) >> 8) & 0xff00)
          | ((n) >> 24)
        );
    }
#endif
}


inline void
Foam::SHA1::set_uint32(unsigned char *cp, uint32_t v)
{
    memcpy(cp, &v, sizeof(uint32_t));
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //


void Foam::SHA1::processBytes(const void *data, size_t len)
{
    // already finalized, thus need to restart from nothing
    if (finalized_)
    {
        clear();
    }

    // complete filling of internal buffer
    if (bufLen_)
    {
        size_t remaining = bufLen_;
        size_t add =
        (
            sizeof(buffer_) - remaining > len
          ? len
          : sizeof(buffer_) - remaining
        );

        unsigned char* bufp = reinterpret_cast<unsigned char*>(buffer_);

        memcpy(&bufp[remaining], data, add);
        bufLen_ += add;

        if (bufLen_ > 64)
        {
            processBlock(buffer_, bufLen_ & ~63);

            bufLen_ &= 63;
            // The regions in the following copy operation do not
            // (cannot) overlap
            memcpy(buffer_, &bufp[(remaining + add) & ~63], bufLen_);
        }

        data = reinterpret_cast<const unsigned char*>(data) + add;
        len -= add;
    }

    // Process available complete blocks
//    if (len >= 64)
//    {
//#if !_STRING_ARCH_unaligned
//# define alignof(type) offsetof (struct { char c; type x; }, x)
//# define UNALIGNED_P(p) (((size_t) p) % alignof (uint32_t) != 0)
//        if (UNALIGNED_P (data))
//        {
//            while (len > 64)
            while (len >= 64)
            {
                processBlock(memcpy(buffer_, data, 64), 64);
                data = reinterpret_cast<const unsigned char*>(data) + 64;
                len -= 64;
            }
//        }
//        else
//#endif
//        {
//            processBlock(data, len & ~63);
//            data = reinterpret_cast<const unsigned char*>(data) + (len & ~63);
//            len &= 63;
//        }
//    }

    // Move remaining bytes in internal buffer.
    if (len > 0)
    {
        unsigned char* bufp = reinterpret_cast<unsigned char*>(buffer_);
        size_t remaining = bufLen_;

        memcpy (&bufp[remaining], data, len);
        remaining += len;
        if (remaining >= 64)
        {
            processBlock(buffer_, 64);
            remaining -= 64;
            memcpy(buffer_, &buffer_[16], remaining);
        }
        bufLen_ = remaining;
    }
}


// SHA1 round constants
#define K1 0x5a827999
#define K2 0x6ed9eba1
#define K3 0x8f1bbcdc
#define K4 0xca62c1d6

// Round functions.  Note that F2 is the same as F4.
#define F1(B,C,D) ( D ^ ( B & ( C ^ D ) ) )
#define F2(B,C,D) (B ^ C ^ D)
#define F3(B,C,D) ( ( B & C ) | ( D & ( B | C ) ) )
#define F4(B,C,D) (B ^ C ^ D)

// Process LEN bytes of BUFFER, it is assumed that LEN % 64 == 0.
// Most of this code comes from GnuPG's cipher/sha1.c

void
Foam::SHA1::processBlock(const void *data, size_t len)
{
    const uint32_t *words = reinterpret_cast<const uint32_t*>(data);
    size_t nwords = len / sizeof(uint32_t);
    const uint32_t *endp = words + nwords;

    // calculate with sixteen words of 32-bits
    uint32_t x[16];
    uint32_t a = hashsumA_;
    uint32_t b = hashsumB_;
    uint32_t c = hashsumC_;
    uint32_t d = hashsumD_;
    uint32_t e = hashsumE_;

    // First increment the byte count.
    // RFC 1321 specifies the possible length of the file up to 2^64 bits.
    // Here we only compute the number of bytes.  Do a double word increment.
    bufTotal_[0] += len;
    if (bufTotal_[0] < len)
    {
        ++bufTotal_[1];
    }

    // rotate left uint32_t by n bits
#define rol_uint32(x, nbits)  (((x) << (nbits)) | ((x) >> (32 - (nbits))))

#define M(I) ( tm = x[I & 0x0F] ^ x[(I-14) & 0x0F]                            \
               ^ x[(I-8) & 0x0F] ^ x[(I-3) & 0x0F]                            \
               , (x[I & 0x0F] = rol_uint32(tm, 1)) )


#define R(A,B,C,D,E,F,K,M)                                                    \
    do                                                                        \
    {                                                                         \
        E += rol_uint32(A, 5) + F(B, C, D) + K + M;                           \
        B = rol_uint32(B, 30);                                                \
    } while (0)

    while (words < endp)
    {
        uint32_t tm;
        for (int t = 0; t < 16; ++t)
        {
            x[t] = swapBytes(*words);
            ++words;
        }

        R( a, b, c, d, e, F1, K1, x[ 0] );
        R( e, a, b, c, d, F1, K1, x[ 1] );
        R( d, e, a, b, c, F1, K1, x[ 2] );
        R( c, d, e, a, b, F1, K1, x[ 3] );
        R( b, c, d, e, a, F1, K1, x[ 4] );
        R( a, b, c, d, e, F1, K1, x[ 5] );
        R( e, a, b, c, d, F1, K1, x[ 6] );
        R( d, e, a, b, c, F1, K1, x[ 7] );
        R( c, d, e, a, b, F1, K1, x[ 8] );
        R( b, c, d, e, a, F1, K1, x[ 9] );
        R( a, b, c, d, e, F1, K1, x[10] );
        R( e, a, b, c, d, F1, K1, x[11] );
        R( d, e, a, b, c, F1, K1, x[12] );
        R( c, d, e, a, b, F1, K1, x[13] );
        R( b, c, d, e, a, F1, K1, x[14] );
        R( a, b, c, d, e, F1, K1, x[15] );
        R( e, a, b, c, d, F1, K1, M(16) );
        R( d, e, a, b, c, F1, K1, M(17) );
        R( c, d, e, a, b, F1, K1, M(18) );
        R( b, c, d, e, a, F1, K1, M(19) );
        R( a, b, c, d, e, F2, K2, M(20) );
        R( e, a, b, c, d, F2, K2, M(21) );
        R( d, e, a, b, c, F2, K2, M(22) );
        R( c, d, e, a, b, F2, K2, M(23) );
        R( b, c, d, e, a, F2, K2, M(24) );
        R( a, b, c, d, e, F2, K2, M(25) );
        R( e, a, b, c, d, F2, K2, M(26) );
        R( d, e, a, b, c, F2, K2, M(27) );
        R( c, d, e, a, b, F2, K2, M(28) );
        R( b, c, d, e, a, F2, K2, M(29) );
        R( a, b, c, d, e, F2, K2, M(30) );
        R( e, a, b, c, d, F2, K2, M(31) );
        R( d, e, a, b, c, F2, K2, M(32) );
        R( c, d, e, a, b, F2, K2, M(33) );
        R( b, c, d, e, a, F2, K2, M(34) );
        R( a, b, c, d, e, F2, K2, M(35) );
        R( e, a, b, c, d, F2, K2, M(36) );
        R( d, e, a, b, c, F2, K2, M(37) );
        R( c, d, e, a, b, F2, K2, M(38) );
        R( b, c, d, e, a, F2, K2, M(39) );
        R( a, b, c, d, e, F3, K3, M(40) );
        R( e, a, b, c, d, F3, K3, M(41) );
        R( d, e, a, b, c, F3, K3, M(42) );
        R( c, d, e, a, b, F3, K3, M(43) );
        R( b, c, d, e, a, F3, K3, M(44) );
        R( a, b, c, d, e, F3, K3, M(45) );
        R( e, a, b, c, d, F3, K3, M(46) );
        R( d, e, a, b, c, F3, K3, M(47) );
        R( c, d, e, a, b, F3, K3, M(48) );
        R( b, c, d, e, a, F3, K3, M(49) );
        R( a, b, c, d, e, F3, K3, M(50) );
        R( e, a, b, c, d, F3, K3, M(51) );
        R( d, e, a, b, c, F3, K3, M(52) );
        R( c, d, e, a, b, F3, K3, M(53) );
        R( b, c, d, e, a, F3, K3, M(54) );
        R( a, b, c, d, e, F3, K3, M(55) );
        R( e, a, b, c, d, F3, K3, M(56) );
        R( d, e, a, b, c, F3, K3, M(57) );
        R( c, d, e, a, b, F3, K3, M(58) );
        R( b, c, d, e, a, F3, K3, M(59) );
        R( a, b, c, d, e, F4, K4, M(60) );
        R( e, a, b, c, d, F4, K4, M(61) );
        R( d, e, a, b, c, F4, K4, M(62) );
        R( c, d, e, a, b, F4, K4, M(63) );
        R( b, c, d, e, a, F4, K4, M(64) );
        R( a, b, c, d, e, F4, K4, M(65) );
        R( e, a, b, c, d, F4, K4, M(66) );
        R( d, e, a, b, c, F4, K4, M(67) );
        R( c, d, e, a, b, F4, K4, M(68) );
        R( b, c, d, e, a, F4, K4, M(69) );
        R( a, b, c, d, e, F4, K4, M(70) );
        R( e, a, b, c, d, F4, K4, M(71) );
        R( d, e, a, b, c, F4, K4, M(72) );
        R( c, d, e, a, b, F4, K4, M(73) );
        R( b, c, d, e, a, F4, K4, M(74) );
        R( a, b, c, d, e, F4, K4, M(75) );
        R( e, a, b, c, d, F4, K4, M(76) );
        R( d, e, a, b, c, F4, K4, M(77) );
        R( c, d, e, a, b, F4, K4, M(78) );
        R( b, c, d, e, a, F4, K4, M(79) );

        a = hashsumA_ += a;
        b = hashsumB_ += b;
        c = hashsumC_ += c;
        d = hashsumD_ += d;
        e = hashsumE_ += e;
    }
}


void Foam::SHA1::calcDigest(SHA1Digest& dig) const
{
    if (bufTotal_[0] || bufTotal_[1])
    {
        unsigned char *r = dig.v_;

        set_uint32(r + 0 * sizeof(uint32_t), swapBytes(hashsumA_));
        set_uint32(r + 1 * sizeof(uint32_t), swapBytes(hashsumB_));
        set_uint32(r + 2 * sizeof(uint32_t), swapBytes(hashsumC_));
        set_uint32(r + 3 * sizeof(uint32_t), swapBytes(hashsumD_));
        set_uint32(r + 4 * sizeof(uint32_t), swapBytes(hashsumE_));
    }
    else
    {
        // no data!
        dig.clear();
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


void Foam::SHA1::clear()
{
    hashsumA_ = 0x67452301;
    hashsumB_ = 0xefcdab89;
    hashsumC_ = 0x98badcfe;
    hashsumD_ = 0x10325476;
    hashsumE_ = 0xc3d2e1f0;

    bufTotal_[0] = bufTotal_[1] = 0;
    bufLen_ = 0;

    finalized_ = false;
}


bool Foam::SHA1::finalize()
{
    if (!finalized_)
    {
        finalized_ = true;

        // account for unprocessed bytes
        uint32_t bytes = bufLen_;
        size_t size = (bytes < 56 ? 64 : 128) / sizeof(uint32_t);

        // count remaining bytes.
        bufTotal_[0] += bytes;
        if (bufTotal_[0] < bytes)
        {
            ++bufTotal_[1];
        }

        // finalized, but no data!
        if (!bufTotal_[0] && !bufTotal_[1])
        {
            return false;
        }

        // place the 64-bit file length in *bits* at the end of the buffer.
        buffer_[size-2] = swapBytes((bufTotal_[1] << 3) | (bufTotal_[0] >> 29));
        buffer_[size-1] = swapBytes(bufTotal_[0] << 3);

        unsigned char* bufp = reinterpret_cast<unsigned char *>(buffer_);

        memcpy(&bufp[bytes], fillbuf, (size-2) * sizeof(uint32_t) - bytes);

        // Process remaining bytes
        processBlock(buffer_, size * sizeof(uint32_t));
    }

    return true;
}


Foam::SHA1Digest Foam::SHA1::digest() const
{
    SHA1Digest dig;

    if (finalized_)
    {
        calcDigest(dig);
    }
    else
    {
        // avoid disturbing our data - use a copy
        SHA1 sha(*this);
        if (sha.finalize())
        {
            sha.calcDigest(dig);
        }
    }

    return dig;
}


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

// void Foam::SHA1::operator=(const SHA1& rhs)
// {
//     // Check for assignment to self
//     if (this == &rhs)
//     {
//         FatalErrorIn("Foam::SHA1::operator=(const Foam::SHA1&)")
//             << "Attempted assignment to self"
//             << abort(FatalError);
//     }
// }


// ************************************************************************* //
