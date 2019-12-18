/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "PackedList.H"
#include "IOstreams.H"


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

#if (UINT_MAX == 0xFFFFFFFF)
// 32-bit counting, Hamming weight method
#define COUNT_PACKEDBITS(sum, x)                                               \
{                                                                              \
    x -= (x >> 1) & 0x55555555;                                                \
    x = (x & 0x33333333) + ((x >> 2) & 0x33333333);                            \
    sum += (((x + (x >> 4)) & 0x0F0F0F0F) * 0x01010101) >> 24;                 \
}
#elif (UINT_MAX == 0xFFFFFFFFFFFFFFFF)
// 64-bit counting, Hamming weight method
#define COUNT_PACKEDBITS(sum, x)                                               \
{                                                                              \
    x -= (x >> 1) & 0x5555555555555555;                                        \
    x = (x & 0x3333333333333333) + ((x >> 2) & 0x3333333333333333);            \
    sum += (((x + (x >> 4)) & 0x0F0F0F0F0F0F0F0F) * 0x0101010101010101) >> 56;\
}
#else
// Arbitrary number of bits, Brian Kernighan's method
    #define COUNT_PACKEDBITS(sum, x)    for (; x; ++sum) { x &= x - 1; }
#endif


template<unsigned nBits>
unsigned int Foam::PackedList<nBits>::count() const
{
    unsigned int c = 0;

    if (size_)
    {
        const label packLen = packedLength();
        for (label i = 0; i < packLen; ++i)
        {
            unsigned int bits = StorageList::operator[](i);
            COUNT_PACKEDBITS(c, bits);
        }
    }

    return c;
}


template<unsigned nBits>
bool Foam::PackedList<nBits>::trim()
{
    if (!size_)
    {
        return false;
    }

    const label oldSize = size_;
    for (label storeI = packedLength()-1; storeI >= 0; --storeI)
    {
        size_ = storeI * packing();
        unsigned int bits = StorageList::operator[](storeI);

        // found some bits
        if (bits)
        {
            while (bits)
            {
                bits >>= nBits;
                ++size_;
            }
            break;
        }
    }

    return (size_ != oldSize);
}


template<unsigned nBits>
void Foam::PackedList<nBits>::flip()
{
    if (!size_)
    {
        return;
    }

    // mask value for complete segments
    const unsigned int mask = maskLower(packing());

    const label packLen = packedLength();
    for (label i=0; i < packLen; ++i)
    {
        StorageList::operator[](i) = mask & ~StorageList::operator[](i);
    }

    // mask off the final partial segment
    {
        const unsigned int off = size_ % packing();
        if (off)
        {
            const unsigned int seg = size_ / packing();

            StorageList::operator[](seg) &= maskLower(off);
        }
    }
}


template<unsigned nBits>
Foam::labelList Foam::PackedList<nBits>::values() const
{
    labelList elems(size_);

    forAll(*this, i)
    {
        elems[i] = get(i);
    }

    return elems;
}


template<unsigned nBits>
Foam::Ostream& Foam::PackedList<nBits>::iteratorBase::printInfo
(
    Ostream& os
) const
{
    os  << "iterator<"  << label(nBits) << "> ["
        << this->index_ << "]"
        << " segment:"  << label(this->index_ / packing())
        << " offset:"   << label(this->index_ % packing())
        << " value:"    << this->get()
        << nl;

    return os;
}


template<unsigned nBits>
Foam::Ostream& Foam::PackedList<nBits>::printBits
(
    Ostream& os,
    const bool fullOutput
) const
{
    const label packLen = packedLength();

    // mask value for complete segments
    unsigned int mask = maskLower(packing());
    const label outputLen = fullOutput ? StorageList::size() : packLen;

    os  << "(\n";
    for (label i=0; i < outputLen; ++i)
    {
        const StorageType& rawBits = StorageList::operator[](i);

        // the final segment may not be full, modify mask accordingly
        if (i == packLen-1)
        {
            const unsigned int off = size_ % packing();

            if (off)
            {
                mask = maskLower(off);
            }
        }
        else if (i == packLen)
        {
            // no mask for unaddressed bit
            mask = 0u;
        }


        for (unsigned int testBit = (1u << max_bits()); testBit; testBit >>= 1)
        {
            if (mask & testBit)
            {
                // addressable region
                if (rawBits & testBit)
                {
                    os  << '1';
                }
                else
                {
                    os  << '-';
                }
            }
            else
            {
                if (rawBits & testBit)
                {
                    os  << '!';
                }
                else
                {
                    os  << '.';
                }
            }
        }
        os  << '\n';
    }
    os  << ")\n";

    return os;
}


template<unsigned nBits>
Foam::Ostream& Foam::PackedList<nBits>::printInfo
(
    Ostream& os,
    const bool fullOutput
) const
{
    os  << "PackedList<" << nBits << ">"
        << " max_value:" << max_value()
        << " packing:"   << packing() << nl
        << " count: "     << count() << nl
        << " size/capacity: " << size_ << "/" << capacity() << nl
        << " storage/capacity: "
        << packedLength() << "/" << StorageList::size()
        << "\n";

    return printBits(os, fullOutput);
}


template<unsigned nBits>
Foam::Istream& Foam::PackedList<nBits>::read(Istream& is)
{
    PackedList<nBits>& lst = *this;

    lst.clear();
    is.fatalCheck("PackedList<nBits>::read(Istream&)");

    token firstTok(is);
    is.fatalCheck
    (
        "PackedList<nBits>::read(Istream&) : "
        "reading first token"
    );

    if (firstTok.isLabel())
    {
        const label sz = firstTok.labelToken();

        // Set list length to that read
        lst.resize(sz);

        // Read list contents depending on data format
        if (is.format() == IOstream::ASCII)
        {
            // Read beginning of contents
            const char delimiter = is.readBeginList("PackedList<nBits>");

            if (sz)
            {
                if (delimiter == token::BEGIN_LIST)
                {
                    for (label i=0; i<sz; ++i)
                    {
                        lst[i] = lst.readValue(is);

                        is.fatalCheck
                        (
                            "PackedList<nBits>::read(Istream&) : "
                            "reading entry"
                        );
                    }
                }
                else if (delimiter == token::BEGIN_BLOCK)
                {
                    // assign for all entries
                    lst = lst.readValue(is);

                    is.fatalCheck
                    (
                        "PackedList<nBits>::read(Istream&) : "
                        "reading the single entry"
                    );
                }
                else
                {
                    FatalIOErrorInFunction(is)
                        << "incorrect list token, expected '(' or '{', found "
                        << firstTok.info()
                        << exit(FatalIOError);
                }
            }

            // Read end of contents
            is.readEndList("PackedList<nBits>");
        }
        else
        {
            if (sz)
            {
                is.read
                (
                    reinterpret_cast<char*>(lst.storage().data()),
                    lst.byteSize()
                );

                is.fatalCheck
                (
                    "PackedList<nBits>::read(Istream&) : "
                    "reading the binary block"
                );
            }
        }
    }
    else if (firstTok.isPunctuation())
    {
        if (firstTok.pToken() == token::BEGIN_LIST)
        {
            token nextTok(is);
            is.fatalCheck("PackedList<nBits>::read(Istream&)");

            while
            (
                !(   nextTok.isPunctuation()
                  && nextTok.pToken() == token::END_LIST
                 )
            )
            {
                is.putBack(nextTok);
                lst.append(lst.readValue(is));

                is  >> nextTok;
                is.fatalCheck("PackedList<nBits>::read(Istream&)");
            }
        }
        else if (firstTok.pToken() == token::BEGIN_BLOCK)
        {
            token nextTok(is);
            is.fatalCheck("PackedList<nBits>::read(Istream&)");

            while
            (
                !(   nextTok.isPunctuation()
                  && nextTok.pToken() == token::END_BLOCK
                 )
            )
            {
                is.putBack(nextTok);
                lst.setPair(is);

                is  >> nextTok;
                is.fatalCheck("PackedList<nBits>::read(Istream&)");
            }
        }
        else
        {
            FatalIOErrorInFunction(is)
                << "incorrect first token, expected '(', found "
                << firstTok.info()
                << exit(FatalIOError);
        }
    }
    else
    {
        FatalIOErrorInFunction(is)
            << "incorrect first token, expected <int>, '(' or '{', found "
            << firstTok.info()
            << exit(FatalIOError);
    }

    return is;
}


template<unsigned nBits>
Foam::Ostream& Foam::PackedList<nBits>::write
(
    Ostream& os,
    const bool indexedOutput
) const
{
    const PackedList<nBits>& lst = *this;
    const label sz = lst.size();

    // Write list contents depending on data format
    if (os.format() == IOstream::ASCII)
    {
        bool uniform = false;

        if (sz > 1 && !indexedOutput)
        {
            uniform = true;

            forAll(lst, i)
            {
                if (lst[i] != lst[0])
                {
                    uniform = false;
                    break;
                }
            }
        }

        if (uniform)
        {
            // uniform values:
            os  << sz << token::BEGIN_BLOCK << lst[0] << token::END_BLOCK;
        }
        else if (indexedOutput)
        {
            // indexed output
            os  << nl << token::BEGIN_BLOCK << nl;

            for
            (
                typename PackedList<nBits>::const_iterator iter = lst.cbegin();
                iter != lst.cend();
                ++iter
            )
            {
                if (iter.writeIfSet(os))
                {
                    os  << nl;
                }
            }

            os  << token::END_BLOCK << nl;
        }
        else if (sz < 11)
        {
            // short list:
            os  << sz << token::BEGIN_LIST;
            forAll(lst, i)
            {
                if (i)
                {
                    os  << token::SPACE;
                }
                os  << lst[i];
            }
            os  << token::END_LIST;
        }
        else
        {
            // longer list:
            os  << nl << sz << nl << token::BEGIN_LIST;
            forAll(lst, i)
            {
                os  << nl << lst[i];
            }
            os  << nl << token::END_LIST << nl;
        }
    }
    else
    {
        os  << nl << sz << nl;
        if (sz)
        {
            os.write
            (
                reinterpret_cast<const char*>(lst.storage().cdata()),
                lst.byteSize()
            );
        }
    }

    return os;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<unsigned nBits>
void Foam::PackedList<nBits>::operator=(const PackedList<nBits>& lst)
{
    StorageList::operator=(lst);
    size_ = lst.size();
}


template<unsigned nBits>
void Foam::PackedList<nBits>::operator=(PackedList<nBits>&& lst)
{
    transfer(lst);
}


template<unsigned nBits>
void Foam::PackedList<nBits>::operator=(const labelUList& lst)
{
    setCapacity(lst.size());
    size_ = lst.size();

    forAll(lst, i)
    {
        set(i, lst[i]);
    }
}


template<unsigned nBits>
void Foam::PackedList<nBits>::operator=(const UIndirectList<label>& lst)
{
    setCapacity(lst.size());
    size_ = lst.size();

    forAll(lst, i)
    {
        set(i, lst[i]);
    }
}


// * * * * * * * * * * * * * * * IOstream Functions  * * * * * * * * * * * * //

template<unsigned nBits>
void Foam::writeEntry(Ostream& os, const PackedList<nBits>& l)
{
    os << l;
}


// * * * * * * * * * * * * * *  Friend Operators * * * * * * * * * * * * * * //

template<unsigned nBits>
Foam::Istream& Foam::operator>>(Istream& is, PackedList<nBits>& lst)
{
    return lst.read(is);
}


template<unsigned nBits>
Foam::Ostream& Foam::operator<<(Ostream& os, const PackedList<nBits>& lst)
{
    return lst.write(os, false);
}


// ************************************************************************* //
