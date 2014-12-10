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

\*---------------------------------------------------------------------------*/


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::labelRange::labelRange()
:
    start_(0),
    size_(0)
{}


inline Foam::labelRange::labelRange(const label start, const label size)
:
    start_(start),
    size_(size)
{
    // disallow invalid sizes
    if (size_ <= 0)
    {
        this->clear();
    }
}


// * * * * * * * * * * * * * * * * Iterators * * * * * * * * * * * * * * * * //

inline Foam::labelRange::const_iterator::const_iterator()
:
   range_(*reinterpret_cast< Foam::labelRange* >(0)),
   index_(-1)
{}


inline Foam::labelRange::const_iterator::const_iterator(const labelRange& range)
:
   range_(range),
   index_(range_.empty() ? -1 : 0)
{}


inline bool Foam::labelRange::const_iterator::operator==
(
    const const_iterator& iter
) const
{
    return (this->index_ == iter.index_);
}


inline bool Foam::labelRange::const_iterator::operator!=
(
    const const_iterator& iter
) const
{
    return !(this->operator==(iter));
}


inline Foam::label Foam::labelRange::const_iterator::operator*()
{
    return range_[index_];
}


inline Foam::label Foam::labelRange::const_iterator::operator()()
{
    return range_[index_];
}


inline Foam::labelRange::const_iterator&
Foam::labelRange::const_iterator::operator++()
{
    if (++index_ >= range_.size())
    {
        // equivalent to end iterator
        index_ = -1;
    }

    return *this;
}


inline Foam::labelRange::const_iterator
Foam::labelRange::const_iterator::operator++(int)
{
    const_iterator old = *this;
    this->operator++();
    return old;
}


inline Foam::labelRange::const_iterator Foam::labelRange::cbegin() const
{
    return const_iterator(*this);
}


inline const Foam::labelRange::const_iterator& Foam::labelRange::cend() const
{
    return endIter_;
}


inline Foam::labelRange::const_iterator Foam::labelRange::begin() const
{
    return const_iterator(*this);
}


inline const Foam::labelRange::const_iterator& Foam::labelRange::end() const
{
    return endIter_;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline void Foam::labelRange::clear()
{
    start_ = size_ = 0;
}


inline bool Foam::labelRange::empty() const
{
    return !size_;
}


inline Foam::label Foam::labelRange::size() const
{
    return size_;
}


inline Foam::label Foam::labelRange::first() const
{
    return start_;
}


inline Foam::label Foam::labelRange::last() const
{
    return start_ + size_ - 1;
}


inline bool Foam::labelRange::contains(const label value) const
{
    return value >= this->first() && value <= this->last();
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline Foam::label Foam::labelRange::operator[](const label i) const
{
    return start_ + i;
}


inline bool Foam::labelRange::operator<(const labelRange& rhs) const
{
    return this->first() < rhs.first();
}


inline bool Foam::labelRange::operator==(const labelRange& rhs) const
{
    return start_ == rhs.start_ && size_ == rhs.size_;
}


inline bool Foam::labelRange::operator!=(const labelRange& rhs) const
{
    return !(operator==(rhs));
}


// ************************************************************************* //
