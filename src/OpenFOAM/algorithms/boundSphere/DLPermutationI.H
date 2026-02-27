/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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

#include "DLPermutation.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class LabelList>
inline Foam::DLPermutation<LabelList>::const_iterator::const_iterator
(
    const DLPermutation<LabelList>& p,
    const label i
)
:
    p_(p),
    i_(i)
{}


template<class LabelList>
inline Foam::DLPermutation<LabelList>::DLPermutation
(
    const List<label>& permutation
)
:
    head_(-1),
    next_(permutation.size(), -1),
    tail_(-1),
    previous_(permutation.size(), -1)
{
    if (permutation.empty()) return;

    head_ = permutation.first();
    for (label i = 1; i < permutation.size(); ++ i)
    {
        next_[permutation[i - 1]] = permutation[i];
    }
    next_[permutation.last()] = -1;

    tail_ = permutation.last();
    for (label i = permutation.size() - 2; i >= 0; -- i)
    {
        previous_[permutation[i + 1]] = permutation[i];
    }
    previous_[permutation.first()] = -2;
}


template<class LabelList>
inline Foam::DLPermutation<LabelList>::DLPermutation
(
    const label n,
    randomGenerator& rndGen
)
:
    DLPermutation<LabelList>(rndGen.permutation(n))
{}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class LabelList>
inline Foam::label
Foam::DLPermutation<LabelList>::const_iterator::operator*() const
{
    return i_;
}


template<class LabelList>
inline typename Foam::DLPermutation<LabelList>::const_iterator&
Foam::DLPermutation<LabelList>::const_iterator::operator=
(
    const const_iterator& iter
)
{
    i_ = iter.i_;
    return *this;
}


template<class LabelList>
inline typename Foam::DLPermutation<LabelList>::const_iterator&
Foam::DLPermutation<LabelList>::const_iterator::operator++()
{
    i_ = p_.next_[i_];
    return *this;
}


template<class LabelList>
inline typename Foam::DLPermutation<LabelList>::const_iterator
Foam::DLPermutation<LabelList>::const_iterator::operator++(int)
{
    const_iterator iter(*this);
    this->operator++();
    return iter;
}


template<class LabelList>
inline typename Foam::DLPermutation<LabelList>::const_iterator&
Foam::DLPermutation<LabelList>::const_iterator::operator--()
{
    i_ = p_.previous_[i_];
    return *this;
}


template<class LabelList>
inline typename Foam::DLPermutation<LabelList>::const_iterator
Foam::DLPermutation<LabelList>::const_iterator::operator--(int)
{
    const_iterator iter(*this);
    this->operator--();
    return iter;
}


template<class LabelList>
inline bool Foam::DLPermutation<LabelList>::const_iterator::operator==
(
    const const_iterator& iter
) const
{
    return i_ == iter.i_;
}


template<class LabelList>
inline bool Foam::DLPermutation<LabelList>::const_iterator::operator!=
(
    const const_iterator& iter
) const
{
    return !(*this == iter);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class LabelList>
inline typename Foam::DLPermutation<LabelList>::const_iterator
Foam::DLPermutation<LabelList>::const_iterator::next() const
{
    if (i_ == -1)
    {
        return const_iterator(p_, p_.head_);
    }

    const_iterator iter(*this);
    ++ iter;
    return iter;
}


template<class LabelList>
inline typename Foam::DLPermutation<LabelList>::const_iterator
Foam::DLPermutation<LabelList>::const_iterator::previous() const
{
    if (i_ == -1)
    {
        return const_iterator(p_, p_.tail_);
    }

    const_iterator iter(*this);
    -- iter;
    return iter;
}


template<class LabelList>
inline Foam::List<Foam::label> Foam::DLPermutation<LabelList>::permutation
(
    const bool reverse
)
{
    List<label> result(next_.size());

    if (!reverse)
    {
        label i = 0;
        forAllConstIter(DLPermutation<LabelList>, *this, iter)
        {
            result[i ++] = *iter;
        }
    }
    else
    {
        label i = result.size();
        forAllConstIterReverse(DLPermutation<LabelList>, *this, iter)
        {
            result[-- i] = *iter;
        }
    }

    return result;
}


template<class LabelList>
inline typename Foam::DLPermutation<LabelList>::const_iterator
Foam::DLPermutation<LabelList>::begin() const
{
    return const_iterator(*this, head_);
}


template<class LabelList>
inline typename Foam::DLPermutation<LabelList>::const_iterator
Foam::DLPermutation<LabelList>::end() const
{
    return const_iterator(*this, -1);
}


template<class LabelList>
inline typename Foam::DLPermutation<LabelList>::const_iterator
Foam::DLPermutation<LabelList>::rbegin() const
{
    return const_iterator(*this, tail_);
}


template<class LabelList>
inline typename Foam::DLPermutation<LabelList>::const_iterator
Foam::DLPermutation<LabelList>::rend() const
{
    return const_iterator(*this, -2);
}


template<class LabelList>
void Foam::DLPermutation<LabelList>::append()
{
    const label i = next_.size();

    next_.append(-1);
    previous_.append(tail_);
    next_[tail_] = i;
    tail_ = i;
}


template<class LabelList>
void Foam::DLPermutation<LabelList>::prepend()
{
    const label i = next_.size();

    next_.append(head_);
    previous_.append(-2);
    previous_[head_] = i;
    head_ = i;
}


template<class LabelList>
void Foam::DLPermutation<LabelList>::raise(const label i)
{
    if (i == head_) return;

    const label i0 = previous_[i];
    const label i1 = next_[i];

    (i1 != -1 ? previous_[i1] : tail_) = previous_[i];
    previous_[i] = -2;
    previous_[head_] = i;

    next_[i0] = next_[i];
    next_[i] = head_;
    head_ = i;
}


template<class LabelList>
void Foam::DLPermutation<LabelList>::sink(const label i)
{
    if (i == tail_) return;

    const label i0 = previous_[i];
    const label i1 = next_[i];

    (i0 != -2 ? next_[i0] : head_) = next_[i];
    next_[i] = -1;
    next_[tail_] = i;

    previous_[i1] = previous_[i];
    previous_[i] = tail_;
    tail_ = i;
}


// ************************************************************************* //
