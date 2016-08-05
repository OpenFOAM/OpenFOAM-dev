/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "error.H"

// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

inline Foam::DLListBase::link::link()
:
    prev_(0),
    next_(0)
{}


inline Foam::DLListBase::DLListBase()
:
    first_(0),
    last_(0),
    nElmts_(0)
{}


inline Foam::DLListBase::DLListBase(link* a)
:
    first_(a),
    last_(a),
    nElmts_(1)
{
    a->prev_ = a;
    a->next_ = a;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

inline Foam::DLListBase::~DLListBase()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline bool Foam::DLListBase::link::registered() const
{
    return prev_ != 0 && next_ != 0;
}


inline void Foam::DLListBase::link::deregister()
{
    prev_ = 0;
    next_ = 0;
}


inline Foam::label Foam::DLListBase::size() const
{
    return nElmts_;
}


inline bool Foam::DLListBase::empty() const
{
    return !nElmts_;
}


inline Foam::DLListBase::link*
Foam::DLListBase::first()
{
    if (!nElmts_)
    {
        FatalErrorInFunction
            << "list is empty"
            << abort(FatalError);
    }
    return first_;
}


inline const Foam::DLListBase::link*
Foam::DLListBase::first() const
{
    if (!nElmts_)
    {
        FatalErrorInFunction
            << "list is empty"
            << abort(FatalError);
    }
    return first_;
}


inline Foam::DLListBase::link*
Foam::DLListBase::last()
{
    if (!nElmts_)
    {
        FatalErrorInFunction
            << "list is empty"
            << abort(FatalError);
    }
    return last_;
}


inline const Foam::DLListBase::link*
Foam::DLListBase::last() const
{
    if (!nElmts_)
    {
        FatalErrorInFunction
            << "list is empty"
            << abort(FatalError);
    }
    return last_;
}


inline void Foam::DLListBase::clear()
{
    first_ = 0;
    last_  = 0;
    nElmts_ = 0;
}


inline void Foam::DLListBase::transfer(DLListBase& lst)
{
    first_  = lst.first_;
    last_   = lst.last_;
    nElmts_ = lst.nElmts_;

    lst.clear();
}


inline Foam::DLListBase::link*
Foam::DLListBase::remove
(
    DLListBase::iterator& it
)
{
    return remove(it.curElmt_);
}


inline Foam::DLListBase::link*
Foam::DLListBase::replace
(
    DLListBase::iterator& oldIter,
    DLListBase::link* newLink
)
{
    return replace(oldIter.curElmt_, newLink);
}


// * * * * * * * * * * * * * * * STL iterator  * * * * * * * * * * * * * * * //

inline Foam::DLListBase::iterator::iterator(DLListBase& s, link* elmt)
:
    curList_(s),
    curElmt_(elmt),
    curLink_(*curElmt_)
{}


inline Foam::DLListBase::iterator::iterator(DLListBase& s)
:
    curList_(s),
    curElmt_(nullptr),
    curLink_()
{}


inline void Foam::DLListBase::iterator::operator=(const iterator& iter)
{
    curElmt_ = iter.curElmt_;
    curLink_ = iter.curLink_;
}


inline bool Foam::DLListBase::iterator::operator==(const iterator& iter) const
{
    return curElmt_ == iter.curElmt_;
}


inline bool Foam::DLListBase::iterator::operator!=(const iterator& iter) const
{
    return curElmt_ != iter.curElmt_;
}


inline Foam::DLListBase::link&
Foam::DLListBase::iterator::operator*()
{
    return *curElmt_;
}


inline Foam::DLListBase::iterator&
Foam::DLListBase::iterator::operator++()
{
    // Check if the curElmt_ is the last element (if it points to itself)
    // or if the list is empty because the last element may have been removed
    if (curLink_.next_ == curElmt_ || curList_.last_ == 0)
    {
        curElmt_ = 0;
    }
    else
    {
        curElmt_ = curLink_.next_;
        curLink_ = *curElmt_;
    }

    return *this;
}


inline Foam::DLListBase::iterator
Foam::DLListBase::iterator::operator++(int)
{
    iterator tmp = *this;
    ++*this;
    return tmp;
}


inline Foam::DLListBase::iterator
Foam::DLListBase::begin()
{
    if (size())
    {
        return iterator(*this, first());
    }
    else
    {
        return endIter_;
    }
}


inline const Foam::DLListBase::iterator& Foam::DLListBase::end()
{
    return endIter_;
}


// * * * * * * * * * * * * * * STL const_iterator  * * * * * * * * * * * * * //

inline Foam::DLListBase::const_iterator::const_iterator
(
    const DLListBase& s,
    const link* elmt
)
:
    curList_(s),
    curElmt_(elmt)
{}


inline Foam::DLListBase::const_iterator::const_iterator(const iterator& iter)
:
    curList_(iter.curList_),
    curElmt_(iter.curElmt_)
{}


inline void Foam::DLListBase::const_iterator::operator=
(
    const const_iterator& iter
)
{
    curElmt_ = iter.curElmt_;
}


inline bool Foam::DLListBase::const_iterator::operator==
(
    const const_iterator& iter
) const
{
    return curElmt_ == iter.curElmt_;
}


inline bool Foam::DLListBase::const_iterator::operator!=
(
    const const_iterator& iter
) const
{
    return curElmt_ != iter.curElmt_;
}


inline const Foam::DLListBase::link&
Foam::DLListBase::const_iterator::operator*()
{
    return *curElmt_;
}


inline Foam::DLListBase::const_iterator&
Foam::DLListBase::const_iterator::operator++()
{
    if (curElmt_ == curList_.last_)
    {
        curElmt_ = 0;
    }
    else
    {
        curElmt_ = curElmt_->next_;
    }

    return *this;
}


inline Foam::DLListBase::const_iterator
Foam::DLListBase::const_iterator::operator++(int)
{
    const_iterator tmp = *this;
    ++*this;
    return tmp;
}


inline Foam::DLListBase::const_iterator
Foam::DLListBase::cbegin() const
{
    if (size())
    {
        return const_iterator(*this, first());
    }
    else
    {
        return endConstIter_;
    }
}


inline const Foam::DLListBase::const_iterator&
Foam::DLListBase::cend() const
{
    return endConstIter_;
}


inline Foam::DLListBase::const_iterator
Foam::DLListBase::begin() const
{
    return this->cbegin();
}


inline const Foam::DLListBase::const_iterator&
Foam::DLListBase::end() const
{
    return endConstIter_;
}


// * * * * * * * * * * STL const_reverse_iterator  * * * * * * * * * * * * * //

inline Foam::DLListBase::const_reverse_iterator::const_reverse_iterator
(
    const DLListBase& s,
    const link* elmt
)
:
    curList_(s),
    curElmt_(elmt)
{}


inline void Foam::DLListBase::const_reverse_iterator::operator=
(
    const const_reverse_iterator& iter
)
{
    curElmt_ = iter.curElmt_;
}


inline bool Foam::DLListBase::const_reverse_iterator::operator==
(
    const const_reverse_iterator& iter
) const
{
    return curElmt_ == iter.curElmt_;
}


inline bool Foam::DLListBase::const_reverse_iterator::operator!=
(
    const const_reverse_iterator& iter
) const
{
    return curElmt_ != iter.curElmt_;
}


inline const Foam::DLListBase::link&
Foam::DLListBase::const_reverse_iterator::operator*()
{
    return *curElmt_;
}


inline Foam::DLListBase::const_reverse_iterator&
Foam::DLListBase::const_reverse_iterator::operator++()
{
    if (curElmt_ == curList_.first_)
    {
        curElmt_ = 0;
    }
    else
    {
        curElmt_ = curElmt_->prev_;
    }

    return *this;
}


inline Foam::DLListBase::const_reverse_iterator
Foam::DLListBase::const_reverse_iterator::operator++(int)
{
    const_reverse_iterator tmp = *this;
    ++*this;
    return tmp;
}


inline Foam::DLListBase::const_reverse_iterator
Foam::DLListBase::crbegin() const
{
    if (size())
    {
        return const_reverse_iterator(*this, last());
    }
    else
    {
        return endConstRevIter_;
    }
}


inline const Foam::DLListBase::const_reverse_iterator&
Foam::DLListBase::crend() const
{
    return endConstRevIter_;
}


inline Foam::DLListBase::const_reverse_iterator
Foam::DLListBase::rbegin() const
{
    return this->crbegin();
}


inline const Foam::DLListBase::const_reverse_iterator&
Foam::DLListBase::rend() const
{
    return endConstRevIter_;
}


// ************************************************************************* //
