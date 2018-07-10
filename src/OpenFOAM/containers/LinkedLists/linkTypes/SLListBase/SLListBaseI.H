/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    Base singly-linked list.

\*---------------------------------------------------------------------------*/

#include "error.H"

// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

inline Foam::SLListBase::link::link()
:
    next_(0)
{}


inline Foam::SLListBase::link::link(link* p)
:
    next_(p)
{}


inline Foam::SLListBase::SLListBase()
:
    last_(0),
    nElmts_(0)
{}


inline Foam::SLListBase::SLListBase(link* a)
:
    last_(a->next_ = a),
    nElmts_(1)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

inline Foam::SLListBase::~SLListBase()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::label Foam::SLListBase::size() const
{
    return nElmts_;
}


inline bool Foam::SLListBase::empty() const
{
    return !nElmts_;
}


inline Foam::SLListBase::link*
Foam::SLListBase::first()
{
    if (!nElmts_)
    {
        FatalErrorInFunction
            << "list is empty"
            << abort(FatalError);
    }
    return last_->next_;
}


inline const Foam::SLListBase::link*
Foam::SLListBase::first() const
{
    if (!nElmts_)
    {
        FatalErrorInFunction
            << "list is empty"
            << abort(FatalError);
    }
    return last_->next_;
}


inline Foam::SLListBase::link*
Foam::SLListBase::last()
{
    if (!nElmts_)
    {
        FatalErrorInFunction
            << "list is empty"
            << abort(FatalError);
    }
    return last_;
}


inline const Foam::SLListBase::link*
Foam::SLListBase::last() const
{
    if (!nElmts_)
    {
        FatalErrorInFunction
            << "list is empty"
            << abort(FatalError);
    }
    return last_;
}


inline void Foam::SLListBase::clear()
{
    last_ = 0;
    nElmts_ = 0;
}


inline void Foam::SLListBase::transfer(SLListBase& lst)
{
    last_   = lst.last_;
    nElmts_ = lst.nElmts_;

    lst.clear();
}


inline Foam::SLListBase::link* Foam::SLListBase::remove
(
    SLListBase::iterator& it
)
{
    return remove(it.curElmt_);
}


// * * * * * * * * * * * * * * * STL iterator  * * * * * * * * * * * * * * * //

inline Foam::SLListBase::iterator::iterator(SLListBase& s, link* elmt)
:
    curList_(s),
    curElmt_(elmt),
    curLink_(*curElmt_)
{}


inline Foam::SLListBase::iterator::iterator(SLListBase& s)
:
    curList_(s),
    curElmt_(nullptr),
    curLink_()
{}


inline void Foam::SLListBase::iterator::operator=(const iterator& iter)
{
    curElmt_ = iter.curElmt_;
    curLink_ = iter.curLink_;
}


inline bool Foam::SLListBase::iterator::operator==(const iterator& iter) const
{
    return curElmt_ == iter.curElmt_;
}


inline bool Foam::SLListBase::iterator::operator!=(const iterator& iter) const
{
    return curElmt_ != iter.curElmt_;
}


inline Foam::SLListBase::link& Foam::SLListBase::iterator::operator*()
{
    return *curElmt_;
}


inline Foam::SLListBase::iterator& Foam::SLListBase::iterator::operator++()
{
    if (curElmt_ == curList_.last_ || curList_.last_ == 0)
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


inline Foam::SLListBase::iterator
Foam::SLListBase::iterator::operator++(int)
{
    iterator tmp = *this;
    ++*this;
    return tmp;
}


inline Foam::SLListBase::iterator
Foam::SLListBase::begin()
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


inline const Foam::SLListBase::iterator&
Foam::SLListBase::end()
{
    return endIter_;
}


// * * * * * * * * * * * * * * STL const_iterator  * * * * * * * * * * * * * //

inline Foam::SLListBase::const_iterator::const_iterator
(
    const SLListBase& s,
    const link* elmt
)
:
    curList_(s),
    curElmt_(elmt)
{}


inline Foam::SLListBase::const_iterator::const_iterator(const iterator& iter)
:
    curList_(iter.curList_),
    curElmt_(iter.curElmt_)
{}


inline void Foam::SLListBase::const_iterator::operator=
(
    const const_iterator& iter
)
{
    curElmt_ = iter.curElmt_;
}


inline bool Foam::SLListBase::const_iterator::operator==
(
    const const_iterator& iter
) const
{
    return curElmt_ == iter.curElmt_;
}


inline bool Foam::SLListBase::const_iterator::operator!=
(
    const const_iterator& iter
) const
{
    return curElmt_ != iter.curElmt_;
}


inline const Foam::SLListBase::link&
Foam::SLListBase::const_iterator::operator*()
{
    return *curElmt_;
}


inline Foam::SLListBase::const_iterator&
Foam::SLListBase::const_iterator::operator++()
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


inline Foam::SLListBase::const_iterator
Foam::SLListBase::const_iterator::operator++(int)
{
    const_iterator tmp = *this;
    ++*this;
    return tmp;
}


inline Foam::SLListBase::const_iterator
Foam::SLListBase::cbegin() const
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


inline const Foam::SLListBase::const_iterator&
Foam::SLListBase::cend() const
{
    return endConstIter_;
}


inline Foam::SLListBase::const_iterator
Foam::SLListBase::begin() const
{
    return this->cbegin();
}


inline const Foam::SLListBase::const_iterator&
Foam::SLListBase::end() const
{
    return endConstIter_;
}


// ************************************************************************* //
