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

#include "error.H"
#include "SLListBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::SLListBase::iterator Foam::SLListBase::endIter_
(
    const_cast<SLListBase&>(static_cast<const SLListBase&>(SLListBase()))
);

Foam::SLListBase::const_iterator Foam::SLListBase::endConstIter_
(
    static_cast<const SLListBase&>(SLListBase()),
    reinterpret_cast<const link*>(0)
);


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::SLListBase::insert(SLListBase::link* a)
{
    nElmts_++;

    if (last_)
    {
        a->next_ = last_->next_;
    }
    else
    {
        last_ = a;
    }

    last_->next_ = a;
}


void Foam::SLListBase::append(SLListBase::link* a)
{
    nElmts_++;

    if (last_)
    {
        a->next_ = last_->next_;
        last_ = last_->next_ = a;
    }
    else
    {
        last_ = a->next_ = a;
    }
}


Foam::SLListBase::link* Foam::SLListBase::removeHead()
{
    nElmts_--;

    if (last_ == 0)
    {
        FatalErrorIn("SLListBase::remove()")
            << "remove from empty list"
            << abort(FatalError);
    }

    SLListBase::link* f = last_->next_;

    if (f == last_)
    {
        last_ = 0;
    }
    else
    {
        last_->next_ = f->next_;
    }

    return f;
}


Foam::SLListBase::link* Foam::SLListBase::remove(SLListBase::link* it)
{
    SLListBase::iterator iter = begin();
    SLListBase::link *prev = &(*iter);

    if (it == prev)
    {
        return removeHead();
    }

    nElmts_--;

    for (++iter; iter != end(); ++iter)
    {
        SLListBase::link *p = &(*iter);

        if (p == it)
        {
            prev->next_ = p->next_;

            if (p == last_)
            {
                last_ = prev;
            }

            return it;
        }

        prev = p;
    }

    return 0;
}


// ************************************************************************* //
