/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "ListOps.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::globalIndex::globalIndex()
{}


Foam::globalIndex::globalIndex(const Xfer<labelList>& offsets)
:
    offsets_(offsets)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::labelList& Foam::globalIndex::offsets()
{
    return offsets_;
}


inline Foam::label Foam::globalIndex::offset(const label procI) const
{
    return offsets_[procI];
}


inline Foam::label Foam::globalIndex::localSize(const label procI) const
{
    return offsets_[procI+1] - offsets_[procI];
}


inline Foam::label Foam::globalIndex::localSize() const
{
    return localSize(Pstream::myProcNo());
}


inline Foam::label Foam::globalIndex::size() const
{
    return offsets_.last();
}


inline Foam::label Foam::globalIndex::toGlobal
(
    const label procI,
    const label i
) const
{
    return i + offsets_[procI];
}


inline Foam::label Foam::globalIndex::toGlobal(const label i) const
{
    return toGlobal(Pstream::myProcNo(), i);
}


//- Is on local processor
inline bool Foam::globalIndex::isLocal(const label procI, const label i) const
{
    return i >= offsets_[procI] && i < offsets_[procI+1];
}


inline bool Foam::globalIndex::isLocal(const label i) const
{
    return isLocal(Pstream::myProcNo(), i);
}


inline Foam::label Foam::globalIndex::toLocal(const label procI, const label i)
const
{
    label localI = i - offsets_[procI];

    if (localI < 0 || i >= offsets_[procI+1])
    {
        FatalErrorIn("globalIndex::toLocal(const label, const label)")
            << "Global " << i << " does not belong on processor "
            << procI << endl << "Offsets:" << offsets_
            << abort(FatalError);
    }
    return localI;
}


inline Foam::label Foam::globalIndex::toLocal(const label i) const
{
    return toLocal(Pstream::myProcNo(), i);
}


inline Foam::label Foam::globalIndex::whichProcID(const label i) const
{
    if (i < 0 || i >= size())
    {
        FatalErrorIn("globalIndex::whichProcID(const label)")
            << "Global " << i << " does not belong on any processor."
            << " Offsets:" << offsets_
            << abort(FatalError);
    }

    return findLower(offsets_, i+1);
}


// ************************************************************************* //
