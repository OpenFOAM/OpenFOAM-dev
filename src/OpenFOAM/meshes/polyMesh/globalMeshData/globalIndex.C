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

\*---------------------------------------------------------------------------*/

#include "globalIndex.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::globalIndex::globalIndex
(
    const label localSize,
    const int tag,
    const label comm,
    const bool parallel
)
:
    offsets_(Pstream::nProcs(comm)+1)
{
    labelList localSizes(Pstream::nProcs(comm), 0);
    localSizes[Pstream::myProcNo(comm)] = localSize;
    if (parallel)
    {
        Pstream::gatherList(localSizes, tag, comm);
        Pstream::scatterList(localSizes, tag, comm);
    }

    label offset = 0;
    offsets_[0] = 0;
    for (label proci = 0; proci < Pstream::nProcs(comm); proci++)
    {
        label oldOffset = offset;
        offset += localSizes[proci];

        if (offset < oldOffset)
        {
            FatalErrorInFunction
                << "Overflow : sum of sizes " << localSizes
                << " exceeds capability of label (" << labelMax
                << "). Please recompile with larger datatype for label."
                << exit(FatalError);
        }
        offsets_[proci+1] = offset;
    }
}


Foam::globalIndex::globalIndex(const label localSize)
:
    offsets_(Pstream::nProcs()+1)
{
    labelList localSizes(Pstream::nProcs(), 0);
    localSizes[Pstream::myProcNo()] = localSize;
    Pstream::gatherList(localSizes, Pstream::msgType());
    Pstream::scatterList(localSizes, Pstream::msgType());

    label offset = 0;
    offsets_[0] = 0;
    for (label proci = 0; proci < Pstream::nProcs(); proci++)
    {
        label oldOffset = offset;
        offset += localSizes[proci];

        if (offset < oldOffset)
        {
            FatalErrorInFunction
                << "Overflow : sum of sizes " << localSizes
                << " exceeds capability of label (" << labelMax
                << "). Please recompile with larger datatype for label."
                << exit(FatalError);
        }
        offsets_[proci+1] = offset;
    }
}


Foam::globalIndex::globalIndex(const labelList& offsets)
:
    offsets_(offsets)
{}


Foam::globalIndex::globalIndex(Istream& is)
{
    is >> offsets_;
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, globalIndex& gi)
{
    return is >> gi.offsets_;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const globalIndex& gi)
{
    return os << gi.offsets_;
}


// ************************************************************************* //
