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

InNamespace
    Foam

Description
    Functions to operate on lists of boolean values

\*---------------------------------------------------------------------------*/

#ifndef boolListOps_H
#define boolListOps_H

#include "PstreamReduceOps.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class BoolList>
inline bool any(const BoolList& l)
{
    forAll(l, i)
    {
        if (l[i])
        {
            return true;
        }
    }
    return false;
}

template<class BoolList>
inline bool gAny(const BoolList& l, const label comm = UPstream::worldComm)
{
    return returnReduce(any(l), orOp(), Pstream::msgType(), comm);
}

template<class BoolList>
inline bool all(const BoolList& l)
{
    forAll(l, i)
    {
        if (!l[i])
        {
            return false;
        }
    }
    return true;
}

template<class BoolList>
inline bool gAll(const BoolList& l, const label comm = UPstream::worldComm)
{
    return returnReduce(all(l), andOp(), Pstream::msgType(), comm);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
