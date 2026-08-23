/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2026 OpenFOAM Foundation
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
    Inter-processor communication reduction functions.

\*---------------------------------------------------------------------------*/

#ifndef PstreamReduceOps_H
#define PstreamReduceOps_H

#include "Pstream.H"
#include "ops.H"
#include "vector2D.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Reduce operation with user specified communication schedule
template<class T, class BinaryOp>
void reduce
(
    const List<UPstream::commsStruct>& comms,
    T& Value,
    const BinaryOp& bop,
    const int tag,
    const label comm
)
{
    if (UPstream::warnComm != -1 && comm != UPstream::warnComm)
    {
        Pout<< "** reducing:" << Value << " with comm:" << comm
            << endl;
        error::printStack(Pout);
    }
    Pstream::gather(comms, Value, bop, tag, comm);
    Pstream::scatter(comms, Value, tag, comm);
}


// Reduce using either linear or tree communication schedule
template<class T, class BinaryOp>
void reduce
(
    T& Value,
    const BinaryOp& bop,
    const int tag = Pstream::msgType(),
    const label comm = UPstream::worldComm
)
{
    if (UPstream::nProcs(comm) < UPstream::nProcsSimpleSum)
    {
        reduce(UPstream::linearCommunication(comm), Value, bop, tag, comm);
    }
    else
    {
        reduce(UPstream::treeCommunication(comm), Value, bop, tag, comm);
    }
}


// Reduce using either linear or tree communication schedule
template<class T, class BinaryOp>
T returnReduce
(
    const T& Value,
    const BinaryOp& bop,
    const int tag = Pstream::msgType(),
    const label comm = UPstream::worldComm
)
{
    T WorkValue(Value);

    if (UPstream::nProcs(comm) < UPstream::nProcsSimpleSum)
    {
        reduce
        (
            UPstream::linearCommunication(comm),
            WorkValue,
            bop,
            tag,
            comm
        );
    }
    else
    {
        reduce
        (
            UPstream::treeCommunication(comm),
            WorkValue,
            bop,
            tag,
            comm
        );
    }

    return WorkValue;
}


// Non-blocking version of reduce. Sets request.
template<class T, class BinaryOp>
void reduce
(
    T& Value,
    const BinaryOp& bop,
    const int tag,
    const label comm,
    label& request
)
{
    NotImplemented;
}


// Insist there are specialisations for the common reductions of scalar(s)
void reduce
(
    scalar& Value,
    const sumOp& bop,
    const int tag = Pstream::msgType(),
    const label comm = UPstream::worldComm
);

void reduce
(
    scalar& Value,
    const minOp& bop,
    const int tag = Pstream::msgType(),
    const label comm = UPstream::worldComm
);

void reduce
(
    vector2D& Value,
    const sumOp& bop,
    const int tag = Pstream::msgType(),
    const label comm = UPstream::worldComm
);

void sumReduce
(
    scalar& Value,
    label& Count,
    const int tag = Pstream::msgType(),
    const label comm = UPstream::worldComm
);

void reduce
(
    scalar& Value,
    const sumOp& bop,
    const int tag,
    const label comm,
    label& request
);


//- Return the largest element of the global list
// or pTraits<ListType::value_type>::min if the list is empty
template<class Type>
Type gMax(const UList<Type>& l, const label comm = UPstream::worldComm)
{
    Type m = max(l);
    reduce(m, maxOp(), Pstream::msgType(), comm);
    return m;
}


//- Return the smallest element of the global list
// or pTraits<ListType::value_type>::max if the list is empty
template<class Type>
Type gMin(const UList<Type>& l, const label comm = UPstream::worldComm)
{
    Type m = min(l);
    reduce(m, minOp(), Pstream::msgType(), comm);
    return m;
}


//- Return the sum of the element of the global list
// or zero if the list is empty
template<class Type>
Type gSum(const UList<Type>& l, const label comm = UPstream::worldComm)
{
    Type s = sum(l);
    reduce(s, sumOp(), Pstream::msgType(), comm);
    return s;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
