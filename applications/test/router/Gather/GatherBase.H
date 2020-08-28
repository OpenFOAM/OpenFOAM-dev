/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

Class
    Foam::AddOp

Description

SourceFiles
    GatherBase.C

\*---------------------------------------------------------------------------*/

#ifndef GatherBase_H
#define GatherBase_H

#include "List.H"
#include "labelList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class GatherBase Declaration
\*---------------------------------------------------------------------------*/

class GatherBase
{

public:

        //- Flatten: appends all elements of list into one single list.
        //  Used to collapse 'Gathered' data.
        template<class T>
        static T flatten(const List<T>);

        //- Flatten and offset 'Gathered' indices (into value) so they
        //  remain valid with respect to values (after they have been flattened)
        template<class DataType, class IndexType, class AddOp>
        static IndexType offset
        (
            const List<DataType>& values,
            const List<IndexType>& indices,
            AddOp aop
        );
};


template<class T>
class AddOp
{

public:

    T operator()
    (
        const T& x,
        const label offset
    ) const
    {
        return x + offset;
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "GatherBase.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
