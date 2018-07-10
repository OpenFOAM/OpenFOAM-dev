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

Class
    Foam::EdgeMap

Description
    Map from edge (expressed as its endpoints) to value

\*---------------------------------------------------------------------------*/

#ifndef EdgeMap_H
#define EdgeMap_H

#include "HashTable.H"
#include "edge.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class EdgeMap Declaration
\*---------------------------------------------------------------------------*/

template<class T>
class EdgeMap
:
    public HashTable<T, edge, Hash<edge>>
{

public:

    // Constructors

        //- Construct given initial map size
        EdgeMap(const label size = 128)
        :
            HashTable<T, edge, Hash<edge>>(size)
        {}

        //- Construct from Istream
        EdgeMap(Istream& is)
        :
            HashTable<T, edge, Hash<edge>>(is)
        {}

        //- Construct as copy
        EdgeMap(const EdgeMap<T>& map)
        :
            HashTable<T, edge, Hash<edge>>(map)
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
