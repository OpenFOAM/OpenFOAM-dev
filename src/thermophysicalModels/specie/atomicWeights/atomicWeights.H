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
    Foam::atomicWeightTable

Description
    A table of atomic weights for all the elements

SourceFiles
    atomicWeights.C

\*---------------------------------------------------------------------------*/

#ifndef atomicWeights_H
#define atomicWeights_H

#include "scalar.H"
#include "HashTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class atomicWeights Declaration
\*---------------------------------------------------------------------------*/

class atomicWeightTable
:
    public HashTable<scalar>
{

public:

    // Public types and data

        //- Structure to hold the element name and atomic weight pair
        struct atomicWeight
        {
            char name[3];
            scalar weight;
        };

        static const int nElements = 104;

        //- Static table of the weights of all known elements
        static const atomicWeight atomicWeights[nElements];


    // Constructors

        //- Construct from atomicWeights_
        atomicWeightTable();
};


// * * * * * * * * * * * * * * * * Global data  * * * * * * * * * * * * * * //

// Atomic weights table for every element in the periodic table
extern atomicWeightTable atomicWeights;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
