/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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
    Foam::turbGen

Description
    Generate a turbulent velocity field conforming to a given
    energy spectrum and being divergence free.

SourceFiles
    turbGen.C
    turbGenIO.C

\*---------------------------------------------------------------------------*/

#ifndef turbGen_H
#define turbGen_H

#include "randomGenerator.H"
#include "primitiveFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class Kmesh;

/*---------------------------------------------------------------------------*\
                           Class turbGen Declaration
\*---------------------------------------------------------------------------*/

class turbGen
{
    // Private Data

        const Kmesh& K_;
        const scalar Ea_;
        const scalar k0_;
        randomGenerator rndGen_;


public:

    // Constructors

        //- Construct from components
        turbGen(const Kmesh& k, const scalar EA, const scalar K0);


    // Member Functions

        //- Generate and return a velocity field
        vectorField U();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
