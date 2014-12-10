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

Class
    Foam::energyScalingFunctions::shifted

Description


SourceFiles
    shifted.C

\*---------------------------------------------------------------------------*/

#ifndef shifted_H
#define shifted_H

#include "energyScalingFunction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace energyScalingFunctions
{

/*---------------------------------------------------------------------------*\
                           Class shifted Declaration
\*---------------------------------------------------------------------------*/

class shifted
:
    public energyScalingFunction
{
    // Private data

        scalar e_at_rCut_;


public:

    //- Runtime type information
    TypeName("shifted");


    // Constructors

        //- Construct from components
        shifted
        (
            const word& name,
            const dictionary& energyScalingFunctionProperties,
            const pairPotential& pairPot
        );


    //- Destructor
    ~shifted()
    {}


    // Member Functions

        void scaleEnergy(scalar& e, const scalar r) const;

        //- Read dictionary
        bool read(const dictionary& energyScalingFunctionProperties);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace energyScalingFunctions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
