/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2019 OpenFOAM Foundation
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
    Foam::wallDampingModels::linear

Description

SourceFiles
    linearWallDamping.C

\*---------------------------------------------------------------------------*/

#ifndef linearWallDamping_H
#define linearWallDamping_H

#include "interpolatedWallDamping.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class phasePair;

namespace wallDampingModels
{

/*---------------------------------------------------------------------------*\
                           Class linear Declaration
\*---------------------------------------------------------------------------*/

class linear
:
    public interpolated
{
    // Private Data

        //- Diameter coefficient
        const dimensionedScalar Cd_;


protected:

    // Protected member functions

        //- Return the force limiter field
        virtual tmp<volScalarField> limiter() const;


public:

    //- Runtime type information
    TypeName("linear");


    // Constructors

        //- Construct from components
        linear
        (
            const dictionary& dict,
            const phasePair& pair
        );


    //- Destructor
    virtual ~linear();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace wallDampingModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
