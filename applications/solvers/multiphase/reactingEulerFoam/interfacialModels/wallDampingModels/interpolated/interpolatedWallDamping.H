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
    Foam::wallDampingModels::interpolated

Description

SourceFiles
    interpolatedWallDamping.C

\*---------------------------------------------------------------------------*/

#ifndef interpolatedWallDamping_H
#define interpolatedWallDamping_H

#include "wallDampingModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class phasePair;

namespace wallDampingModels
{

/*---------------------------------------------------------------------------*\
                           Class interpolated Declaration
\*---------------------------------------------------------------------------*/

class interpolated
:
    public wallDampingModel
{
    // Private Data

        //- Set limiter to zero in near wall cells
        const Switch zeroInNearWallCells_;

protected:

    // Protected member functions

        //- Return the force limiter field
        virtual tmp<volScalarField> limiter() const = 0;

        //- Set limiter to zero in near wall cells
        virtual tmp<volScalarField> zeroNearWallCells
        (
            const tmp<volScalarField>&
        ) const;


public:

    //- Runtime type information
    TypeName("interpolated");


    // Constructors

        //- Construct from components
        interpolated
        (
            const dictionary& dict,
            const phasePair& pair
        );


    //- Destructor
    virtual ~interpolated();


    // Member Functions

        //- Return damped coefficient
        virtual tmp<volScalarField> damp
        (
            const tmp<volScalarField>&
        ) const;

        //- Return damped force
        virtual tmp<volVectorField> damp
        (
            const tmp<volVectorField>&
        ) const;

        //- Return damped face force
        virtual tmp<surfaceScalarField> damp
        (
            const tmp<surfaceScalarField>&
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace wallDampingModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
