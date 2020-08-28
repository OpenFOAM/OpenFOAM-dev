/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2020 OpenFOAM Foundation
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
    Foam::liftModels::wallDamped

Description

SourceFiles
    wallDamped.C

\*---------------------------------------------------------------------------*/

#ifndef wallDampedLift_H
#define wallDampedLift_H

#include "liftModel.H"
#include "wallDampingModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class phasePair;

namespace liftModels
{

/*---------------------------------------------------------------------------*\
                         Class wallDamped Declaration
\*---------------------------------------------------------------------------*/

class wallDamped
:
    public liftModel
{
    // Private Data

        //- The lift model to damp
        autoPtr<liftModel> liftModel_;

        //- The wall-damping model
        autoPtr<wallDampingModel> wallDampingModel_;


public:

    //- Runtime type information
    TypeName("wallDamped");


    // Constructors

        //- Construct from a dictionary and a phase pair
        wallDamped
        (
            const dictionary& dict,
            const phasePair& pair
        );


    //- Destructor
    virtual ~wallDamped();


    // Member Functions

        //- Return lift coefficient
        virtual tmp<volScalarField> Cl() const;

        //- Return phase-intensive lift force
        virtual tmp<volVectorField> Fi() const;

        //- Return lift force
        virtual tmp<volVectorField> F() const;

        //- Return face lift force
        virtual tmp<surfaceScalarField> Ff() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace liftModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
