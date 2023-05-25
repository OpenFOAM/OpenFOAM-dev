/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2023 OpenFOAM Foundation
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
    Foam::wallDampingModel

Description
    Wall damping models can be used to filter interfacial models near the walls.
    This is particularly useful for the lift force because of its dependence on
    the velocity gradient.

    All damping functions accept the following parameters:

      - Cd: A coefficient for filtering the distance from the wall based on the
        dispersed phase diameter. This can be useful to correct gradient
        sampling error when the dispersed phase diameter is significantly
        larger than near wall mesh resolution.
      - zeroWallDist: A constant offset from the wall for the zero point of
        the damping function. Below this distance, the damping will reduce the
        value to zero.
      - zeroInNearWallCells: A switch which sets the value to zero in near wall
        cells regardless of the other parameters. This is recommended to be set
        if a lift force is applied together with turbulent wall functions.

Usage
    \table
        Property     | Description                 | Required | Default value
        Cd           | Diameter coefficient        | yes      | none
        zeroWallDist | Offset from wall            | no       | 0
        zeroInNearWallCells | Zero near wall cells | no       | no
    \endtable

SourceFiles
    wallDampingModel.C
    wallDampingModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef wallDampingModel_H
#define wallDampingModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "wallDependentModel.H"
#include "volFields.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"
#include "dispersedPhaseInterface.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class wallDampingModel Declaration
\*---------------------------------------------------------------------------*/

class wallDampingModel
:
    public wallDependentModel
{
protected:

    // Protected data

        //- Interface
        const dispersedPhaseInterface interface_;

        //- Diameter coefficient
        const dimensionedScalar Cd_;

        //- Distance from wall below which the field is damped to zero
        const dimensionedScalar zeroWallDist_;

        //- Set the value to zero in wall-adjacent cells
        const Switch zeroInNearWallCells_;


    // Protected member functions

        //- Return the force limiter field
        virtual tmp<volScalarField> limiter() const = 0;


public:

    //- Runtime type information
    TypeName("wallDampingModel");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            wallDampingModel,
            dictionary,
            (
                const dictionary& dict,
                const phaseInterface& interface
            ),
            (dict, interface)
        );


    // Constructors

        //- Construct from a dictionary and an interface
        wallDampingModel
        (
            const dictionary& dict,
            const phaseInterface& interface
        );


    //- Destructor
    virtual ~wallDampingModel();


    // Selectors

        static autoPtr<wallDampingModel> New
        (
            const dictionary& dict,
            const phaseInterface& interface
        );


    // Member Functions

        //- Return damped coefficient
        virtual tmp<volScalarField> damping() const;

        //- Return damped face coefficient
        virtual tmp<surfaceScalarField> dampingf() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
