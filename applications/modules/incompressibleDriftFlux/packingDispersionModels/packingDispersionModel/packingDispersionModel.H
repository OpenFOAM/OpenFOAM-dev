/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024 OpenFOAM Foundation
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
    Foam::packingDispersionModel

Description
    Packing dispersion model

SourceFiles
    packingDispersionModel.C

\*---------------------------------------------------------------------------*/

#ifndef packingDispersionModel_H
#define packingDispersionModel_H

#include "relativeVelocityModel.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class packingDispersionModel Declaration
\*---------------------------------------------------------------------------*/

class packingDispersionModel
{

protected:

    // Protected data

        //- Relative velocity model
        const relativeVelocityModel& relativeVelocity_;

        //- Mixture properties
        const incompressibleDriftFluxMixture& mixture_;


public:

    //- Runtime type information
    TypeName("packingDispersionModel");

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        packingDispersionModel,
        dictionary,
        (
            const dictionary& dict,
            const relativeVelocityModel& relativeVelocity
        ),
        (dict, relativeVelocity)
    );


    // Constructors

        //- Construct from components
        packingDispersionModel
        (
            const relativeVelocityModel& relativeVelocity
        );

        //- Disallow default bitwise copy construction
        packingDispersionModel(const packingDispersionModel&) = delete;


    // Selector
    static autoPtr<packingDispersionModel> New
    (
        const dictionary& dict,
        const relativeVelocityModel& relativeVelocity
    );


    //- Destructor
    virtual ~packingDispersionModel();


    // Member Functions

        //- Return the derivative of the packing stress w.r.t. phase-fraction
        virtual tmp<volScalarField> sigmaPrime() const = 0;

        //- Return the diffusion coefficient of the dispersed phase
        virtual tmp<volScalarField> Dd() const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const packingDispersionModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
