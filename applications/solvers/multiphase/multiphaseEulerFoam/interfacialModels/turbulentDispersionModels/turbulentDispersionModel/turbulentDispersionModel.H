/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2020 OpenFOAM Foundation
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
    Foam::turbulentDispersionModel

Description

SourceFiles
    turbulentDispersionModel.C
    turbulentDispersionModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef turbulentDispersionModel_H
#define turbulentDispersionModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "volFields.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"
#include "phaseCompressibleMomentumTransportModelFwd.H"

namespace Foam
{

class phasePair;

/*---------------------------------------------------------------------------*\
                  Class turbulentDispersionModel Declaration
\*---------------------------------------------------------------------------*/

class turbulentDispersionModel
{
protected:

    // Protected data

        //- Phase pair
        const phasePair& pair_;


public:

    //- Runtime type information
    TypeName("turbulentDispersionModel");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            turbulentDispersionModel,
            dictionary,
            (
                const dictionary& dict,
                const phasePair& pair
            ),
            (dict, pair)
        );

    // Static Data Members

        //- Diffusivity dimensions
        static const dimensionSet dimD;

        //- Force dimensions
        static const dimensionSet dimF;


    // Constructors

        //- Construct from a dictionary and a phase pair
        turbulentDispersionModel
        (
            const dictionary& dict,
            const phasePair& pair
        );


    //- Destructor
    virtual ~turbulentDispersionModel();


    // Selectors

        static autoPtr<turbulentDispersionModel> New
        (
            const dictionary& dict,
            const phasePair& pair
        );


    // Member Functions

        //- Return a reference to the turbulence model for the continuous phase
        const phaseCompressibleMomentumTransportModel&
            continuousTurbulence() const;

        //- Turbulent diffusivity
        //  multiplying the gradient of the phase-fraction
        virtual tmp<volScalarField> D() const = 0;

        //- Turbulent dispersion force
        virtual tmp<volVectorField> F() const;

        //- Turbulent dispersion force on faces
        virtual tmp<surfaceScalarField> Ff() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
