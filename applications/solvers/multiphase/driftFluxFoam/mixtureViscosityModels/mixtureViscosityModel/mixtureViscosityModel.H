/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2022 OpenFOAM Foundation
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

Namespace
    Foam::mixtureViscosityModels

Description
    A namespace for incompressible mixtureViscosityModel implementations.

Class
    Foam::mixtureViscosityModel

Description
    An abstract base class for incompressible mixtureViscosityModels.

SourceFiles
    mixtureViscosityModel.C
    mixtureViscosityModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef mixtureViscosityModel_H
#define mixtureViscosityModel_H

#include "incompressibleTwoPhaseInteractingMixture.H"
#include "viscosityModel.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class mixtureViscosityModel Declaration
\*---------------------------------------------------------------------------*/

class mixtureViscosityModel
:
    public viscosityModel
{

protected:

    // Protected data

        //- Mixture properties
        const incompressibleTwoPhaseInteractingMixture& mixture_;


public:

    //- Runtime type information
    TypeName("mixtureViscosityModel");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            mixtureViscosityModel,
            dictionary,
            (
                const incompressibleTwoPhaseInteractingMixture& mixture
            ),
            (mixture)
        );


    // Constructors

        //- Construct from mixture
        mixtureViscosityModel
        (
            const incompressibleTwoPhaseInteractingMixture& mixture
        );

        //- Disallow default bitwise copy construction
        mixtureViscosityModel(const mixtureViscosityModel&) = delete;


    // Selectors

        //- Return a reference to the selected viscosity model
        static autoPtr<mixtureViscosityModel> New
        (
            const incompressibleTwoPhaseInteractingMixture& mixture
        );


    //- Destructor
    virtual ~mixtureViscosityModel()
    {}


    // Member Functions

        //- Return the laminar viscosity
        virtual tmp<volScalarField> nu() const
        {
            NotImplemented;
            return tmp<volScalarField>(nullptr);
        }

        //- Return the laminar viscosity for patch
        virtual tmp<scalarField> nu(const label patchi) const
        {
            NotImplemented;
            return tmp<scalarField>(nullptr);
        }

        //- Correct the laminar viscosity
        virtual void correct()
        {}

        //- Return the mixture viscosity
        //  given the viscosity and velocity of the continuous phase
        virtual tmp<volScalarField> mu
        (
            const volScalarField& muc,
            const volVectorField& U
        ) const = 0;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const mixtureViscosityModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
