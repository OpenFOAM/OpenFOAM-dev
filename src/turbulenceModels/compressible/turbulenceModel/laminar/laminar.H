/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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
    Foam::compressible::laminar

Description
    Turbulence model for laminar compressible flow.

SourceFiles
    laminar.C

\*---------------------------------------------------------------------------*/

#ifndef compressibleLaminar_H
#define compressibleLaminar_H

#include "turbulenceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

/*---------------------------------------------------------------------------*\
                           Class laminar Declaration
\*---------------------------------------------------------------------------*/

class laminar
:
    public turbulenceModel
{

public:

    //- Runtime type information
    TypeName("laminar");

    // Constructors

        //- Construct from components
        laminar
        (
            const volScalarField& rho,
            const volVectorField& U,
            const surfaceScalarField& phi,
            const fluidThermo& thermophysicalModel,
            const word& turbulenceModelName = turbulenceModel::typeName
        );


    // Selectors

        //- Return a reference to the selected turbulence model
        static autoPtr<laminar> New
        (
            const volScalarField& rho,
            const volVectorField& U,
            const surfaceScalarField& phi,
            const fluidThermo& thermophysicalModel,
            const word& turbulenceModelName = turbulenceModel::typeName
        );


    //- Destructor
    virtual ~laminar()
    {}


    // Member Functions

        //- Const access to the coefficients dictionary
        virtual const dictionary& coeffDict() const;

        //- Return the turbulence viscosity, i.e. 0 for laminar flow
        virtual tmp<volScalarField> mut() const;

        //- Return the effective viscosity, i.e. the laminar viscosity
        virtual tmp<volScalarField> muEff() const
        {
            return tmp<volScalarField>(new volScalarField("muEff", mu()));
        }

        //- Return the turbulence thermal diffusivity, i.e. 0 for laminar flow
        virtual tmp<volScalarField> alphat() const;

        //- Return the effective turbulent temperature diffusivity,
        //  i.e. the laminar thermal diffusivity
        virtual tmp<volScalarField> kappaEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    "kappaEff",
                    thermo().kappa()
                )
            );
        }

        //- Return the effective turbulent temperature diffusivity for a patch,
        //  i.e. the laminar thermal diffusivity
        virtual tmp<scalarField> kappaEff(const label patchI) const
        {
            return thermo().kappa(patchI);
        }

        //- Return the effective turbulent thermal diffusivity,
        //  i.e. the laminar thermal diffusivity
        virtual tmp<volScalarField> alphaEff() const
        {
            return tmp<volScalarField>(new volScalarField("alphaEff", alpha()));
        }

        //- Return the effective turbulent thermal diffusivity for a patch,
        //  i.e. the laminar thermal diffusivity
        virtual tmp<scalarField> alphaEff(const label patchI) const
        {
            return alpha().boundaryField()[patchI];
        }

        //- Return the turbulence kinetic energy, i.e. 0 for laminar flow
        virtual tmp<volScalarField> k() const;

        //- Return the turbulence kinetic energy dissipation rate,
        //  i.e. 0 for laminar flow
        virtual tmp<volScalarField> epsilon() const;

        //- Return the Reynolds stress tensor, i.e. 0 for laminar flow
        virtual tmp<volSymmTensorField> R() const;

        //- Return the effective stress tensor, i.e. the laminar stress
        virtual tmp<volSymmTensorField> devRhoReff() const;

        //- Return the source term for the momentum equation
        virtual tmp<fvVectorMatrix> divDevRhoReff(volVectorField& U) const;

        //- Correct the laminar viscosity
        virtual void correct();

        //- Read turbulenceProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
