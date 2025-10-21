/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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
    Foam::RASModels::ubMomentumTransportModel

Description

SourceFiles
    ubMomentumTransportModel.C

\*---------------------------------------------------------------------------*/

#ifndef ubMomentumTransportModel_H
#define ubMomentumTransportModel_H

#include "phaseCompressibleMomentumTransportModel.H"
#include "RASModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                       Class ubMomentumTransportModel Declaration
\*---------------------------------------------------------------------------*/

class ubMomentumTransportModel
:
    public RASModel<phaseCompressibleMomentumTransportModel>
{
    // Private Data

        const compressible::momentumTransportModel& mixtureMomentumTransport_;


public:

    //- Runtime type information
    TypeName("ubMomentumTransportModel");


    // Constructors

        //- Construct from components
        ubMomentumTransportModel
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const viscosity& viscosity,
            const compressible::momentumTransportModel& mixtureMomentumTransport
        );

        //- Disallow default bitwise copy construction
        ubMomentumTransportModel(const ubMomentumTransportModel&) = delete;


    //- Destructor
    virtual ~ubMomentumTransportModel()
    {}


    // Member Functions

        //- Re-read model coefficients if they have changed
        virtual bool read();

        //- Return the turbulence kinetic energy
        virtual tmp<volScalarField> k() const
        {
            return mixtureMomentumTransport_.k();
        }

        //- Return the turbulence kinetic energy dissipation rate
        virtual tmp<volScalarField> epsilon() const
        {
            return mixtureMomentumTransport_.epsilon();
        }

        //- Return the turbulence specific dissipation rate
        virtual tmp<volScalarField> omega() const
        {
            return mixtureMomentumTransport_.omega();
        }

        //- Return the turbulence viscosity
        virtual tmp<volScalarField> nut() const
        {
            return mixtureMomentumTransport_.nut();
        }

        //- Return the turbulence viscosity on patch
        virtual tmp<scalarField> nut(const label patchi) const
        {
            return mixtureMomentumTransport_.nut(patchi);
        }

        //- Return the Reynolds stress tensor [m^2/s^2]
        virtual tmp<volSymmTensorField> R() const
        {
            return mixtureMomentumTransport_.R();
        }

        //- Return the effective surface stress
        virtual tmp<surfaceVectorField> devTau() const
        {
            NotImplemented;
            return tmp<surfaceVectorField>(nullptr);
        }

        //- Return the source term for the momentum equation
        virtual tmp<fvVectorMatrix> divDevTau(volVectorField& U) const;

        //- Return the source term for the momentum equation
        virtual tmp<fvVectorMatrix> divDevTau
        (
            const volScalarField& rho,
            volVectorField& U
        ) const;

        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const ubMomentumTransportModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
