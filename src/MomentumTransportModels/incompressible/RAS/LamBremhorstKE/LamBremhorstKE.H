/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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
    Foam::incompressible::RASModels::LamBremhorstKE

Description
    Lam and Bremhorst low-Reynolds number k-epsilon turbulence model
    for incompressible flows

    This turbulence model is described in:
    \verbatim
        Lam, C. K. G., & Bremhorst, K. (1981).
        A modified form of the k-ε model for predicting wall turbulence.
        Journal of Fluids Engineering, 103(3), 456-460.
    \endverbatim

SourceFiles
    LamBremhorstKE.C

\*---------------------------------------------------------------------------*/

#ifndef LamBremhorstKE_H
#define LamBremhorstKE_H

#include "incompressibleMomentumTransportModels.H"
#include "eddyViscosity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                       Class LamBremhorstKE Declaration
\*---------------------------------------------------------------------------*/

class LamBremhorstKE
:
    public eddyViscosity<incompressible::RASModel>
{
    // Private Member Functions

        tmp<volScalarField> Rt() const;
        tmp<volScalarField> fMu(const volScalarField& Rt) const;
        tmp<volScalarField> f1(const volScalarField& fMu) const;
        tmp<volScalarField> f2(const volScalarField& Rt) const;
        void correctNut(const volScalarField& fMu);


protected:

    // Protected data

        dimensionedScalar Cmu_;
        dimensionedScalar Ceps1_;
        dimensionedScalar Ceps2_;
        dimensionedScalar sigmaEps_;

        volScalarField k_;
        volScalarField epsilon_;

        //- Wall distance
        //  Note: different to wall distance in parent RASModel
        //  which is for near-wall cells only
        const volScalarField& y_;


    // Protected Member Functions

        virtual void correctNut();


public:

    //- Runtime type information
    TypeName("LamBremhorstKE");


    // Constructors

        //- Construct from components
        LamBremhorstKE
        (
            const geometricOneField& alpha,
            const geometricOneField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const viscosity& viscosity,
            const word& type = typeName
        );

        //- Disallow default bitwise copy construction
        LamBremhorstKE(const LamBremhorstKE&) = delete;


    //- Destructor
    virtual ~LamBremhorstKE()
    {}


    // Member Functions

        //- Read RASProperties dictionary
        virtual bool read();

        //- Return the effective diffusivity for k
        tmp<volScalarField> DkEff() const
        {
            return volScalarField::New
            (
                "DkEff",
                nut_ + nu()
            );
        }

        //- Return the effective diffusivity for epsilon
        tmp<volScalarField> DepsilonEff() const
        {
            return volScalarField::New
            (
                "DepsilonEff",
                nut_/sigmaEps_ + nu()
            );
        }

        //- Return the turbulence kinetic energy
        virtual tmp<volScalarField> k() const
        {
            return k_;
        }

        //- Return the turbulence kinetic energy dissipation rate
        virtual tmp<volScalarField> epsilon() const
        {
            return epsilon_;
        }

        //- Return the turbulence specific dissipation rate
        virtual tmp<volScalarField> omega() const
        {
            return volScalarField::New
            (
                "omega",
                epsilon_/(Cmu_*k_),
                epsilon_.boundaryField().types()
            );
        }

        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const LamBremhorstKE&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
