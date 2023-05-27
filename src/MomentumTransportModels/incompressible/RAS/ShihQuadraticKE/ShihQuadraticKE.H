/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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
    Foam::incompressible::RASModels::ShihQuadraticKE

Description
    Shih's quadratic algebraic Reynolds stress k-epsilon turbulence model for
    incompressible flows

    This turbulence model is described in:
    \verbatim
        Shih, T. H., Zhu, J., & Lumley, J. L. (1993).
        A realizable Reynolds stress algebraic equation model.
        NASA technical memorandum 105993.
    \endverbatim

    Implemented according to the specification in:
    <a href=
    "https://personalpages.manchester.ac.uk/staff/david.d.apsley/turbmod.pdf"
    >Apsley: Turbulence Models 2002</a>

SourceFiles
    ShihQuadraticKE.C

\*---------------------------------------------------------------------------*/

#ifndef ShihQuadraticKE_H
#define ShihQuadraticKE_H

#include "incompressibleMomentumTransportModels.H"
#include "nonlinearEddyViscosity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                       Class ShihQuadraticKE Declaration
\*---------------------------------------------------------------------------*/

class ShihQuadraticKE
:
    public nonlinearEddyViscosity<incompressible::RASModel>
{

protected:

    // Protected data

        // Model coefficients

            dimensionedScalar Ceps1_;
            dimensionedScalar Ceps2_;
            dimensionedScalar sigmak_;
            dimensionedScalar sigmaEps_;
            dimensionedScalar Cmu1_;
            dimensionedScalar Cmu2_;
            dimensionedScalar Cbeta_;
            dimensionedScalar Cbeta1_;
            dimensionedScalar Cbeta2_;
            dimensionedScalar Cbeta3_;


        // Fields

            volScalarField k_;
            volScalarField epsilon_;


    // Protected Member Functions

        virtual void correctNut();
        virtual void correctNonlinearStress(const volTensorField& gradU);


public:

    //- Runtime type information
    TypeName("ShihQuadraticKE");


    // Constructors

        //- Construct from components
        ShihQuadraticKE
        (
            const geometricOneField& alpha,
            const geometricOneField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const viscosity& viscosity,
            const word& type = typeName
        );


    //- Destructor
    virtual ~ShihQuadraticKE()
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
                nut_/sigmak_ + nu()
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
                epsilon_/(0.09*k_)
            );
        }

        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
