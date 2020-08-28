/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
    Foam::incompressible::RASModels::LienLeschziner

Description
    Lien and Leschziner low-Reynolds number k-epsilon turbulence model for
    incompressible flows.

    This turbulence model is described in:
    \verbatim
        Lien, F. S., & Leschziner, M. A. (1993).
        A pressure-velocity solution strategy for compressible flow
        and its application to shock/boundary-layer interaction
        using second-moment turbulence closure.
        Journal of fluids engineering, 115(4), 717-725.
    \endverbatim

    Implemented according to the specification in:
    <a href=
    "https://personalpages.manchester.ac.uk/staff/david.d.apsley/turbmod.pdf"
    >Apsley: Turbulence Models 2002</a>

    In addition to the low-Reynolds number damping functions support for
    wall-functions is also included to allow for low- and high-Reynolds number
    operation.

SourceFiles
    LienLeschziner.C

\*---------------------------------------------------------------------------*/

#ifndef LienLeschziner_H
#define LienLeschziner_H

#include "kinematicMomentumTransportModel.H"
#include "eddyViscosity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                       Class LienLeschziner Declaration
\*---------------------------------------------------------------------------*/

class LienLeschziner
:
    public eddyViscosity<incompressible::RASModel>
{

protected:

    // Protected data

        // Model coefficients

            dimensionedScalar Ceps1_;
            dimensionedScalar Ceps2_;
            dimensionedScalar sigmak_;
            dimensionedScalar sigmaEps_;

            dimensionedScalar Cmu_;
            dimensionedScalar kappa_;

            dimensionedScalar Anu_;
            dimensionedScalar Aeps_;
            dimensionedScalar AE_;


        // Fields

            volScalarField k_;
            volScalarField epsilon_;

            //- Wall distance
            //  Note: different to wall distance in parent RASModel
            //  which is for near-wall cells only
            const volScalarField& y_;


    // Protected Member Functions

        tmp<volScalarField> fMu() const;
        tmp<volScalarField> f2() const;
        tmp<volScalarField> E(const volScalarField& f2) const;

        virtual void correctNut();


public:

    TypeName("LienLeschziner");

    // Constructors

        //- Construct from components
        LienLeschziner
        (
            const geometricOneField& alpha,
            const geometricOneField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& type = typeName
        );


    //- Destructor
    virtual ~LienLeschziner()
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
