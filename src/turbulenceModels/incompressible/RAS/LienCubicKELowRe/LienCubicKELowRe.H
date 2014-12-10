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
    Foam::incompressible::RASModels::LienCubicKELowRe

Group
    grpIcoRASTurbulence

Description
    Lien cubic non-linear low-Reynolds k-epsilon turbulence models for
    incompressible flows.

    References:
    \verbatim
        Lien, F.S., Chen, W.L., Leschziner, M.A.,
        "Low-Reynolds-number eddy-viscosity modeling based on non-linear
         stress-strain/vorticity relations"
        Engineering Turbulence Modelling and Experiments 3
        (Edited by Rodi, W. and Bergeles, G.), 91-100. 1996.
        Elsevier Science Publishers.

        Etemad, S., et al.,
        "Turbulent flow and heat transfer in a square-sectioned U bend"
        Progress in compuational fluid dynamics 6, 89-100. 2006.
    \endverbatim

SourceFiles
    LienCubicKELowRe.C

\*---------------------------------------------------------------------------*/

#ifndef LienCubicKELowRe_H
#define LienCubicKELowRe_H

#include "RASModel.H"
#include "wallDist.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                           Class LienCubicKELowRe Declaration
\*---------------------------------------------------------------------------*/

class LienCubicKELowRe
:
    public RASModel
{

protected:

    // Protected data

        // Model coefficients

            dimensionedScalar C1_;
            dimensionedScalar C2_;
            dimensionedScalar sigmak_;
            dimensionedScalar sigmaEps_;
            dimensionedScalar A1_;
            dimensionedScalar A2_;
            dimensionedScalar Ctau1_;
            dimensionedScalar Ctau2_;
            dimensionedScalar Ctau3_;
            dimensionedScalar alphaKsi_;

            dimensionedScalar CmuWall_;
            dimensionedScalar kappa_;

            dimensionedScalar Am_;
            dimensionedScalar Aepsilon_;
            dimensionedScalar Amu_;


        // Fields

            volScalarField k_;
            volScalarField epsilon_;

            wallDist y_;

            volScalarField eta_;
            volScalarField ksi_;
            volScalarField Cmu_;
            volScalarField fEta_;
            volScalarField C5viscosity_;

            volScalarField yStar_;

            volScalarField nut_;

            volSymmTensorField nonlinearStress_;


public:

    //- Runtime type information
    TypeName("LienCubicKELowRe");

    // Constructors

        //- Construct from components
        LienCubicKELowRe
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            transportModel& transport,
            const word& turbulenceModelName = turbulenceModel::typeName,
            const word& modelName = typeName
        );


    //- Destructor
    virtual ~LienCubicKELowRe()
    {}


    // Member Functions

        //- Return the turbulence viscosity
        virtual tmp<volScalarField> nut() const
        {
            return nut_;
        }

        //- Return the effective diffusivity for k
        tmp<volScalarField> DkEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField("DkEff", nut_/sigmak_ + nu())
            );
        }

        //- Return the effective diffusivity for epsilon
        tmp<volScalarField> DepsilonEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField("DepsilonEff", nut_/sigmaEps_ + nu())
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

        //- Return the Reynolds stress tensor
        virtual tmp<volSymmTensorField> R() const;

        //- Return the effective stress tensor including the laminar stress
        virtual tmp<volSymmTensorField> devReff() const;

        //- Return the source term for the momentum equation
        virtual tmp<fvVectorMatrix> divDevReff(volVectorField& U) const;

        //- Return the source term for the momentum equation
        virtual tmp<fvVectorMatrix> divDevRhoReff
        (
            const volScalarField& rho,
            volVectorField& U
        ) const;

        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct();

        //- Read RASProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // Edn namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
