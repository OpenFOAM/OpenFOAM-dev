/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2013 OpenFOAM Foundation
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
    Foam::incompressible::RASModels::v2f

Group
    grpIcoRASTurbulence

Description
    Lien and Kalitzin's v2-f turbulence model for incompressible flows, with
    a limit imposed on the turbulent viscosity given by Davidson et al.

    The model solves for turbulence k and epsilon, with additional equations
    for the turbulence stress normal to streamlines, v2, and elliptic damping
    function, f.  The variant implemented employs N=6, such that f=0 on walls.

    Wall boundary conditions are:

        k       = kLowReWallFunction
        epsilon = epsilonLowReWallFunction
        v2      = v2WallFunction
        f       = fWallFunction

    These are applicable to both low- and high-Reynolds number flows.

    Inlet values can be approximated by:

        v2      = 2/3 k
        f       = zero-gradient


    References:

    Lien F-S, Kalitzin G, 2001.  Computations of transonic flow with the v2-f
    turbulence model. Int. J. Heat Fluid Flow 22, pp 53-61

    Davidson L, Nielsen P, Sveningsson A, 2003.  Modifications of the v2-f
    model for computing the flow in a 3D wall jet.  Turbulence, Heat and Mass
    Transfer 4, pp 577-584

    The default model coefficients are given as:
    \verbatim
        v2fCoeffs
        {
            Cmu         0.22;
            CmuKEps     0.09;
            C1          1.4;
            C2          0.3;
            CL          0.23;
            Ceta        70;
            Ceps2       1.9;
            sigmaEps    1.3;
            sigmaK      1;
        }
    \endverbatim

Note
    If the kLowReWallFunction is employed, a velocity variant of the turbulent
    viscosity wall function should be used, e.g. nutUWallFunction.  Turbulence
    k variants (nutk...) for this case will not behave correctly.

SeeAlso
    Foam::kEpsilon
    Foam::kLowReWallFunctionFvPatchScalarField
    Foam::epsilonLowReWallFunctionFvPatchScalarField
    Foam::v2WallFunctionFvPatchScalarField
    Foam::fWallFunctionFvPatchScalarField

SourceFiles
    v2f.C

\*---------------------------------------------------------------------------*/

#ifndef v2f_H
#define v2f_H

#include "RASModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                             Class v2f Declaration
\*---------------------------------------------------------------------------*/

class v2f
:
    public RASModel
{

protected:

    // Protected data

        // Model coefficients

            dimensionedScalar Cmu_;
            dimensionedScalar CmuKEps_;
            dimensionedScalar C1_;
            dimensionedScalar C2_;
            dimensionedScalar CL_;
            dimensionedScalar Ceta_;
            dimensionedScalar Ceps2_;
            dimensionedScalar sigmaK_;
            dimensionedScalar sigmaEps_;


        // Fields

            //- Turbulence kinetic energy
            volScalarField k_;

            //- Turbulence dissipation
            volScalarField epsilon_;

            //- Turbulence stress normal to streamlines
            volScalarField v2_;

            //- Damping function
            volScalarField f_;

            //- Turbulence viscosity
            volScalarField nut_;


        // Bounding values

            dimensionedScalar v2Min_;
            dimensionedScalar fMin_;


    // Protected Member Functions

        //- Return boundary type names for the R field
        wordList RBoundaryTypes() const;

        //- Apply Davidson correction to nut
        tmp<volScalarField> davidsonCorrectNut
        (
            const tmp<volScalarField>& value
        ) const;

        //- Return time scale, Ts
        tmp<volScalarField> Ts() const;

        //- Return length scale, Ls
        tmp<volScalarField> Ls() const;


public:

    //- Runtime type information
    TypeName("v2f");

    // Constructors

        //- Construct from components
        v2f
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            transportModel& transport,
            const word& turbulenceModelName = turbulenceModel::typeName,
            const word& modelName = typeName
        );


    //- Destructor
    virtual ~v2f()
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
                new volScalarField("DkEff", nut_/sigmaK_ + nu())
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

        //- Return turbulence stress normal to streamlines
        virtual tmp<volScalarField> v2() const
        {
            return v2_;
        }

        //- Return the damping function
        virtual tmp<volScalarField> f() const
        {
            return f_;
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
} // End namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
