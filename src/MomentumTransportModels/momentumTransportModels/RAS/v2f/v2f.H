/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2020 OpenFOAM Foundation
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
    Foam::RASModels::v2f

Description
    Lien and Kalitzin's v2-f turbulence model for incompressible and
    compressible flows, with a limit imposed on the turbulent viscosity given
    by Davidson et al.

    The model solves for turbulence kinetic energy k and turbulence dissipation
    rate epsilon, with additional equations for the turbulence stress normal to
    streamlines, v2, and elliptic damping function, f.

    The variant implemented employs N=6, such that f=0 on walls.

    Wall boundary conditions are:

        k       = kLowReWallFunction
        epsilon = epsilonWallFunction
        v2      = v2WallFunction
        f       = fWallFunction

    These are applicable to both low- and high-Reynolds number flows.

    Inlet values can be approximated by:

        v2      = 2/3 k
        f       = zero-gradient

    References:
    \verbatim
        Lien, F. S., & Kalitzin, G. (2001).
        Computations of transonic flow with the v2f turbulence model.
        International Journal of Heat and Fluid Flow, 22(1), 53-61.

        Davidson, L., Nielsen, P., & Sveningsson, A. (2003).
        Modifications of the v2-f model for computing the flow in a
        3D wall jet.
        Turbulence, Heat and Mass Transfer, 4, 577-584
    \endverbatim

    The default model coefficients are
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
            Ceps3       -0.33;
            sigmaEps    1.3;
            sigmaK      1;
        }
    \endverbatim

Note
    If the kLowReWallFunction is employed, a velocity variant of the turbulent
    viscosity wall function should be used, e.g. nutUWallFunction.  Turbulence
    k variants (nutk...) for this case will not behave correctly.

See also
    Foam::RASModels::v2fBase
    Foam::RASModels::kEpsilon
    Foam::kLowReWallFunctionFvPatchScalarField
    Foam::epsilonWallFunctionFvPatchScalarField
    Foam::v2WallFunctionFvPatchScalarField
    Foam::fWallFunctionFvPatchScalarField

SourceFiles
    v2f.C

\*---------------------------------------------------------------------------*/

#ifndef v2f_H
#define v2f_H

#include "v2fBase.H"
#include "RASModel.H"
#include "eddyViscosity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                             Class v2f Declaration
\*---------------------------------------------------------------------------*/

template<class BasicMomentumTransportModel>
class v2f
:
    public eddyViscosity<RASModel<BasicMomentumTransportModel>>,
    public v2fBase
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
            dimensionedScalar Ceps3_;
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


        // Bounding values

            dimensionedScalar v2Min_;
            dimensionedScalar fMin_;


    // Protected Member Functions

        virtual void correctNut();

        //- Return time scale, Ts
        tmp<volScalarField> Ts() const;

        //- Return length scale, Ls
        tmp<volScalarField> Ls() const;


public:

    typedef typename BasicMomentumTransportModel::alphaField alphaField;
    typedef typename BasicMomentumTransportModel::rhoField rhoField;
    typedef typename BasicMomentumTransportModel::transportModel transportModel;


    //- Runtime type information
    TypeName("v2f");


    // Constructors

        //- Construct from components
        v2f
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& type = typeName
        );


    //- Destructor
    virtual ~v2f()
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
                this->nut_/sigmaK_ + this->nu()
            );
        }

        //- Return the effective diffusivity for epsilon
        tmp<volScalarField> DepsilonEff() const
        {
            return volScalarField::New
            (
                "DepsilonEff",
                this->nut_/sigmaEps_ + this->nu()
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

        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "v2f.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
