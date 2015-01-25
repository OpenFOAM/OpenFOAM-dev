/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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
    Foam::incompressible::RASModels::RNGkEpsilon

Group
    grpIcoRASTurbulence

Description
    Renormalisation group k-epsilon turbulence model for incompressible flows.

    References:
    \verbatim
        "Development of Turbulence Models for Shear Flow by a
         Double Expansion Technique",
        Yakhot, V., Orszag, S.A., Thangam, S., Gatski, T.B. and Speziale, C.G.",
        Physics of Fluids A, Vol. 4, No. 7, 1992, pp 1510-1520.
    \endverbatim

    The default model coefficients are
    \verbatim
        RNGkEpsilonCoeffs
        {
            Cmu         0.0845;
            C1          1.42;
            C2          1.68;
            sigmak      0.71942;
            sigmaEps    0.71942;
            eta0        4.38;
            beta        0.012;
        }
    \endverbatim

SourceFiles
    RNGkEpsilon.C

\*---------------------------------------------------------------------------*/

#ifndef RNGkEpsilon_H
#define RNGkEpsilon_H

#include "turbulentTransportModel.H"
#include "eddyViscosity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                           Class RNGkEpsilon Declaration
\*---------------------------------------------------------------------------*/

class RNGkEpsilon
:
    public eddyViscosity<incompressible::RASModel>
{

protected:

    // Protected data

        // Model coefficients

            dimensionedScalar Cmu_;
            dimensionedScalar C1_;
            dimensionedScalar C2_;
            dimensionedScalar sigmak_;
            dimensionedScalar sigmaEps_;
            dimensionedScalar eta0_;
            dimensionedScalar beta_;


        // Fields

            volScalarField k_;
            volScalarField epsilon_;


    // Protected Member Functions

        virtual void correctNut();


public:

    //- Runtime type information
    TypeName("RNGkEpsilon");

    // Constructors

        //- Construct from components
        RNGkEpsilon
        (
            const geometricOneField& alpha,
            const geometricOneField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& propertiesName = turbulenceModel::propertiesName,
            const word& type = typeName
        );


    //- Destructor
    virtual ~RNGkEpsilon()
    {}


    // Member Functions

        //- Read RASProperties dictionary
        virtual bool read();

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
