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
    Foam::incompressible::RASModels::LaunderGibsonRSTM

Group
    grpIcoRASTurbulence

Description
    Launder-Gibson Reynolds stress turbulence model for incompressible flows.

    The default model coefficients correspond to the following:
    \verbatim
        LaunderGibsonRSTMCoeffs
        {
            Cmu         0.09;
            kappa       0.41;
            Clg1        1.8;
            Clg2        0.6;
            C1          1.44;
            C2          1.92;
            C1Ref       0.5;
            C2Ref       0.3;
            Cs          0.25;
            Ceps        0.15;
            sigmaEps    1.3;
            sigmaR      0.81967;
            couplingFactor  0.0;
        }
    \endverbatim

SourceFiles
    LaunderGibsonRSTM.C

\*---------------------------------------------------------------------------*/

#ifndef LaunderGibsonRSTM_H
#define LaunderGibsonRSTM_H

#include "RASModel.H"
#include "wallDistReflection.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                           Class LaunderGibsonRSTM Declaration
\*---------------------------------------------------------------------------*/

class LaunderGibsonRSTM
:
    public RASModel
{

protected:

    // Protected data

        // Model coefficients

            dimensionedScalar Cmu_;
            dimensionedScalar kappa_;

            dimensionedScalar Clg1_;
            dimensionedScalar Clg2_;

            dimensionedScalar C1_;
            dimensionedScalar C2_;
            dimensionedScalar Cs_;
            dimensionedScalar Ceps_;
            dimensionedScalar sigmaR_;
            dimensionedScalar sigmaEps_;

            dimensionedScalar C1Ref_;
            dimensionedScalar C2Ref_;

            dimensionedScalar couplingFactor_;


        // Fields

            wallDistReflection yr_;

            volSymmTensorField R_;
            volScalarField k_;
            volScalarField epsilon_;
            volScalarField nut_;


public:

    //- Runtime type information
    TypeName("LaunderGibsonRSTM");

    // Constructors

        //- Construct from components
        LaunderGibsonRSTM
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            transportModel& transport,
            const word& turbulenceModelName = turbulenceModel::typeName,
            const word& modelName = typeName
        );


    //- Destructor
    virtual ~LaunderGibsonRSTM()
    {}


    // Member Functions

        //- Return the turbulence viscosity
        virtual tmp<volScalarField> nut() const
        {
            return nut_;
        }

        //- Return the effective diffusivity for R
        tmp<volScalarField> DREff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField("DREff", nut_/sigmaR_ + nu())
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
        virtual tmp<volSymmTensorField> R() const
        {
            return R_;
        }

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
