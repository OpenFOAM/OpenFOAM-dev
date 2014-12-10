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
    Foam::compressible::RASModels::realizableKE

Group
    grpCmpRASTurbulence

Description
    Realizable k-epsilon turbulence model for compressible flows.

    Model described in the paper:
    \verbatim
        "A New k-epsilon Eddy Viscosity Model for High Reynolds Number
        Turbulent Flows"

        Tsan-Hsing Shih, William W. Liou, Aamir Shabbir, Zhigang Tang and
        Jiang Zhu

        Computers and Fluids Vol. 24, No. 3, pp. 227-238, 1995
    \endverbatim

    The default model coefficients correspond to the following:
    \verbatim
        realizableKECoeffs
        {
            Cmu         0.09;
            A0          4.0;
            C2          1.9;
            sigmak      1.0;
            sigmaEps    1.2;
            Prt         1.0;    // only for compressible
        }
    \endverbatim

SourceFiles
    realizableKE.C

\*---------------------------------------------------------------------------*/

#ifndef realizableKE_H
#define realizableKE_H

#include "RASModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                        Class realizableKE Declaration
\*---------------------------------------------------------------------------*/

class realizableKE
:
    public RASModel
{

protected:

    // Protected data

        // Model coefficients

            dimensionedScalar Cmu_;
            dimensionedScalar A0_;
            dimensionedScalar C2_;
            dimensionedScalar sigmak_;
            dimensionedScalar sigmaEps_;
            dimensionedScalar Prt_;


        // Fields

            volScalarField k_;
            volScalarField epsilon_;
            volScalarField mut_;
            volScalarField alphat_;


   // Protected Member Functions

        tmp<volScalarField> rCmu
        (
            const volTensorField& gradU,
            const volScalarField& S2,
            const volScalarField& magS
        );

        tmp<volScalarField> rCmu(const volTensorField& gradU);


public:

    //- Runtime type information
    TypeName("realizableKE");

    // Constructors

        //- Construct from components
        realizableKE
        (
            const volScalarField& rho,
            const volVectorField& U,
            const surfaceScalarField& phi,
            const fluidThermo& thermophysicalModel,
            const word& turbulenceModelName = turbulenceModel::typeName,
            const word& modelName = typeName
        );


    //- Destructor
    virtual ~realizableKE()
    {}


    // Member Functions

        //- Return the effective diffusivity for k
        tmp<volScalarField> DkEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField("DkEff", mut_/sigmak_ + mu())
            );
        }

        //- Return the effective diffusivity for epsilon
        tmp<volScalarField> DepsilonEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField("DepsilonEff", mut_/sigmaEps_ + mu())
            );
        }

        //- Return the turbulence viscosity
        virtual tmp<volScalarField> mut() const
        {
            return mut_;
        }

        //- Return the turbulence thermal diffusivity
        virtual tmp<volScalarField> alphat() const
        {
            return alphat_;
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
        virtual tmp<volSymmTensorField> devRhoReff() const;

        //- Return the source term for the momentum equation
        virtual tmp<fvVectorMatrix> divDevRhoReff(volVectorField& U) const;

        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct();

        //- Read RASProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace compressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
