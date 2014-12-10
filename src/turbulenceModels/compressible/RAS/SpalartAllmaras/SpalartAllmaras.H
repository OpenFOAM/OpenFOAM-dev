/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
    Foam::compressible::RASModels::SpalartAllmaras

Group
    grpCmpRASTurbulence

Description
    Spalart-Allmaras one-eqn mixing-length model for compressible
    external flows.

    Reference:
    \verbatim
        "A One-Equation Turbulence Model for Aerodynamic Flows"
        P.R. Spalart,
        S.R. Allmaras,
        La Recherche Aerospatiale, No. 1, 1994, pp. 5-21.
    \endverbatim

    Extended according to:
    \verbatim
        "An Unstructured Grid Generation and Adaptive Solution Technique
        for High Reynolds Number Compressible Flows"
        G.A. Ashford,
        Ph.D. thesis, University of Michigan, 1996.
    \endverbatim
    by using the optional flag \c ashfordCorrection

    The default model coefficients correspond to the following:
    \verbatim
        SpalartAllmarasCoeffs
        {
            Cb1         0.1355;
            Cb2         0.622;
            Cw2         0.3;
            Cw3         2.0;
            Cv1         7.1;
            Cv2         5.0;
            sigmaNut    0.66666;
            Prt         1.0;    // only for compressible
            kappa       0.41;
        }
    \endverbatim

SourceFiles
    SpalartAllmaras.C

\*---------------------------------------------------------------------------*/

#ifndef compressibleSpalartAllmaras_H
#define compressibleSpalartAllmaras_H

#include "RASModel.H"
#include "wallDist.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                      Class SpalartAllmaras Declaration
\*---------------------------------------------------------------------------*/

class SpalartAllmaras
:
    public RASModel
{

protected:

    // Protected data

        // Model coefficients

            dimensionedScalar sigmaNut_;
            dimensionedScalar kappa_;
            dimensionedScalar Prt_;

            dimensionedScalar Cb1_;
            dimensionedScalar Cb2_;
            dimensionedScalar Cw1_;
            dimensionedScalar Cw2_;
            dimensionedScalar Cw3_;
            dimensionedScalar Cv1_;
            dimensionedScalar Cv2_;


        //- Optional flag to activate the Ashford correction
        Switch ashfordCorrection_;


        // Fields

            volScalarField nuTilda_;
            volScalarField mut_;
            volScalarField alphat_;


        //- Wall distance
        wallDist d_;


    // Protected Member Functions

        tmp<volScalarField> chi() const;
        tmp<volScalarField> fv1(const volScalarField& chi) const;
        tmp<volScalarField> fv2
        (
            const volScalarField& chi,
            const volScalarField& fv1
        ) const;
        tmp<volScalarField> fv3
        (
            const volScalarField& chi,
            const volScalarField& fv1
        ) const;
        tmp<volScalarField> fw(const volScalarField& Stilda) const;


public:

    //- Runtime type information
    TypeName("SpalartAllmaras");


    // Constructors

        //- Construct from components
        SpalartAllmaras
        (
            const volScalarField& rho,
            const volVectorField& U,
            const surfaceScalarField& phi,
            const fluidThermo& thermophysicalModel,
            const word& turbulenceModelName = turbulenceModel::typeName,
            const word& modelName = typeName
        );


    //- Destructor
    virtual ~SpalartAllmaras()
    {}


    // Member Functions

        //- Return the effective diffusivity for nuTilda
        tmp<volScalarField> DnuTildaEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    "DnuTildaEff",
                    (rho_*nuTilda_ + mu())/sigmaNut_
                )
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
        virtual tmp<volScalarField> k() const;

        //- Return the turbulence kinetic energy dissipation rate
        virtual tmp<volScalarField> epsilon() const;

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
