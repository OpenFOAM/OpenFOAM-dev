/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2020 OpenFOAM Foundation
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
    Foam::RASModels::mixtureKEpsilon

Description
    Mixture k-epsilon turbulence model for two-phase gas-liquid systems

    The basic structure of the model is based on:
    \verbatim
        Behzadi, A., Issa, R. I., & Rusche, H. (2004).
        Modelling of dispersed bubble and droplet flow at high phase fractions.
        Chemical Engineering Science, 59(4), 759-770.
    \endverbatim

    but an effective density for the gas is used in the averaging and an
    alternative model for bubble-generated turbulence from:
    \verbatim
        Lahey Jr, R. T. (2005).
        The simulation of multidimensional multiphase flows.
        Nuclear Engineering and Design, 235(10), 1043-1060.
    \endverbatim

    The default model coefficients are
    \verbatim
        mixtureKEpsilonCoeffs
        {
            Cmu         0.09;
            C1          1.44;
            C2          1.92;
            C3          C2;
            Cp          0.25; // Bubble generated turbulence
            alphap      0.3;  // Gas phase fraction below which
                              // bubble generated turbulence is included
            sigmak      1.0;
            sigmaEps    1.3;
        }
    \endverbatim

SourceFiles
    mixtureKEpsilon.C

\*---------------------------------------------------------------------------*/

#ifndef mixtureKEpsilon_H
#define mixtureKEpsilon_H

#include "RASModel.H"
#include "eddyViscosity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                       Class mixtureKEpsilon Declaration
\*---------------------------------------------------------------------------*/

template<class BasicMomentumTransportModel>
class mixtureKEpsilon
:
    public eddyViscosity<RASModel<BasicMomentumTransportModel>>
{
    // Private Data

        mutable mixtureKEpsilon<BasicMomentumTransportModel>
            *liquidTurbulencePtr_;


    // Private Member Functions

        //- Return the turbulence model for the other phase
        mixtureKEpsilon<BasicMomentumTransportModel>& liquidTurbulence() const;


protected:

    // Protected data

        // Model coefficients

            dimensionedScalar Cmu_;
            dimensionedScalar C1_;
            dimensionedScalar C2_;
            dimensionedScalar C3_;
            dimensionedScalar Cp_;
            dimensionedScalar alphap_;
            dimensionedScalar sigmak_;
            dimensionedScalar sigmaEps_;

        // Fields

            volScalarField k_;
            volScalarField epsilon_;

        // Mixture fields

            autoPtr<volScalarField> Ct2_;
            autoPtr<volScalarField> rhom_;
            autoPtr<volScalarField> km_;
            autoPtr<volScalarField> epsilonm_;


    // Protected Member Functions

        wordList epsilonBoundaryTypes(const volScalarField& epsilon) const;

        void correctInletOutlet
        (
            volScalarField& vsf,
            const volScalarField& refVsf
        ) const;

        void initMixtureFields();

        virtual void correctNut();

        tmp<volScalarField> Ct2() const;

        tmp<volScalarField> rholEff() const;
        tmp<volScalarField> rhogEff() const;
        tmp<volScalarField> rhom() const;

        tmp<volScalarField> mix
        (
            const volScalarField& fc,
            const volScalarField& fd
        ) const;

        tmp<volScalarField> mixU
        (
            const volScalarField& fc,
            const volScalarField& fd
        ) const;

        tmp<surfaceScalarField> mixFlux
        (
            const surfaceScalarField& fc,
            const surfaceScalarField& fd
        ) const;

        tmp<volScalarField> bubbleG() const;
        virtual tmp<fvScalarMatrix> kSource() const;
        virtual tmp<fvScalarMatrix> epsilonSource() const;

        //- Return the effective diffusivity for k
        tmp<volScalarField> DkEff(const volScalarField& nutm) const
        {
            return volScalarField::New
            (
                "DkEff",
                nutm/sigmak_
            );
        }

        //- Return the effective diffusivity for epsilon
        tmp<volScalarField> DepsilonEff(const volScalarField& nutm) const
        {
            return volScalarField::New
            (
                "DepsilonEff",
                nutm/sigmaEps_
            );
        }

public:

    typedef typename BasicMomentumTransportModel::alphaField alphaField;
    typedef typename BasicMomentumTransportModel::rhoField rhoField;
    typedef typename BasicMomentumTransportModel::transportModel transportModel;


    //- Runtime type information
    TypeName("mixtureKEpsilon");


    // Constructors

        //- Construct from components
        mixtureKEpsilon
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& type = typeName
        );

        //- Disallow default bitwise copy construction
        mixtureKEpsilon(const mixtureKEpsilon&) = delete;


    //- Destructor
    virtual ~mixtureKEpsilon()
    {}


    // Member Functions

        //- Re-read model coefficients if they have changed
        virtual bool read();

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


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const mixtureKEpsilon&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "mixtureKEpsilon.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
