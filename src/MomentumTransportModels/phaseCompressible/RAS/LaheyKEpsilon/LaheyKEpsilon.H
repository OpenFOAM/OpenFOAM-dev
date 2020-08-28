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
    Foam::RASModels::LaheyKEpsilon

Description
    Continuous-phase k-epsilon model including bubble-generated turbulence.

    Reference:
    \verbatim
        Lahey Jr, R. T. (2005).
        The simulation of multidimensional multiphase flows.
        Nuclear Engineering and Design, 235(10), 1043-1060.
    \endverbatim

    The default model coefficients are
    \verbatim
        LaheyKEpsilonCoeffs
        {
            Cmu             0.09;
            C1              1.44;
            C2              1.92;
            C3              -0.33;
            sigmak          1.0;
            sigmaEps        1.3;
            Cp              0.25;
            Cmub            0.6;
            alphaInversion  0.3;
        }
    \endverbatim

SourceFiles
    LaheyKEpsilon.C

\*---------------------------------------------------------------------------*/

#ifndef LaheyKEpsilon_H
#define LaheyKEpsilon_H

#include "kEpsilon.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                        Class LaheyKEpsilon Declaration
\*---------------------------------------------------------------------------*/

template<class BasicMomentumTransportModel>
class LaheyKEpsilon
:
    public kEpsilon<BasicMomentumTransportModel>
{
    // Private Data

        mutable const PhaseCompressibleMomentumTransportModel
        <
            typename BasicMomentumTransportModel::transportModel
        > *gasTurbulencePtr_;


    // Private Member Functions

        //- Return the turbulence model for the gas phase
        const PhaseCompressibleMomentumTransportModel
        <
            typename BasicMomentumTransportModel::transportModel
        >&
        gasTurbulence() const;


protected:

    // Protected data

        // Model coefficients

            dimensionedScalar alphaInversion_;
            dimensionedScalar Cp_;
            dimensionedScalar C3_;
            dimensionedScalar Cmub_;


    // Protected Member Functions

        virtual void correctNut();
        tmp<volScalarField> bubbleG() const;
        tmp<volScalarField> phaseTransferCoeff() const;
        virtual tmp<fvScalarMatrix> kSource() const;
        virtual tmp<fvScalarMatrix> epsilonSource() const;


public:

    typedef typename BasicMomentumTransportModel::alphaField alphaField;
    typedef typename BasicMomentumTransportModel::rhoField rhoField;
    typedef typename BasicMomentumTransportModel::transportModel transportModel;


    //- Runtime type information
    TypeName("LaheyKEpsilon");


    // Constructors

        //- Construct from components
        LaheyKEpsilon
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
        LaheyKEpsilon(const LaheyKEpsilon&) = delete;


    //- Destructor
    virtual ~LaheyKEpsilon()
    {}


    // Member Functions

        //- Read model coefficients if they have changed
        virtual bool read();

        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const LaheyKEpsilon&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "LaheyKEpsilon.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
