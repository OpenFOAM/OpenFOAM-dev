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
    Foam::LESModels::NicenoKEqn

Description
    One-equation SGS model for the continuous phase in a two-phase system
    including bubble-generated turbulence.

    Reference:
    \verbatim
        Niceno, B., Dhotre, M. T., & Deen, N. G. (2008).
        One-equation sub-grid scale (SGS) modelling for
        Euler–Euler large eddy simulation (EELES) of dispersed bubbly flow.
        Chemical Engineering Science, 63(15), 3923-3931.
    \endverbatim

    The default model coefficients are:
    \verbatim
        NicenoKEqnCoeffs
        {
            Ck              0.094;
            Ce              1.048;
            alphaInversion  0.3;
            Cp              Ck;
            Cmub            0.6;
        }
    \endverbatim

SourceFiles
    NicenoKEqn.C

\*---------------------------------------------------------------------------*/

#ifndef NicenoKEqn_H
#define NicenoKEqn_H

#include "kEqn.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

/*---------------------------------------------------------------------------*\
                         Class NicenoKEqn Declaration
\*---------------------------------------------------------------------------*/

template<class BasicMomentumTransportModel>
class NicenoKEqn
:
    public kEqn<BasicMomentumTransportModel>
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
            dimensionedScalar Cmub_;


    // Protected Member Functions

        virtual void correctNut();
        tmp<volScalarField> bubbleG() const;
        tmp<volScalarField> phaseTransferCoeff() const;
        virtual tmp<fvScalarMatrix> kSource() const;


public:

    typedef typename BasicMomentumTransportModel::alphaField alphaField;
    typedef typename BasicMomentumTransportModel::rhoField rhoField;
    typedef typename BasicMomentumTransportModel::transportModel transportModel;


    //- Runtime type information
    TypeName("NicenoKEqn");


    // Constructors

        //- Construct from components
        NicenoKEqn
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
        NicenoKEqn(const NicenoKEqn&) = delete;


    //- Destructor
    virtual ~NicenoKEqn()
    {}


    // Member Functions

        //- Read model coefficients if they have changed
        virtual bool read();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const NicenoKEqn&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "NicenoKEqn.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
