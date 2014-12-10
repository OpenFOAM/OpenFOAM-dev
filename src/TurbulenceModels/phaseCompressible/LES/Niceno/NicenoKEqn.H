/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2014 OpenFOAM Foundation
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

Group
    grpLESTurbulence

Description
    One-equation SGS model for the continuous phase in a two-phase system
    including bubble-generated turbulence.

    Reference:
    \verbatim
        "One-equation sub-grid scale (SGS) modelling for Euler-Euler
        large eddy simulation (EELES) of dispersed bubbly flow"
        B. Niceno,
        M.T. Dhotre,
        N.G. Dee
        Chemical Engineering Science 63 (2008) pp. 3923-3931.
    \endverbatim

    The default model coefficients correspond to the following:
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

template<class BasicTurbulenceModel>
class NicenoKEqn
:
    public kEqn<BasicTurbulenceModel>
{
    // Private data

        mutable const PhaseCompressibleTurbulenceModel
        <
            typename BasicTurbulenceModel::transportModel
        > *gasTurbulencePtr_;


    // Private Member Functions

        //- Return the turbulence model for the gas phase
        const PhaseCompressibleTurbulenceModel
        <
            typename BasicTurbulenceModel::transportModel
        >&
        gasTurbulence() const;

        // Disallow default bitwise copy construct and assignment
        NicenoKEqn(const NicenoKEqn&);
        NicenoKEqn& operator=(const NicenoKEqn&);


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

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;


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
            const word& propertiesName = turbulenceModel::propertiesName,
            const word& type = typeName
        );


    //- Destructor
    virtual ~NicenoKEqn()
    {}


    // Member Functions

        //- Read model coefficients if they have changed
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "NicenoKEqn.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
