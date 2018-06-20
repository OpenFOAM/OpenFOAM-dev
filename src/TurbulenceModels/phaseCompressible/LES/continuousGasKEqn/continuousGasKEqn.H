/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2018 OpenFOAM Foundation
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
    Foam::LESModels::continuousGasKEqn

Description
    One-equation SGS model for the gas-phase in a two-phase system
    supporting phase-inversion.

    In the limit that the gas-phase fraction approaches zero a contribution from
    the other phase is blended into the k-equation up to the phase-fraction of
    alphaInversion at which point phase-inversion is considered to have occurred
    and the model reverts to the pure single-phase form.

    This model is unpublished and is provided as a stable numerical framework
    on which a more physical model may be built.

    The default model coefficients are
    \verbatim
        continuousKEqnCoeffs
        {
            Ck              0.094;
            Ce              1.048;
            alphaInversion  0.7;
        }
    \endverbatim

SourceFiles
    continuousGasKEqn.C

\*---------------------------------------------------------------------------*/

#ifndef continuousGasKEqn_H
#define continuousGasKEqn_H

#include "kEqn.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

/*---------------------------------------------------------------------------*\
                           Class continuousGasKEqn Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class continuousGasKEqn
:
    public kEqn<BasicTurbulenceModel>
{
    // Private data

        mutable const turbulenceModel *liquidTurbulencePtr_;


    // Private Member Functions

        // Disallow default bitwise copy construct and assignment
        continuousGasKEqn(const continuousGasKEqn&);
        void operator=(const continuousGasKEqn&);


protected:

    // Protected data

        // Model coefficients

            dimensionedScalar alphaInversion_;


    // Protected Member Functions

        //- Return the turbulence model for the liquid phase
        const turbulenceModel& liquidTurbulence() const;

        tmp<volScalarField> phaseTransferCoeff() const;
        virtual tmp<fvScalarMatrix> kSource() const;


public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;


    //- Runtime type information
    TypeName("continuousGasKEqn");


    // Constructors

        //- Construct from components
        continuousGasKEqn
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
    virtual ~continuousGasKEqn()
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
    #include "continuousGasKEqn.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
