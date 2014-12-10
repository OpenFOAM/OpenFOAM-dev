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
    Foam::incompressible::LESModels::homogeneousDynSmagorinsky

Group
    grpIcoLESTurbulence

Description
    The Isochoric homogeneous dynamic Smagorinsky Model for
    incompressible flows.

    Algebraic eddy viscosity SGS model founded on the assumption that
    local equilibrium prevails.
    Thus,
    \verbatim
        B = 2/3*k*I - 2*nuSgs*dev(D)
        Beff = 2/3*k*I - 2*nuEff*dev(D)

    where

        k = cI*delta^2*||D||^2
        nuSgs = cD*delta^2*||D||
        nuEff = nuSgs + nu

    In the dynamic version of the choric  Smagorinsky model
    the coefficients cI and cD are calculated during the simulation,

        cI=<K*m>/<m*m>

    and

        cD=1/2*<L.M>/<M.M>,

    where

        K = 0.5*(F(U.U) - F(U).F(U))
        m = delta^2*(4*||F(D)||^2 - F(||D||^2))
        L = dev(F(U*U) - F(U)*F(U))
        M = delta^2*(F(||D||*dev(D)) - 4*||F(D)||*F(dev(D)))

    The averaging <...> is over the whole domain, i.e. homogeneous turbulence
    is assumed
    \endverbatim

SourceFiles
    homogeneousDynSmagorinsky.C

\*---------------------------------------------------------------------------*/

#ifndef homogeneousDynSmagorinsky_H
#define homogeneousDynSmagorinsky_H

#include "Smagorinsky.H"
#include "LESfilter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace LESModels
{

/*---------------------------------------------------------------------------*\
                           Class homogeneousDynSmagorinsky Declaration
\*---------------------------------------------------------------------------*/

class homogeneousDynSmagorinsky
:
    public GenEddyVisc
{
    // Private data

        volScalarField k_;

        autoPtr<LESfilter> filterPtr_;
        LESfilter& filter_;


    // Private Member Functions

        //- Update sub-grid scale fields
        void updateSubGridScaleFields(const volSymmTensorField& D);

        //- Calculate coefficients cD, cI from filtering velocity field
        dimensionedScalar cD(const volSymmTensorField& D) const;
        dimensionedScalar cI(const volSymmTensorField& D) const;

        // Disallow default bitwise copy construct and assignment
        homogeneousDynSmagorinsky(const homogeneousDynSmagorinsky&);
        homogeneousDynSmagorinsky& operator=(const homogeneousDynSmagorinsky&);


public:

    //- Runtime type information
    TypeName("homogeneousDynSmagorinsky");

    // Constructors

        //- Construct from components
        homogeneousDynSmagorinsky
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            transportModel& transport,
            const word& turbulenceModelName = turbulenceModel::typeName,
            const word& modelName = typeName
        );


    //- Destructor
    virtual ~homogeneousDynSmagorinsky()
    {}


    // Member Functions

        //- Return SGS kinetic energy
        virtual tmp<volScalarField> k() const
        {
            return k_;
        }

        //- Correct Eddy-Viscosity and related properties
        virtual void correct(const tmp<volTensorField>& gradU);

        //- Read LESProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
