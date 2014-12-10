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
    Foam::incompressible::LESModels::homogeneousDynOneEqEddy

Group
    grpIcoLESTurbulence

Description
    One Equation Eddy Viscosity Model for incompressible flows.

    Eddy viscosity SGS model using a modeled balance equation to simulate
    the behaviour of k.

    Thus
    \verbatim
        d/dt(k) + div(U*k) - div(nuSgs*grad(k))
        =
        -B*L - ce*k^3/2/delta

    and

        B = 2/3*k*I - 2*nuSgs*dev(D)
        Beff = 2/3*k*I - 2*nuEff*dev(D)

    where

        D = symm(grad(U));
        nuSgs = ck*sqrt(k)*delta
        nuEff = nuSgs + nu
    \endverbatim

SourceFiles
    homogeneousDynOneEqEddy.C

\*---------------------------------------------------------------------------*/

#ifndef homogeneousDynOneEqEddy_H
#define homogeneousDynOneEqEddy_H

#include "GenEddyVisc.H"
#include "LESfilter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace LESModels
{

/*---------------------------------------------------------------------------*\
                           Class homogeneousDynOneEqEddy Declaration
\*---------------------------------------------------------------------------*/

class homogeneousDynOneEqEddy
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

        //- Calculate ck, ce by filtering the velocity field U.
        dimensionedScalar ck(const volSymmTensorField& D) const;
        dimensionedScalar ce(const volSymmTensorField& D) const;

        // Disallow default bitwise copy construct and assignment
        homogeneousDynOneEqEddy(const homogeneousDynOneEqEddy&);
        homogeneousDynOneEqEddy& operator=(const homogeneousDynOneEqEddy&);


public:

    //- Runtime type information
    TypeName("homogeneousDynOneEqEddy");

    // Constructors

        //- Construct from components
        homogeneousDynOneEqEddy
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            transportModel& transport,
            const word& turbulenceModelName = turbulenceModel::typeName,
            const word& modelName = typeName
        );


    //- Destructor
    virtual ~homogeneousDynOneEqEddy()
    {}


    // Member Functions

        //- Return SGS kinetic energy
        virtual tmp<volScalarField> k() const
        {
            return k_;
        }

        //- Return the effective diffusivity for k
        tmp<volScalarField> DkEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField("DkEff", nuSgs_ + nu())
            );
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
