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
    Foam::incompressible::LESModels::LRRDiffStress

Group
    grpIcoLESTurbulence

Description
    Differential SGS Stress Equation Model for incompressible flows.

    The DSEM uses a model version of the full balance equation for the SGS
    stress tensor to simulate the behaviour of B, hence,
    \verbatim
        d/dt(B) + div(U*B) - div(nuSgs*grad(B))
        =
        P - c1*e/k*B - 0.667*(1 - c1)*e*I - c2*(P - 0.333*trP*I)
    where
        k = 0.5*trB,
        epsilon = ce*k^3/2/delta
        epsilon/k = ce*k^1/2/delta
        P = -(B'L + L'B)
        nuEff = ck*sqrt(k)*delta + nu
    \endverbatim

    This version from Launder, Rece & Rodi 1975

SourceFiles
    LRRDiffStress.C

\*---------------------------------------------------------------------------*/

#ifndef LRRDiffStress_H
#define LRRDiffStress_H

#include "GenSGSStress.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace LESModels
{

/*---------------------------------------------------------------------------*\
                           Class LRRDiffStress Declaration
\*---------------------------------------------------------------------------*/

class LRRDiffStress
:
    public GenSGSStress
{
    // Private data

        dimensionedScalar ck_;
        dimensionedScalar c1_;
        dimensionedScalar c2_;


    // Private Member Functions

        //- Update sub-grid scale fields
        void updateSubGridScaleFields(const volScalarField& K);

        // Disallow default bitwise copy construct and assignment
        LRRDiffStress(const LRRDiffStress&);
        LRRDiffStress& operator=(const LRRDiffStress&);


public:

    //- Runtime type information
    TypeName("LRRDiffStress");

    // Constructors

        //- Construct from components
        LRRDiffStress
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            transportModel& transport,
            const word& turbulenceModelName = turbulenceModel::typeName,
            const word& modelName = typeName
        );


    //- Destructor
    virtual ~LRRDiffStress()
    {}


    // Member Functions

        //- Return the effective diffusivity for B
        tmp<volScalarField> DBEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField("DBEff", nuSgs_ + nu())
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
