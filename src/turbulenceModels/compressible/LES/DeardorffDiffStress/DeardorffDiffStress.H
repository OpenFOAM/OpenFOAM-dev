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
    Foam::compressible::LESModels::DeardorffDiffStress

Group
    grpCmpLESTurbulence

Description
    Differential SGS Stress Equation Model for compressible flows

    The DSEM uses a model version of the full balance equation for the SGS
    stress tensor to simulate the behaviour of B.
    Thus,
    \verbatim
        d/dt(rho*B) + div(rho*U*B) - div(muSgs*grad(B))
          =
            P
          - c1*rho*epsilon/k*B
          - 0.667*(1 - c1)*rho*epsilon*I
          - c2*(P - 0.333*trP*I)

    where

        k = 0.5*trB,
        epsilon = ce*k^3/2/delta,
        epsilon/k = ce*k^1/2/delta
        P = -rho*(B'L + L'B)
        muSgs = ck*rho*sqrt(k)*delta
        muEff = muSgs + mu
    \endverbatim

SourceFiles
    DeardorffDiffStress.C

\*---------------------------------------------------------------------------*/

#ifndef compressibleDeardorffDiffStress_H
#define compressibleDeardorffDiffStress_H

#include "GenSGSStress.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{
namespace LESModels
{

/*---------------------------------------------------------------------------*\
                           Class DeardorffDiffStress Declaration
\*---------------------------------------------------------------------------*/

class DeardorffDiffStress
:
    public GenSGSStress
{
    // Private data

        dimensionedScalar ck_;
        dimensionedScalar cm_;


    // Private Member Functions

        //- Update sub-grid scale fields
        void updateSubGridScaleFields(const volScalarField& K);

        // Disallow default bitwise copy construct and assignment
        DeardorffDiffStress(const DeardorffDiffStress&);
        DeardorffDiffStress& operator=(const DeardorffDiffStress&);


public:

    //- Runtime type information
    TypeName("DeardorffDiffStress");

    // Constructors

        //- Constructor from components
        DeardorffDiffStress
        (
            const volScalarField& rho,
            const volVectorField& U,
            const surfaceScalarField& phi,
            const fluidThermo& thermoPhysicalModel,
            const word& turbulenceModelName = turbulenceModel::typeName,
            const word& modelName = typeName
        );


    //- Destructor
    virtual ~DeardorffDiffStress()
    {}


    // Member Functions

        //- Return the effective diffusivity for B
        tmp<volScalarField> DBEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField("DBEff", muSgs_ + mu())
            );
        }

        //- Correct Eddy-Viscosity and related properties
        virtual void correct(const tmp<volTensorField>& gradU);

        //- Read LESProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace compressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
