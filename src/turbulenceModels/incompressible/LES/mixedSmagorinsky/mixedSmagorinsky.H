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
    Foam::incompressible::LESModels::mixedSmagorinsky

Group
    grpIcoLESTurbulence

Description
    The mixed Isochoric Smagorinsky Model for incompressible flows.

    The mixed model is a linear combination of an eddy viscosity model
    (Smagorinsky) with a scale similarity model.  Hence
    \verbatim
        B = (L + C) + R = (F(v*v) - F(v)*F(v)) + R
    \endverbatim

    The algebraic eddy viscosity SGS model is founded on the assumption
    that local equilibrium prevails, hence
    \verbatim
        R = 2/3*k*I - 2*nuEff*dev(D)
    where
        k = cI*delta^2*||D||^2
        nuEff = ck*sqrt(k)*delta + nu
    \endverbatim

    The Leonard and cross contributions are incorporated
    by adding,
    \verbatim
         + div(((filter(U*U) - filter(U)*filter(U)) -
           0.333*I*tr(filter(U*U) - filter(U)*filter(U))))
         + div((filter(U*epsilon) - filter(U)*filter(epsilon)))
    \endverbatim
    to the rhs. of the equations.

SourceFiles
    mixedSmagorinsky.C

\*---------------------------------------------------------------------------*/

#ifndef mixedSmagorinsky_H
#define mixedSmagorinsky_H

#include "scaleSimilarity.H"
#include "Smagorinsky.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace LESModels
{

/*---------------------------------------------------------------------------*\
                           Class mixedSmagorinsky Declaration
\*---------------------------------------------------------------------------*/

class mixedSmagorinsky
:
    public scaleSimilarity,
    public Smagorinsky
{
    // Private Member Functions

        // Disallow default bitwise copy construct and assignment
        mixedSmagorinsky(const mixedSmagorinsky&);
        mixedSmagorinsky& operator=(const mixedSmagorinsky&);

public:

    //- Runtime type information
    TypeName("mixedSmagorinsky");


    // Constructors

        //- Construct from components
        mixedSmagorinsky
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            transportModel& transport,
            const word& turbulenceModelName = turbulenceModel::typeName,
            const word& modelName = typeName
        );


    //- Destructor
    virtual ~mixedSmagorinsky()
    {}


    // Member Functions

        //- Return the SGS turbulent kinetic energy.
        virtual tmp<volScalarField> k() const;

        //- Return the SGS turbulent disipation rate.
        virtual tmp<volScalarField> epsilon() const;

        //- Return the SGS viscosity.
        virtual tmp<volScalarField> nuSgs() const
        {
            return nuSgs_;
        }

        //- Return the sub-grid stress tensor.
        virtual tmp<volSymmTensorField> B() const;

        //- Return the effective sub-grid turbulence stress tensor
        //  including the laminar stress
        virtual tmp<volSymmTensorField> devReff() const;

        //- Return the deviatoric part of the effective sub-grid
        //  turbulence stress tensor including the laminar stress
        virtual tmp<fvVectorMatrix> divDevReff(volVectorField& U) const;

        //- Return the deviatoric part of the effective sub-grid
        //  turbulence stress tensor including the laminar stress
        virtual tmp<fvVectorMatrix> divDevRhoReff
        (
            const volScalarField& rho,
            volVectorField& U
        ) const;

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
