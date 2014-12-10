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
    Foam::incompressible::LESModels::spectEddyVisc

Group
    grpIcoLESTurbulence

Description
    The Isochoric spectral Eddy Viscosity Model for incompressible flows.

    Algebraic eddy viscosity SGS model founded on the assumption that
    local equilibrium prevail and that viscous effects may be of importance.

    Thus,
    \verbatim
        B = 2/3*k*I - 2*nuSgs*dev(D)
        Beff = 2/3*k*I - 2*nuEff*dev(D)

    where

        k = cK1*delta^(2/3)*eps^(2/3)*exp(-cK2*delta^(4/3)*nu*eps^(-1/3)) -
            ck3*sqrt(eps*nu)*erfc(ck4*delta^(-2/3)*sqrt(nu)*eps^(-1/6)))

        nuEff = nu/(1 - exp(-cB*pow((nu/nuEff), 1.0/3.0)*pow(Re, -2.0/3.0)))
        nuSgs = nuEff - nu

        Re = delta^2*mag(D)/nu
    \endverbatim

SourceFiles
    spectEddyVisc.C

\*---------------------------------------------------------------------------*/

#ifndef spectEddyVisc_H
#define spectEddyVisc_H

#include "GenEddyVisc.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace LESModels
{

/*---------------------------------------------------------------------------*\
                           Class spectEddyVisc Declaration
\*---------------------------------------------------------------------------*/

class spectEddyVisc
:
    public GenEddyVisc
{
    // Private data

        dimensionedScalar cB_;
        dimensionedScalar cK1_;
        dimensionedScalar cK2_;
        dimensionedScalar cK3_;
        dimensionedScalar cK4_;


    // Private Member Functions

        //- Update sub-grid scale fields
        void updateSubGridScaleFields(const volTensorField& gradU);

        // Disallow default bitwise copy construct and assignment
        spectEddyVisc(const spectEddyVisc&);
        spectEddyVisc& operator=(const spectEddyVisc&);


public:

    //- Runtime type information
    TypeName("spectEddyVisc");

    // Constructors

        // Construct from components
        spectEddyVisc
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            transportModel& transport,
            const word& turbulenceModelName = turbulenceModel::typeName,
            const word& modelName = typeName
        );


    //- Destructor
    virtual ~spectEddyVisc()
    {}


    // Member Functions

        //- Return SGS kinetic energy
        virtual tmp<volScalarField> k() const;

        //- Correct Eddy-Viscosity and related properties
        virtual void correct(const tmp<volTensorField>&);

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
