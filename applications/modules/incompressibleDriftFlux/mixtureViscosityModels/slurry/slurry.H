/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2024 OpenFOAM Foundation
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
    Foam::mixtureViscosityModels::slurry

Description
    Thomas' viscosity correction for slurry.

    Reference:
    \verbatim
        "Transport characteristics of suspension:
         VIII. A note on the viscosity of Newtonian suspensions
         of uniform spherical particles".
         D.G. Thomas,
         J. Colloid Sci. 20 (3), 1965, p267.
    \endverbatim

SourceFiles
    slurry.C

\*---------------------------------------------------------------------------*/

#ifndef slurry_H
#define slurry_H

#include "mixtureViscosityModel.H"
#include "dimensionedScalar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class incompressibleDriftFluxMixture;

namespace mixtureViscosityModels
{

/*---------------------------------------------------------------------------*\
                           Class slurry Declaration
\*---------------------------------------------------------------------------*/

class slurry
:
    public mixtureViscosityModel
{

public:

    //- Runtime type information
    TypeName("slurry");


    // Constructors

        //- Construct from mixture
        slurry(const incompressibleDriftFluxMixture& mixture);


    //- Destructor
    virtual ~slurry()
    {}


    // Member Functions

        //- Return the mixture viscosity
        //  given the viscosity of the continuous phase
        virtual tmp<volScalarField> mu
        (
            const volScalarField& muc,
            const volVectorField& U
        ) const;

        //- Read phaseProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace mixtureViscosityModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
