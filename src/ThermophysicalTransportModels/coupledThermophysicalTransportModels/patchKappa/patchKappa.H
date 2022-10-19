/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022 OpenFOAM Foundation
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
    Foam::patchKappa

Description
    Provides the patch thermal conductivity \c kappa for the fluid or solid
    region an appropriate.

    For fluid regions \c kappaEff is obtained from the
    Foam::thermophysicalTransportModel and for solid regions
    \c kappa is obtained from the Foam::solidThermophysicalTransportModel.

SourceFiles
    patchKappa.C

\*---------------------------------------------------------------------------*/

#ifndef patchKappa_H
#define patchKappa_H

#include "fvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class patchKappa Declaration
\*---------------------------------------------------------------------------*/

class patchKappa
{
    // Private data

        //- Reference to the region patch
        const fvPatch& patch_;


public:

    // Constructors

        //- Construct from patch
        patchKappa
        (
            const fvPatch& patch
        );


    // Member Functions

        //- Thermal conductivity for patch [W/m/K]
        tmp<scalarField> kappa() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
