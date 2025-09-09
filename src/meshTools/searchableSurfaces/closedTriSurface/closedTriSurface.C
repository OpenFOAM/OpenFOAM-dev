/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2025 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "closedTriSurface.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace searchableSurfaces
    {
        defineTypeNameAndDebug(closedTriSurface, 0);

        addToRunTimeSelectionTable
        (
            searchableSurface,
            closedTriSurface,
            dictionary
        );

        addBackwardCompatibleToRunTimeSelectionTable
        (
            searchableSurface,
            closedTriSurface,
            dictionary,
            closedTriSurfaceMesh,
            "closedTriSurfaceMesh"
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::searchableSurfaces::closedTriSurface::closedTriSurface
(
    const IOobject& io,
    const triSurface& s
)
:
    triSurface(io, s)
{}


Foam::searchableSurfaces::closedTriSurface::closedTriSurface
(
    const IOobject& io
)
:
    triSurface(io)
{}


Foam::searchableSurfaces::closedTriSurface::closedTriSurface
(
    const IOobject& io,
    const dictionary& dict
)
:
    triSurface(io, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::searchableSurfaces::closedTriSurface::~closedTriSurface()
{}


// ************************************************************************* //
