/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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

#include "LagrangianInjection.H"
#include "LagrangianMeshLocation.H"
#include "timeIOdictionary.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(LagrangianInjection, 0);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::LagrangianInjection::checkLocation
(
    const LagrangianMesh::location location,
    const point& p
) const
{
    switch (location)
    {
        case LagrangianMesh::location::inCell:
            break;

        case LagrangianMesh::location::onBoundary:
            WarningInFunction
                << "Did not accurately locate the injection point "
                << p << " within the mesh" << endl;
            break;

        case LagrangianMesh::location::outsideMesh:
            FatalErrorInFunction
                << "Did not find the injection point "
                << p << " within the mesh"
                << exit(FatalError);
            break;
    }
}


void Foam::LagrangianInjection::checkLocation
(
    const List<LagrangianMesh::location>& locations,
    const List<point>& points
) const
{
    if (!locations.size()) return;

    LagrangianMesh::location l = locations[0];
    point p = points[0];

    for (label i = 1; i < locations.size(); ++ i)
    {
        if (l != LagrangianMeshLocation::furthestOp()(l, locations[i]))
        {
            l = locations[i];
            p = points[i];
        }
    }

    checkLocation(l, p);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::LagrangianInjection::LagrangianInjection
(
    const word& name,
    const LagrangianMesh& mesh
)
:
    LagrangianModel(name, mesh)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::LagrangianInjection::~LagrangianInjection()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::LagrangianInjection::addsSupToField(const word& fieldName) const
{
    return false;
}


// ************************************************************************* //
