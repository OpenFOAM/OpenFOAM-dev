/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2025 OpenFOAM Foundation
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

#include "generatedCellZone.H"
#include "polyMesh.H"
#include "containsPoints.H"

Foam::labelUList Foam::generatedCellZone::identityMap(const label len) const
{
    // Static identity list, resized as required
    static labelList map;

    if (len > map.size())
    {
        map.resize(len);

        forAll(map, i)
        {
            map[i] = i;
        }
    }

    return SubList<label>(map, len);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::generatedCellZone::generatedCellZone(const polyMesh& mesh)
:
    mesh_(mesh),
    all_(true)
{}


Foam::generatedCellZone::generatedCellZone
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    mesh_(mesh),
    all_(true)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::generatedCellZone::~generatedCellZone()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::generatedCellZone::movePoints()
{
    if (!all())
    {
        cellZone_.movePoints();
    }
}


void Foam::generatedCellZone::topoChange(const polyTopoChangeMap& map)
{
    if (!all())
    {
        cellZone_.topoChange(map);
    }
}


void Foam::generatedCellZone::mapMesh(const polyMeshMap& map)
{
    if (!all())
    {
        cellZone_.mapMesh(map);
    }
}


void Foam::generatedCellZone::distribute(const polyDistributionMap& map)
{
    if (!all())
    {
        cellZone_.distribute(map);
    }
}


bool Foam::generatedCellZone::read(const dictionary& dict)
{
    if (dict.found("cellZone"))
    {
        if (!dict.isDict("cellZone"))
        {
            if (dict.lookup<word>("cellZone") == "all")
            {
                all_ = true;
            }
        }
        else
        {
            all_ = false;

            cellZone_.read
            (
                "cellZone",
                zoneGenerator::zoneTypes::cell,
                mesh_,
                dict
            );
        }
    }
    else if (dict.found("points"))
    {
        // For backward compatibility

        IOWarningInFunction(dict)
            << "points is deprecated, use cellZone instead." << nl
            << "    For backward compatibility the points entry "
               "is automatically converted into a cellZone generated "
               "using the containsPoints zoneGenerator."
            << endl;

        all_ = false;

        cellZone_.set
        (
            autoPtr<zoneGenerator>
            (
                new zoneGenerators::containsPoints("points", mesh_, dict)
            )
        );
    }
    else if (dict.found("select"))
    {
        const word selection(dict.lookup("select"));

        IOWarningInFunction(dict)
            << "select " << selection
            << " is deprecated, use cellZone instead."
            << endl;

        if (selection == "all")
        {
            all_ = true;
        }
    }
    else
    {
        FatalIOErrorInFunction(dict)
            << "cellZone not specified"
            << exit(FatalIOError);
    }

    return true;
}


// ************************************************************************* //
