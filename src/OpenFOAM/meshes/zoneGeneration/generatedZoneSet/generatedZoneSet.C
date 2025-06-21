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

#include "generatedZoneSet.H"
#include "lookup.H"
#include "polyMesh.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::generatedZoneSet::generatedZoneSet
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
{
    read(name, mesh, dict);
}


Foam::generatedZoneSet::generatedZoneSet
(
    const word& name,
    const zoneTypes& zoneType,
    const polyMesh& mesh,
    const dictionary& dict
)
{
    read(name, zoneType, mesh, dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::generatedZoneSet::read
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
{
    if (dict.isDict(name))
    {
        zoneGenerator_ = zoneGenerator::New(name, mesh, dict.subDict(name));
    }
    else
    {
        Istream& zoneStream = dict.lookup(name);
        const word zoneName(zoneStream);

        // If an empty keyword is present assume it is a zone name
        // and add a zone lookup
        dictionary zoneDict(name, dict);
        zoneDict.add
        (
            primitiveEntry
            (
                "type",
                zoneGenerators::lookup::typeName,
                zoneStream.lineNumber()
            )
        );

        zoneGenerator_ = new zoneGenerators::lookup(zoneName, mesh, zoneDict);
    }

    zoneSet::operator=(zoneGenerator_->generate());
}


void Foam::generatedZoneSet::read
(
    const word& name,
    const zoneTypes& zoneType,
    const polyMesh& mesh,
    const dictionary& dict
)
{
    if (dict.isDict(name))
    {
        zoneGenerator_ = zoneGenerator::New
        (
            name,
            zoneType,
            mesh,
            dict.subDict(name));
    }
    else
    {
        Istream& zoneStream = dict.lookup(name);
        const word zoneName(zoneStream);

        // If an empty keyword is present assume it is a zone name
        // and add a zone lookup
        dictionary zoneDict(name, dict);
        zoneDict.add
        (
            primitiveEntry
            (
                "type",
                zoneGenerators::lookup::typeName,
                zoneStream.lineNumber()
            )
        );
        zoneDict.add("zoneType", zoneTypesNames[zoneType]);

        zoneGenerator_ = new zoneGenerators::lookup(zoneName, mesh, zoneDict);
    }

    zoneSet::operator=(zoneGenerator_->generate());
}


void Foam::generatedZoneSet::set(const autoPtr<zoneGenerator>& zg)
{
    zoneGenerator_ = zg;
    zoneSet::operator=(zoneGenerator_->generate());
}


bool Foam::generatedZoneSet::movePoints()
{
    if (zoneGenerator_->moveUpdate())
    {
        zoneSet::operator=(zoneGenerator_->generate());
    }

    return true;
}


void Foam::generatedZoneSet::distribute(const polyDistributionMap&)
{
    // Regenerate zones following redistribution
    zoneSet::operator=(zoneGenerator_->generate());
}


void Foam::generatedZoneSet::topoChange(const polyTopoChangeMap&)
{
    // Regenerate zones following topology change
    zoneSet::operator=(zoneGenerator_->generate());
}


void Foam::generatedZoneSet::mapMesh(const polyMeshMap&)
{
    // Zones are automatically mapped but should be regenerated
    zoneSet::operator=(zoneGenerator_->generate());
}


// ************************************************************************* //
