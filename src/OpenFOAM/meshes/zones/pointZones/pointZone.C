/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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

#include "pointZone.H"
#include "pointZones.H"
#include "polyMesh.H"
#include "polyTopoChangeMap.H"
#include "syncTools.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    typedef Zone<pointZone, pointZones> pointZoneType;
    defineTemplateRunTimeSelectionTable(pointZoneType, dictionary);

    defineTypeNameAndDebug(pointZone, 0);
    addToRunTimeSelectionTable(pointZone, pointZone, dictionary);
}

const char* const Foam::pointZone::labelsName = "pointLabels";


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::pointZone::checkDefinition(const bool report) const
{
    return Zone::checkDefinition
    (
        zones_.mesh().points().size(),
        report
    );
}


bool Foam::pointZone::checkParallelSync(const bool report) const
{
    const polyMesh& mesh = zones().mesh();

    const label index = zones_.findIndex(name());

    labelList maxZone(mesh.nPoints(), -1);
    labelList minZone(mesh.nPoints(), labelMax);
    forAll(*this, i)
    {
        const label pointi = operator[](i);
        maxZone[pointi] = index;
        minZone[pointi] = index;
    }
    syncTools::syncPointList(mesh, maxZone, maxEqOp<label>(), label(-1));
    syncTools::syncPointList(mesh, minZone, minEqOp<label>(), labelMax);

    bool error = false;

    forAll(maxZone, pointi)
    {
        // Check point in same (or no) zone on all processors
        if
        (
            (
                maxZone[pointi] != -1
             || minZone[pointi] != labelMax
            )
         && (maxZone[pointi] != minZone[pointi])
        )
        {
            if (report && !error)
            {
                Info<< " ***Problem with pointZone " << name()
                    << ". Point " << pointi
                    << " at " << mesh.points()[pointi]
                    << " is in zone "
                    << (minZone[pointi] == labelMax ? -1 : minZone[pointi])
                    << " on some processors and in zone "
                    << maxZone[pointi]
                    << " on some other processors." << nl
                    << "(suppressing further warnings)"
                    << endl;
            }
            error = true;
        }
    }

    return error;
}


void Foam::pointZone::topoChange(const polyTopoChangeMap& map)
{
    Zone::topoChange(map.pointMap(), map.reversePointMap());
}


void Foam::pointZone::writeDict(Ostream& os) const
{
    os  << nl << name_ << nl << token::BEGIN_BLOCK << nl
        << "    type " << type() << token::END_STATEMENT << nl;

    writeEntry(os, this->labelsName, *this);

    os  << token::END_BLOCK << endl;
}


// ************************************************************************* //
