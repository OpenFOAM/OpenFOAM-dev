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

#include "point_zoneGenerator.H"
#include "polyMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace zoneGenerators
    {
        defineTypeNameAndDebug(point, 0);
        addToRunTimeSelectionTable
        (
            zoneGenerator,
            point,
            dictionary
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::zoneGenerators::point::point
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    zoneGenerator(name, mesh, dict),
    zoneGenerators_(mesh, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::zoneGenerators::point::~point()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::zoneSet Foam::zoneGenerators::point::generate() const
{
    boolList selectedPoints(mesh_.nPoints(), false);

    forAll(zoneGenerators_, i)
    {
        zoneSet zs(zoneGenerators_[i].generate());

        if (zs.pValid() && zs.pZone().name() != zoneName_)
        {
            const labelList& zonePoints = zs.pZone();

            forAll(zonePoints, zpi)
            {
                selectedPoints[zonePoints[zpi]] = true;
            }
        }

        if (zs.cValid())
        {
            const labelList& zoneCells = zs.cZone();

            forAll(zoneCells, zci)
            {
                const labelList& cellFaces = mesh_.cells()[zoneCells[zci]];

                forAll(cellFaces, cFacei)
                {
                    const face& f = mesh_.faces()[cellFaces[cFacei]];

                    forAll(f, fp)
                    {
                        selectedPoints[f[fp]] = true;
                    }
                }
            }
        }

        if (zs.fValid())
        {
            const labelList& zoneFaces = zs.fZone();

            forAll(zoneFaces, zfi)
            {
                const face& f = mesh_.faces()[zoneFaces[zfi]];

                forAll(f, fp)
                {
                    selectedPoints[f[fp]] = true;
                }
            }
        }
    }

    moveUpdate_ = zoneGenerators_.moveUpdate();

    return zoneSet
    (
        new pointZone
        (
            zoneName_,
            indices(selectedPoints),
            mesh_.pointZones(),
            moveUpdate_,
            true
        )
    );
}


// ************************************************************************* //
