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

#include "patchCells_zoneGenerator.H"
#include "polyMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace zoneGenerators
    {
        defineTypeNameAndDebug(patchCells, 0);
        addToRunTimeSelectionTable
        (
            zoneGenerator,
            patchCells,
            dictionary
        );
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::zoneGenerators::patchCells::patchCells
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    zoneGenerator(name, mesh, dict),
    patchSet_(mesh.boundaryMesh().patchSet(dict))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::zoneGenerators::patchCells::~patchCells()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::zoneSet Foam::zoneGenerators::patchCells::generate() const
{
    boolList selectedCells(mesh_.nCells(), false);

    forAllConstIter(labelHashSet, patchSet_, iter)
    {
        const label patchi = iter.key();
        const polyPatch& pp = mesh_.boundaryMesh()[patchi];
        const labelList& patchCells(pp.faceCells());

        forAll(patchCells, facei)
        {
            selectedCells[patchCells[facei]] = true;
        }
    }

    return zoneSet
    (
        new cellZone
        (
            zoneName_,
            indices(selectedCells),
            mesh_.cellZones(),
            moveUpdate_,
            true
        )
    );
}


// ************************************************************************* //
