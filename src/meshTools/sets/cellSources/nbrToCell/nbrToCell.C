/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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

#include "nbrToCell.H"
#include "polyMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(nbrToCell, 0);
    addToRunTimeSelectionTable(topoSetSource, nbrToCell, word);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::nbrToCell::combine(topoSet& set, const bool add) const
{
    const cellList& cells = mesh().cells();
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    boolList isCoupled(mesh_.nFaces()-mesh_.nInternalFaces(), false);

    forAll(patches, patchi)
    {
        const polyPatch& pp = patches[patchi];

        if (pp.coupled())
        {
            label facei = pp.start();
            forAll(pp, i)
            {
                isCoupled[facei-mesh_.nInternalFaces()] = true;
                facei++;
            }
        }
    }

    forAll(cells, celli)
    {
        const cell& cFaces = cells[celli];

        label nNbrCells = 0;

        forAll(cFaces, i)
        {
            label facei = cFaces[i];

            if (mesh_.isInternalFace(facei))
            {
                nNbrCells++;
            }
            else if (isCoupled[facei-mesh_.nInternalFaces()])
            {
                nNbrCells++;
            }
        }

        if (nNbrCells <= minNbrs_)
        {
            addOrDelete(set, celli, add);
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::nbrToCell::nbrToCell
(
    const polyMesh& mesh,
    const label minNbrs
)
:
    topoSetSource(mesh),
    minNbrs_(minNbrs)
{}


Foam::nbrToCell::nbrToCell
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    topoSetSource(mesh),
    minNbrs_(dict.lookup<label>("neighbours"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::nbrToCell::~nbrToCell()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::nbrToCell::applyToSet
(
    const topoSetSource::setAction action,
    topoSet& set
) const
{
    if ((action == topoSetSource::NEW) || (action == topoSetSource::ADD))
    {
        Info<< "    Adding cells with only " << minNbrs_ << " or less"
                " neighbouring cells" << " ..." << endl;

        combine(set, true);
    }
    else if (action == topoSetSource::DELETE)
    {
        Info<< "    Removing cells with only " << minNbrs_ << " or less"
                " neighbouring cells" << " ..." << endl;

        combine(set, false);
    }
}


// ************************************************************************* //
