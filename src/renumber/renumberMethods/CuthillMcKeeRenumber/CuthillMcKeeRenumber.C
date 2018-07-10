/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "CuthillMcKeeRenumber.H"
#include "addToRunTimeSelectionTable.H"
#include "bandCompression.H"
#include "decompositionMethod.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(CuthillMcKeeRenumber, 0);

    addToRunTimeSelectionTable
    (
        renumberMethod,
        CuthillMcKeeRenumber,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::CuthillMcKeeRenumber::CuthillMcKeeRenumber(const dictionary& renumberDict)
:
    renumberMethod(renumberDict),
    reverse_
    (
        renumberDict.optionalSubDict
        (
            typeName + "Coeffs"
        ).lookupOrDefault<Switch>("reverse", false)
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::labelList Foam::CuthillMcKeeRenumber::renumber
(
    const polyMesh& mesh,
    const pointField& points
) const
{
    CompactListList<label> cellCells;
    decompositionMethod::calcCellCells
    (
        mesh,
        identity(mesh.nCells()),
        mesh.nCells(),
        false,                      // local only
        cellCells
    );

    labelList orderedToOld = bandCompression(cellCells());

    if (reverse_)
    {
        reverse(orderedToOld);
    }

    return orderedToOld;
}


Foam::labelList Foam::CuthillMcKeeRenumber::renumber
(
    const labelList& cellCells,
    const labelList& offsets,
    const pointField& cc
) const
{
    labelList orderedToOld = bandCompression(cellCells, offsets);

    if (reverse_)
    {
        reverse(orderedToOld);
    }

    return orderedToOld;
}


Foam::labelList Foam::CuthillMcKeeRenumber::renumber
(
    const labelListList& cellCells,
    const pointField& points
) const
{
    labelList orderedToOld = bandCompression(cellCells);

    if (reverse_)
    {
        reverse(orderedToOld);
    }

    return orderedToOld;
}


// ************************************************************************* //
