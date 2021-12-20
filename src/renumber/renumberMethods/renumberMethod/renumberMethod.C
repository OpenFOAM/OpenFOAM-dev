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

InClass
    renumberMethod

\*---------------------------------------------------------------------------*/

#include "renumberMethod.H"
#include "decompositionMethod.H"
#include "dlLibraryTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(renumberMethod, 0);
    defineRunTimeSelectionTable(renumberMethod, dictionary);
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::autoPtr<Foam::renumberMethod> Foam::renumberMethod::New
(
    const dictionary& renumberDict
)
{
    const word methodType(renumberDict.lookup("method"));

    Info<< "Selecting renumberMethod " << methodType << endl;

    libs.open
    (
        renumberDict,
        "libs",
        dictionaryConstructorTablePtr_
    );

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(methodType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown renumberMethod "
            << methodType << nl << nl
            << "Valid renumberMethods are : " << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<renumberMethod>(cstrIter()(renumberDict));
}


Foam::labelList Foam::renumberMethod::renumber
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

    // Renumber based on agglomerated points
    return renumber(cellCells(), points);
}


Foam::labelList Foam::renumberMethod::renumber
(
    const labelList& cellCells,
    const labelList& offsets,
    const pointField& cc
) const
{
    NotImplemented;
    return labelList();
}


Foam::labelList Foam::renumberMethod::renumber
(
    const polyMesh& mesh,
    const labelList& fineToCoarse,
    const pointField& coarsePoints
) const
{
    CompactListList<label> coarseCellCells;
    decompositionMethod::calcCellCells
    (
        mesh,
        fineToCoarse,
        coarsePoints.size(),
        false,                      // local only
        coarseCellCells
    );

    // Renumber based on agglomerated points
    labelList coarseDistribution
    (
        renumber
        (
            coarseCellCells(),
            coarsePoints
        )
    );

    // Rework back into renumbering for original mesh_
    labelList fineDistribution(fineToCoarse.size());

    forAll(fineDistribution, i)
    {
        fineDistribution[i] = coarseDistribution[fineToCoarse[i]];
    }

    return fineDistribution;
}


// ************************************************************************* //
