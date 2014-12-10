/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

#include "ptscotchDecomp.H"
#include "addToRunTimeSelectionTable.H"
#include "Time.H"

static const char* notImplementedMessage =
"You are trying to use ptscotch but do not have the "
"ptscotchDecomp library loaded."
"\nThis message is from the dummy ptscotchDecomp stub library instead.\n"
"\n"
"Please install ptscotch and make sure that libptscotch.so is in your "
"LD_LIBRARY_PATH.\n"
"The ptscotchDecomp library can then be built in "
"$FOAM_SRC/parallel/decompose/ptscotchDecomp\n";


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(ptscotchDecomp, 0);

    addToRunTimeSelectionTable
    (
        decompositionMethod,
        ptscotchDecomp,
        dictionary
    );
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::ptscotchDecomp::check(const int retVal, const char* str)
{}


Foam::label Foam::ptscotchDecomp::decompose
(
    const fileName& meshPath,
    const List<int>& initxadj,
    const List<int>& initadjncy,
    const scalarField& initcWeights,

    List<int>& finalDecomp
) const
{
    FatalErrorIn
    (
        "label ptscotchDecomp::decompose"
        "("
            "onst fileName&,"
            "const List<int>&, "
            "const List<int>&, "
            "const scalarField&, "
            "List<int>&"
        ")"
    )   << notImplementedMessage << exit(FatalError);

    return -1;
}
Foam::label Foam::ptscotchDecomp::decompose
(
    const fileName& meshPath,
    const int adjncySize,
    const int adjncy[],
    const int xadjSize,
    const int xadj[],
    const scalarField& cWeights,
    List<int>& finalDecomp
) const
{
    FatalErrorIn
    (
        "label ptscotchDecomp::decompose"
        "("
            "const fileName&,"
            "const int,"
            "const int,"
            "const int,"
            "const int,"
            "const scalarField&,"
            "List<int>&"
        ")"
    )   << notImplementedMessage << exit(FatalError);

    return -1;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ptscotchDecomp::ptscotchDecomp
(
    const dictionary& decompositionDict
)
:
    decompositionMethod(decompositionDict)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::labelList Foam::ptscotchDecomp::decompose
(
    const polyMesh& mesh,
    const pointField& points,
    const scalarField& pointWeights
)
{
    FatalErrorIn
    (
        "labelList ptscotchDecomp::decompose"
        "("
            "const pointField&, "
            "const scalarField&"
        ")"
    )   << notImplementedMessage << exit(FatalError);

    return labelList::null();
}


Foam::labelList Foam::ptscotchDecomp::decompose
(
    const polyMesh& mesh,
    const labelList& agglom,
    const pointField& agglomPoints,
    const scalarField& pointWeights
)
{
    FatalErrorIn
    (
        "labelList ptscotchDecomp::decompose"
        "("
            "const labelList&, "
            "const pointField&, "
            "const scalarField&"
        ")"
    )   << notImplementedMessage << exit(FatalError);

    return labelList::null();
}


Foam::labelList Foam::ptscotchDecomp::decompose
(
    const labelListList& globalCellCells,
    const pointField& cellCentres,
    const scalarField& cWeights
)
{
    FatalErrorIn
    (
        "labelList ptscotchDecomp::decompose"
        "("
            "const labelListList&, "
            "const pointField&, "
            "const scalarField&"
        ")"
    )   << notImplementedMessage << exit(FatalError);

    return labelList::null();
}


// ************************************************************************* //
