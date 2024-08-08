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

#include "ptscotch.H"
#include "addToRunTimeSelectionTable.H"
#include "Time.H"

static const char* notImplementedMessage =
"You are trying to use ptscotch but do not have the "
"ptscotch library loaded."
"\nThis message is from the dummy ptscotch stub library instead.\n"
"\n"
"Please install ptscotch and make sure that libptscotch.so is in your "
"LD_LIBRARY_PATH.\n"
"The ptscotch library can then be built in "
"$FOAM_SRC/parallel/decompose/ptscotch\n";


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace decompositionMethods
{
    defineTypeNameAndDebug(ptscotch, 0);

    addToRunTimeSelectionTable
    (
        decompositionMethod,
        ptscotch,
        distributor
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::decompositionMethods::ptscotch::check
(
    const int retVal,
    const char* str
)
{}


Foam::label Foam::decompositionMethods::ptscotch::decompose
(
    const fileName& meshPath,
    const List<label>& initxadj,
    const List<label>& initadjncy,
    const scalarField& initcWeights,

    List<label>& finalDecomp
) const
{
    FatalErrorInFunction
        << notImplementedMessage << exit(FatalError);

    return -1;
}


Foam::label Foam::decompositionMethods::ptscotch::decompose
(
    const fileName& meshPath,
    const label adjncySize,
    const label adjncy[],
    const label xadjSize,
    const label xadj[],
    const scalarField& cWeights,
    List<label>& finalDecomp
) const
{
    FatalErrorInFunction
        << notImplementedMessage << exit(FatalError);

    return -1;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::decompositionMethods::ptscotch::ptscotch
(
    const dictionary& decompositionDict,
    const dictionary& methodDict
)
:
    decompositionMethod(decompositionDict)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::labelList Foam::decompositionMethods::ptscotch::decompose
(
    const polyMesh& mesh,
    const pointField& points,
    const scalarField& pointWeights
)
{
    FatalErrorInFunction
        << notImplementedMessage << exit(FatalError);

    return labelList::null();
}


Foam::labelList Foam::decompositionMethods::ptscotch::decompose
(
    const polyMesh& mesh,
    const labelList& agglom,
    const pointField& agglomPoints,
    const scalarField& pointWeights
)
{
    FatalErrorInFunction
        << notImplementedMessage << exit(FatalError);

    return labelList::null();
}


Foam::labelList Foam::decompositionMethods::ptscotch::decompose
(
    const labelListList& globalCellCells,
    const pointField& cellCentres,
    const scalarField& cWeights
)
{
    FatalErrorInFunction
        << notImplementedMessage << exit(FatalError);

    return labelList::null();
}


// ************************************************************************* //
