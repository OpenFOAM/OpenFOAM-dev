/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2013 OpenFOAM Foundation
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

#include "surfaceCellSizeFunction.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(surfaceCellSizeFunction, 0);
    defineRunTimeSelectionTable(surfaceCellSizeFunction, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfaceCellSizeFunction::surfaceCellSizeFunction
(
    const word& type,
    const dictionary& surfaceCellSizeFunctionDict,
    const searchableSurface& surface,
    const scalar& defaultCellSize
)
:
    dictionary(surfaceCellSizeFunctionDict),
    surface_(surface),
    coeffsDict_(subDict(type + "Coeffs")),
    defaultCellSize_(defaultCellSize),
    refinementFactor_
    (
        lookupOrDefault<scalar>("refinementFactor", 1.0)
    )
{}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::surfaceCellSizeFunction> Foam::surfaceCellSizeFunction::New
(
    const dictionary& surfaceCellSizeFunctionDict,
    const searchableSurface& surface,
    const scalar& defaultCellSize
)
{
    word surfaceCellSizeFunctionTypeName
    (
        surfaceCellSizeFunctionDict.lookup("surfaceCellSizeFunction")
    );

    Info<< indent << "Selecting surfaceCellSizeFunction "
        << surfaceCellSizeFunctionTypeName << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(surfaceCellSizeFunctionTypeName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "surfaceCellSizeFunction::New(dictionary&, "
            "const conformalVoronoiMesh&, const searchableSurface&)"
        )   << "Unknown surfaceCellSizeFunction type "
            << surfaceCellSizeFunctionTypeName
            << endl << endl
            << "Valid surfaceCellSizeFunction types are :" << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<surfaceCellSizeFunction>
    (
        cstrIter()(surfaceCellSizeFunctionDict, surface, defaultCellSize)
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::surfaceCellSizeFunction::~surfaceCellSizeFunction()
{}


// ************************************************************************* //
