/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2019 OpenFOAM Foundation
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

#include "initialPointsMethod.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(initialPointsMethod, 0);
defineRunTimeSelectionTable(initialPointsMethod, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

initialPointsMethod::initialPointsMethod
(
    const word& type,
    const dictionary& initialPointsDict,
    const Time& runTime,
    Random& rndGen,
    const conformationSurfaces& geometryToConformTo,
    const cellShapeControl& cellShapeControls,
    const autoPtr<backgroundMeshDecomposition>& decomposition
)
:
    dictionary(initialPointsDict),
    runTime_(runTime),
    rndGen_(rndGen),
    geometryToConformTo_(geometryToConformTo),
    cellShapeControls_(cellShapeControls),
    decomposition_(decomposition),
    detailsDict_(optionalSubDict(type + "Coeffs")),
    minimumSurfaceDistanceCoeffSqr_
    (
        sqr( initialPointsDict.lookup<scalar>("minimumSurfaceDistanceCoeff"))
    ),
    fixInitialPoints_(Switch(initialPointsDict.lookup("fixInitialPoints")))
{}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

autoPtr<initialPointsMethod> initialPointsMethod::New
(
    const dictionary& initialPointsDict,
    const Time& runTime,
    Random& rndGen,
    const conformationSurfaces& geometryToConformTo,
    const cellShapeControl& cellShapeControls,
    const autoPtr<backgroundMeshDecomposition>& decomposition
)
{
    word initialPointsMethodTypeName
    (
        initialPointsDict.lookup("initialPointsMethod")
    );

    Info<< nl << "Selecting initialPointsMethod "
        << initialPointsMethodTypeName << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(initialPointsMethodTypeName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown initialPointsMethod type "
            << initialPointsMethodTypeName
            << endl << endl
            << "Valid initialPointsMethod types are :" << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return
        autoPtr<initialPointsMethod>
        (
            cstrIter()
            (
                initialPointsDict,
                runTime,
                rndGen,
                geometryToConformTo,
                cellShapeControls,
                decomposition
            )
        );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

initialPointsMethod::~initialPointsMethod()
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
