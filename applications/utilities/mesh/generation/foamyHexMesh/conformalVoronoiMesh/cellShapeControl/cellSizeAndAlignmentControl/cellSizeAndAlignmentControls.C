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

#include "cellSizeAndAlignmentControls.H"
#include "searchableSurfaceControl.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(cellSizeAndAlignmentControls, 0);
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

bool Foam::cellSizeAndAlignmentControls::evalCellSizeFunctions
(
    const point& pt,
    scalar& minSize,
    label& maxPriority
) const
{
    bool anyFunctionFound = false;

    // Regions requesting with the same priority take the smallest

    if (controlFunctions_.size())
    {
        // Maintain priority of current hit. Initialise so it always goes
        // through at least once.
        label previousPriority = labelMin;

        forAll(controlFunctions_, i)
        {
            const cellSizeAndAlignmentControl& cSF = controlFunctions_[i];

            if (isA<searchableSurfaceControl>(cSF))
            {
                const searchableSurfaceControl& sSC =
                    refCast<const searchableSurfaceControl>(cSF);

                anyFunctionFound = sSC.cellSize(pt, minSize, previousPriority);

                if (previousPriority > maxPriority)
                {
                    maxPriority = previousPriority;
                }
            }
        }
    }

    return anyFunctionFound;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cellSizeAndAlignmentControls::cellSizeAndAlignmentControls
(
    const Time& runTime,
    const dictionary& shapeControlDict,
    const conformationSurfaces& geometryToConformTo,
    const scalar& defaultCellSize
)
:
    shapeControlDict_(shapeControlDict),
    geometryToConformTo_(geometryToConformTo),
    controlFunctions_(shapeControlDict_.size()),
    defaultCellSize_(defaultCellSize)
{
    label functionI = 0;

    forAllConstIter(dictionary, shapeControlDict_, iter)
    {
        word shapeControlEntryName = iter().keyword();

        const dictionary& controlFunctionDict
        (
            shapeControlDict_.subDict(shapeControlEntryName)
        );

        Info<< nl << "Shape Control : " << shapeControlEntryName << endl;
        Info<< incrIndent;

        controlFunctions_.set
        (
            functionI,
            cellSizeAndAlignmentControl::New
            (
                runTime,
                shapeControlEntryName,
                controlFunctionDict,
                geometryToConformTo_,
                defaultCellSize_
            )
        );

        Info<< decrIndent;

        functionI++;
    }

    // Sort controlFunctions_ by maxPriority
    SortableList<label> functionPriorities(functionI);

    forAll(controlFunctions_, funcI)
    {
        functionPriorities[funcI] = controlFunctions_[funcI].maxPriority();
    }

    functionPriorities.reverseSort();

    labelList invertedFunctionPriorities =
        invert(functionPriorities.size(), functionPriorities.indices());

    controlFunctions_.reorder(invertedFunctionPriorities);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::cellSizeAndAlignmentControls::~cellSizeAndAlignmentControls()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::scalar Foam::cellSizeAndAlignmentControls::cellSize
(
    const point& pt
) const
{
    scalar size = defaultCellSize_;
    label maxPriority = -1;

    evalCellSizeFunctions(pt, size, maxPriority);

    return size;
}


Foam::scalar Foam::cellSizeAndAlignmentControls::cellSize
(
    const point& pt,
    label& maxPriority
) const
{
    scalar size = defaultCellSize_;
    maxPriority = -1;

    evalCellSizeFunctions(pt, size, maxPriority);

    return size;
}


// ************************************************************************* //
