/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022 OpenFOAM Foundation
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

#include "polyCellSet.H"
#include "polyMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    template<> const char* NamedEnum
    <
        polyCellSet::selectionModeType,
        4
        >::names[] =
    {
        "points",
        "cellSet",
        "cellZone",
        "all"
    };

    const NamedEnum<polyCellSet::selectionModeType, 4>
        polyCellSet::selectionModeTypeNames_;
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::polyCellSet::setCells()
{
    Info<< incrIndent;

    switch (selectionMode_)
    {
        case selectionModeType::points:
        {
            Info<< indent << "- selecting cells using points" << endl;

            labelHashSet selectedCells;

            forAll(points_, i)
            {
                label celli = mesh_.findCell(points_[i]);
                if (celli >= 0)
                {
                    selectedCells.insert(celli);
                }

                label globalCelli = returnReduce(celli, maxOp<label>());
                if (globalCelli < 0)
                {
                    WarningInFunction
                        << "Unable to find owner cell for point " << points_[i]
                        << endl;
                }

            }

            cells_ = selectedCells.toc();

            break;
        }
        case selectionModeType::cellSet:
        {
            Info<< indent
                << "- selecting cells using cellSet " << cellSetName_ << endl;

            cells_ = cellSet(mesh_, cellSetName_).toc();

            break;
        }
        case selectionModeType::cellZone:
        {
            Info<< indent
                << "- selecting cells using cellZone " << cellSetName_ << endl;

            label zoneID = mesh_.cellZones().findZoneID(cellSetName_);
            if (zoneID == -1)
            {
                FatalErrorInFunction
                    << "Cannot find cellZone " << cellSetName_ << endl
                    << "Valid cellZones are " << mesh_.cellZones().names()
                    << exit(FatalError);
            }
            cells_ = mesh_.cellZones()[zoneID];

            break;
        }
        case selectionModeType::all:
        {
            Info<< indent << "- selecting all cells" << endl;
            cells_ = identity(mesh_.nCells());

            break;
        }
    }

    Info<< decrIndent;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::polyCellSet::polyCellSet
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    mesh_(mesh),
    selectionMode_(selectionModeType::all),
    cellSetName_(word::null)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::polyCellSet::~polyCellSet()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::polyCellSet::movePoints()
{
    if (selectionMode_ == selectionModeType::points)
    {
        setCells();
    }
}


void Foam::polyCellSet::topoChange(const polyTopoChangeMap&)
{
    setCells();
}


void Foam::polyCellSet::mapMesh(const polyMeshMap&)
{
    setCells();
}


void Foam::polyCellSet::distribute(const polyDistributionMap&)
{
    setCells();
}


bool Foam::polyCellSet::read(const dictionary& dict)
{
    if (dict.found("selectionMode"))
    {
        selectionMode_ =
            selectionModeTypeNames_.read(dict.lookup("selectionMode"));
    }
    else if (dict.found("points"))
    {
        selectionMode_ = selectionModeType::points;
    }
    else if (dict.found("cellSet"))
    {
        selectionMode_ = selectionModeType::cellSet;
    }
    else if (dict.found("cellZone"))
    {
        selectionMode_ = selectionModeType::cellZone;
    }
    else
    {
        FatalIOErrorInFunction(dict)
            << "selectionMode, cellZone, cellSet or points not specified."
            << nl << "Valid selectionModes:"
            << selectionModeTypeNames_.sortedToc()
            << exit(FatalIOError);
    }

    switch (selectionMode_)
    {
        case selectionModeType::points:
        {
            dict.lookup("points") >> points_;
            break;
        }
        case selectionModeType::cellSet:
        {
            dict.lookup("cellSet") >> cellSetName_;
            break;
        }
        case selectionModeType::cellZone:
        {
            dict.lookup("cellZone") >> cellSetName_;
            break;
        }
        case selectionModeType::all:
        {
            break;
        }
        default:
        {
            FatalErrorInFunction
                << "Unknown selectionMode "
                << selectionModeTypeNames_[selectionMode_]
                << nl << "Valid selectionModes:"
                << selectionModeTypeNames_.toc()
                << exit(FatalError);
        }
    }

    setCells();

    return true;
}


// ************************************************************************* //
