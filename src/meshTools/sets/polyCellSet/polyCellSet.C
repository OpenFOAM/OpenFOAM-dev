/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2024 OpenFOAM Foundation
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
        polyCellSet::selectionTypes,
        4
        >::names[] =
    {
        "points",
        "cellSet",
        "cellZone",
        "all"
    };

    const NamedEnum<polyCellSet::selectionTypes, 4>
        polyCellSet::selectionTypeNames;
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::polyCellSet::setCells()
{
    Info<< incrIndent;

    switch (selectionType_)
    {
        case selectionTypes::points:
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

                const label globalCelli = returnReduce(celli, maxOp<label>());
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
        case selectionTypes::cellSet:
        {
            Info<< indent
                << "- selecting cells using cellSet " << cellSetName_ << endl;

            cells_ = cellSet(mesh_, cellSetName_).toc();

            break;
        }
        case selectionTypes::cellZone:
        {
            Info<< indent
                << "- selecting cells using cellZone " << cellSetName_ << endl;

            const label zoneID = mesh_.cellZones().findIndex(cellSetName_);
            if (zoneID == -1)
            {
                FatalErrorInFunction
                    << "Cannot find cellZone " << cellSetName_ << endl
                    << "Valid cellZones are " << mesh_.cellZones().toc()
                    << exit(FatalError);
            }

            // cells_ not required for cellZone

            break;
        }
        case selectionTypes::all:
        {
            Info<< indent << "- selecting all cells" << endl;

            // cells_ not required for cellZone

            break;
        }
    }

    Info<< decrIndent;
}


Foam::labelUList Foam::polyCellSet::identityMap(const label len) const
{
    // Static identity list, resized as required
    static labelList map;

    if (len > map.size())
    {
        map.resize(len);

        forAll(map, i)
        {
            map[i] = i;
        }
    }

    return SubList<label>(map, len);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::polyCellSet::polyCellSet(const polyMesh& mesh)
:
    mesh_(mesh),
    selectionType_(selectionTypes::all),
    cellSetName_(word::null)
{}


Foam::polyCellSet::polyCellSet(const polyMesh& mesh, const dictionary& dict)
:
    mesh_(mesh),
    selectionType_(selectionTypes::all),
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
    if (selectionType_ == selectionTypes::points)
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
    if (dict.found("select") || dict.found("selectionMode"))
    {
        selectionType_ = selectionTypeNames.read
        (
            dict.lookupBackwardsCompatible({"select", "selectionMode"})
        );
    }
    else if (dict.found("points"))
    {
        selectionType_ = selectionTypes::points;
    }
    else if (dict.found("cellSet"))
    {
        selectionType_ = selectionTypes::cellSet;
    }
    else if (dict.found("cellZone"))
    {
        selectionType_ = selectionTypes::cellZone;
    }
    else
    {
        dict.lookup("select");
    }

    switch (selectionType_)
    {
        case selectionTypes::points:
        {
            dict.lookup("points") >> points_;
            break;
        }
        case selectionTypes::cellSet:
        {
            dict.lookup("cellSet") >> cellSetName_;
            break;
        }
        case selectionTypes::cellZone:
        {
            dict.lookup("cellZone") >> cellSetName_;
            break;
        }
        case selectionTypes::all:
        {
            break;
        }
        default:
        {
            FatalErrorInFunction
                << "Unknown selection type "
                << selectionTypeNames[selectionType_]
                << nl << "Valid selection type:"
                << selectionTypeNames.toc()
                << exit(FatalError);
        }
    }

    setCells();

    return true;
}


// ************************************************************************* //
