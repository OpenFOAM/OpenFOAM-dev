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

#include "fvCellSet.H"
#include "volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    template<> const char* NamedEnum
    <
        fvCellSet::selectionModeType,
        4
        >::names[] =
    {
        "points",
        "cellSet",
        "cellZone",
        "all"
    };

    const NamedEnum<fvCellSet::selectionModeType, 4>
        fvCellSet::selectionModeTypeNames_;
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fvCellSet::setCellSet()
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

            cellSet selectedCells(mesh_, cellSetName_);
            cells_ = selectedCells.toc();

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

    // Set volume information
    V_ = 0;
    forAll(cells_, i)
    {
        V_ += mesh_.V()[cells_[i]];
    }
    reduce(V_, sumOp<scalar>());

    Info<< indent
        << "- selected " << returnReduce(cells_.size(), sumOp<label>())
        << " cell(s) with volume " << V_ << endl;

    Info<< decrIndent;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvCellSet::fvCellSet
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    mesh_(mesh),
    selectionMode_(selectionModeType::all),
    cellSetName_(word::null),
    V_(NaN)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fvCellSet::~fvCellSet()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fvCellSet::updateMesh(const mapPolyMesh&)
{
    setCellSet();
}


bool Foam::fvCellSet::read(const dictionary& dict)
{
    selectionMode_ =
        selectionModeTypeNames_.read(dict.lookup("selectionMode"));

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
                << ". Valid selectionMode types are" << selectionModeTypeNames_
                << exit(FatalError);
        }
    }

    setCellSet();

    return true;
}


// ************************************************************************* //
