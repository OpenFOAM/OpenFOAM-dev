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

#include "cellSetConstraint.H"
#include "volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace fv
    {
        defineTypeNameAndDebug(cellSetConstraint, 0);
    }

    template<> const char* NamedEnum
    <
        fv::cellSetConstraint::selectionModeType,
        4
        >::names[] =
    {
        "points",
        "cellSet",
        "cellZone",
        "all"
    };

    const NamedEnum<fv::cellSetConstraint::selectionModeType, 4>
        fv::cellSetConstraint::selectionModeTypeNames_;
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::cellSetConstraint::readCoeffs()
{
    selectionMode_ =
        selectionModeTypeNames_.read(coeffs().lookup("selectionMode"));

    switch (selectionMode_)
    {
        case selectionModeType::points:
        {
            coeffs().lookup("points") >> points_;
            break;
        }
        case selectionModeType::cellSet:
        {
            coeffs().lookup("cellSet") >> cellSetName_;
            break;
        }
        case selectionModeType::cellZone:
        {
            coeffs().lookup("cellZone") >> cellSetName_;
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
}


void Foam::fv::cellSetConstraint::setCellSet()
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
                label celli = mesh().findCell(points_[i]);
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

            cellSet selectedCells(mesh(), cellSetName_);
            cells_ = selectedCells.toc();

            break;
        }
        case selectionModeType::cellZone:
        {
            Info<< indent
                << "- selecting cells using cellZone " << cellSetName_ << endl;

            label zoneID = mesh().cellZones().findZoneID(cellSetName_);
            if (zoneID == -1)
            {
                FatalErrorInFunction
                    << "Cannot find cellZone " << cellSetName_ << endl
                    << "Valid cellZones are " << mesh().cellZones().names()
                    << exit(FatalError);
            }
            cells_ = mesh().cellZones()[zoneID];

            break;
        }
        case selectionModeType::all:
        {
            Info<< indent << "- selecting all cells" << endl;
            cells_ = identity(mesh().nCells());

            break;
        }
    }

    // Set volume information
    V_ = 0;
    forAll(cells_, i)
    {
        V_ += mesh().V()[cells_[i]];
    }
    reduce(V_, sumOp<scalar>());

    Info<< indent
        << "- selected " << returnReduce(cells_.size(), sumOp<label>())
        << " cell(s) with volume " << V_ << endl;

    Info<< decrIndent;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::cellSetConstraint::cellSetConstraint
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    fvConstraint(name, modelType, dict, mesh),
    selectionMode_(selectionModeType::all),
    cellSetName_(word::null),
    V_(NaN)
{
    readCoeffs();
    setCellSet();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fv::cellSetConstraint::~cellSetConstraint()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::cellSetConstraint::updateMesh(const mapPolyMesh&)
{
    setCellSet();
}


bool Foam::fv::cellSetConstraint::movePoints()
{
    return true;
}


bool Foam::fv::cellSetConstraint::read(const dictionary& dict)
{
    if (fvConstraint::read(dict))
    {
        readCoeffs();
        setCellSet();
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
