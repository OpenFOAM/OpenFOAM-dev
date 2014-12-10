/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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

#include "fvOption.H"
#include "fvMesh.H"
#include "fvMatrices.H"
#include "volFields.H"
#include "fvsPatchFields.H"
#include "ListOps.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace fv
    {
        defineTypeNameAndDebug(option, 0);
        defineRunTimeSelectionTable(option, dictionary);
    }

    template<> const char* NamedEnum
    <
        fv::option::selectionModeType,
        5
        >::names[] =
    {
        "points",
        "cellSet",
        "cellZone",
        "mapRegion",
        "all"
    };

    const NamedEnum<fv::option::selectionModeType, 5>
        fv::option::selectionModeTypeNames_;
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

bool Foam::fv::option::alwaysApply() const
{
    return false;
}


void Foam::fv::option::setSelection(const dictionary& dict)
{
    switch (selectionMode_)
    {
        case smPoints:
        {
            dict.lookup("points") >> points_;
            break;
        }
        case smCellSet:
        {
            dict.lookup("cellSet") >> cellSetName_;
            break;
        }
        case smCellZone:
        {
            dict.lookup("cellZone") >> cellSetName_;
            break;
        }
        case smMapRegion:
        {
            dict.lookup("nbrRegionName") >> nbrRegionName_;
            break;
        }
        case smAll:
        {
            break;
        }
        default:
        {
            FatalErrorIn("option::setSelection(const dictionary&)")
                << "Unknown selectionMode "
                << selectionModeTypeNames_[selectionMode_]
                << ". Valid selectionMode types are" << selectionModeTypeNames_
                << exit(FatalError);
        }
    }
}


void Foam::fv::option::setCellSet()
{
    switch (selectionMode_)
    {
        case smPoints:
        {
            Info<< indent << "- selecting cells using points" << endl;

            labelHashSet selectedCells;

            forAll(points_, i)
            {
                label cellI = mesh_.findCell(points_[i]);
                if (cellI >= 0)
                {
                    selectedCells.insert(cellI);
                }

                label globalCellI = returnReduce(cellI, maxOp<label>());
                if (globalCellI < 0)
                {
                    WarningIn("option::setCellSet()")
                        << "Unable to find owner cell for point " << points_[i]
                        << endl;
                }

            }

            cells_ = selectedCells.toc();

            break;
        }
        case smCellSet:
        {
            Info<< indent
                << "- selecting cells using cellSet " << cellSetName_ << endl;

            cellSet selectedCells(mesh_, cellSetName_);
            cells_ = selectedCells.toc();

            break;
        }
        case smCellZone:
        {
            Info<< indent
                << "- selecting cells using cellZone " << cellSetName_ << endl;

            label zoneID = mesh_.cellZones().findZoneID(cellSetName_);
            if (zoneID == -1)
            {
                FatalErrorIn("option::setCellIds()")
                    << "Cannot find cellZone " << cellSetName_ << endl
                    << "Valid cellZones are " << mesh_.cellZones().names()
                    << exit(FatalError);
            }
            cells_ = mesh_.cellZones()[zoneID];

            break;
        }
        case smMapRegion:
        {
            if (active_ && master_)
            {
                Info<< indent << "- selecting inter region mapping" << endl;

                const fvMesh& nbrMesh =
                    mesh_.time().lookupObject<fvMesh>(nbrRegionName_);

                if (mesh_.name() == nbrMesh.name())
                {
                    FatalErrorIn("option::setCellIds()")
                        << "Inter-region model selected, but local and "
                        << "neighbour regions are the same: " << nl
                        << "    local region: " << mesh_.name() << nl
                        << "    secondary region: " << nbrMesh.name() << nl
                        << exit(FatalError);
                }

                if (mesh_.bounds().overlaps(nbrMesh.bounds()))
                {
                    meshInterpPtr_.reset
                    (
                        new meshToMesh
                        (
                            mesh_,
                            nbrMesh,
                            meshToMesh::interpolationMethodNames_.read
                            (
                                dict_.lookup("interpolationMethod")
                            ),
                            false // not interpolating patches
                        )
                    );
                }
                else
                {
                    FatalErrorIn("option::setCellSet()")
                        << "regions " << mesh_.name() << " and "
                        << nbrMesh.name() <<  " do not intersect"
                        << exit(FatalError);
                }

                V_ = meshInterpPtr_->V();
            }
            break;
        }
        case smAll:
        {
            Info<< indent << "- selecting all cells" << endl;
            cells_ = identity(mesh_.nCells());

            break;
        }
        default:
        {
            FatalErrorIn("option::setCellSet()")
                << "Unknown selectionMode "
                << selectionModeTypeNames_[selectionMode_]
                << ". Valid selectionMode types are" << selectionModeTypeNames_
                << exit(FatalError);
        }
    }

    // Set volume information
    if (selectionMode_ != smMapRegion)
    {
        V_ = 0.0;
        forAll(cells_, i)
        {
            V_ += mesh_.V()[cells_[i]];
        }
        reduce(V_, sumOp<scalar>());

        Info<< indent
            << "- selected " << returnReduce(cells_.size(), sumOp<label>())
            << " cell(s) with volume " << V_ << nl << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::option::option
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh,
    const bool master
)
:
    name_(name),
    modelType_(modelType),
    mesh_(mesh),
    dict_(dict),
    coeffs_(dict.subDict(modelType + "Coeffs")),
    active_(readBool(dict_.lookup("active"))),
    timeStart_(-1.0),
    duration_(0.0),
    selectionMode_(selectionModeTypeNames_.read(dict_.lookup("selectionMode"))),
    cellSetName_("none"),
    V_(0.0),
    meshInterpPtr_(),
    nbrRegionName_("none"),
    master_(master),
    fieldNames_(),
    applied_()
{
    Info<< incrIndent << indent << "Source: " << name_ << endl;

    if (dict_.readIfPresent("timeStart", timeStart_))
    {
        dict_.lookup("duration") >> duration_;
        Info<< indent << "- applying source at time " << timeStart_
            << " for duration " << duration_ << endl;
    }
    else
    {
        Info<< indent << "- applying source for all time" << endl;
    }

    setSelection(dict_);

    setCellSet();

    Info<< decrIndent;
}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::fv::option> Foam::fv::option::New
(
    const word& name,
    const dictionary& coeffs,
    const fvMesh& mesh
)
{
    word modelType(coeffs.lookup("type"));

    Info<< indent
        << "Selecting finite volume options model type " << modelType << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(modelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "option::New(const word&, const dictionary&, const fvMesh&)"
        )   << "Unknown Model type " << modelType << nl << nl
            << "Valid model types are:" << nl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<option>(cstrIter()(name, modelType, coeffs, mesh));
}


Foam::fv::option::~option()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fv::option::isActive()
{
    if (active_ && inTimeLimits(mesh_.time().value()))
    {
        // Update the cell set if the mesh is changing
        if (mesh_.changing())
        {
            setCellSet();
        }
        return true;
    }
    else
    {
        return false;
    }
}


Foam::label Foam::fv::option::applyToField(const word& fieldName) const
{
    if (alwaysApply())
    {
        return 0;
    }

    return findIndex(fieldNames_, fieldName);
}


void Foam::fv::option::checkApplied() const
{
    forAll(applied_, i)
    {
        if (!applied_[i])
        {
            WarningIn("void option::checkApplied() const")
                << "Source " << name_ << " defined for field "
                << fieldNames_[i] << " but never used" << endl;
        }
    }
}


void Foam::fv::option::correct(volScalarField& fld)
{
    // do nothing
}


void Foam::fv::option::correct(volVectorField& fld)
{
    // do nothing
}


void Foam::fv::option::correct(volSphericalTensorField& fld)
{
    // do nothing
}


void Foam::fv::option::correct(volSymmTensorField& fld)
{
    // do nothing
}


void Foam::fv::option::correct(volTensorField& fld)
{
    // do nothing
}


void Foam::fv::option::addSup
(
    fvMatrix<scalar>& eqn,
    const label fieldI
)
{
    // do nothing
}


void Foam::fv::option::addSup
(
    fvMatrix<vector>& eqn,
    const label fieldI
)
{
    // do nothing
}


void Foam::fv::option::addSup
(
    fvMatrix<sphericalTensor>& eqn,
    const label fieldI
)
{
    // do nothing
}


void Foam::fv::option::addSup
(
    fvMatrix<symmTensor>& eqn,
    const label fieldI
)
{
    // do nothing
}


void Foam::fv::option::addSup
(
    fvMatrix<tensor>& eqn,
    const label fieldI
)
{
    // do nothing
}


void Foam::fv::option::addSup
(
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const label fieldI
)
{
    // do nothing
}


void Foam::fv::option::addSup
(
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const label fieldI
)
{
    // do nothing
}


void Foam::fv::option::addSup
(
    const volScalarField& rho,
    fvMatrix<sphericalTensor>& eqn,
    const label fieldI
)
{
    // do nothing
}


void Foam::fv::option::addSup
(
    const volScalarField& rho,
    fvMatrix<symmTensor>& eqn,
    const label fieldI
)
{
    // do nothing
}


void Foam::fv::option::addSup
(
    const volScalarField& rho,
    fvMatrix<tensor>& eqn,
    const label fieldI
)
{
    // do nothing
}


void Foam::fv::option::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const label fieldI
)
{
    addSup(alpha*rho, eqn, fieldI);
}


void Foam::fv::option::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const label fieldI
)
{
    addSup(alpha*rho, eqn, fieldI);
}


void Foam::fv::option::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    fvMatrix<sphericalTensor>& eqn,
    const label fieldI
)
{
    addSup(alpha*rho, eqn, fieldI);
}


void Foam::fv::option::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    fvMatrix<symmTensor>& eqn,
    const label fieldI
)
{
    addSup(alpha*rho, eqn, fieldI);
}


void Foam::fv::option::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    fvMatrix<tensor>& eqn,
    const label fieldI
)
{
    addSup(alpha*rho, eqn, fieldI);
}


void Foam::fv::option::setValue(fvMatrix<scalar>& eqn, const label fieldI)
{
    // do nothing
}


void Foam::fv::option::setValue(fvMatrix<vector>& eqn, const label fieldI)
{
    // do nothing
}


void Foam::fv::option::setValue
(
    fvMatrix<sphericalTensor>& eqn,
    const label fieldI
)
{
    // do nothing
}


void Foam::fv::option::setValue
(
    fvMatrix<symmTensor>& eqn,
    const label fieldI
)
{
    // do nothing
}


void Foam::fv::option::setValue(fvMatrix<tensor>& eqn, const label fieldI)
{
    // do nothing
}


void Foam::fv::option::makeRelative(surfaceScalarField& phi) const
{
    // do nothing
}


void Foam::fv::option::makeRelative
(
    FieldField<fvsPatchField, scalar>& phi
) const
{
    // do nothing
}


void Foam::fv::option::makeRelative
(
    const surfaceScalarField& rho,
    surfaceScalarField& phi
) const
{
    // do nothing
}


void Foam::fv::option::makeAbsolute(surfaceScalarField& phi) const
{
    // do nothing
}


void Foam::fv::option::makeAbsolute
(
    const surfaceScalarField& rho,
    surfaceScalarField& phi
) const
{
    // do nothing
}


// ************************************************************************* //
