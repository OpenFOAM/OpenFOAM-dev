/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "sampledIsoSurface.H"
#include "dictionary.H"
#include "volFields.H"
#include "volPointInterpolation.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(sampledIsoSurface, 0);
    addNamedToRunTimeSelectionTable
    (
        sampledSurface,
        sampledIsoSurface,
        word,
        isoSurface
    );
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::sampledIsoSurface::getIsoFields() const
{
    const fvMesh& fvm = static_cast<const fvMesh&>(mesh());

    // Get volField
    // ~~~~~~~~~~~~

    if (fvm.foundObject<volScalarField>(isoField_))
    {
        if (debug)
        {
            InfoInFunction
                << "Lookup volField " << isoField_ << endl;
        }
        storedVolFieldPtr_.clear();
        volFieldPtr_ = &fvm.lookupObject<volScalarField>(isoField_);
    }
    else
    {
        // Bit of a hack. Read field and store.

        if (debug)
        {
            InfoInFunction
                << "Checking " << isoField_
                << " for same time " << fvm.time().timeName()
                << endl;
        }

        if
        (
            storedVolFieldPtr_.empty()
         || (fvm.time().timeName() != storedVolFieldPtr_().instance())
        )
        {
            if (debug)
            {
                InfoInFunction
                    << "Reading volField " << isoField_
                    << " from time " << fvm.time().timeName()
                    << endl;
            }

            IOobject vfHeader
            (
                isoField_,
                fvm.time().timeName(),
                fvm,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            );

            if (vfHeader.typeHeaderOk<volScalarField>(true))
            {
                storedVolFieldPtr_.reset
                (
                    new volScalarField
                    (
                        vfHeader,
                        fvm
                    )
                );
                volFieldPtr_ = storedVolFieldPtr_.operator->();
            }
            else
            {
                FatalErrorInFunction
                    << "Cannot find isosurface field " << isoField_
                    << " in database or directory " << vfHeader.path()
                    << exit(FatalError);
            }
        }
    }


    // Get pointField
    // ~~~~~~~~~~~~~~

    // In case of multiple iso values we don't want to calculate multiple e.g.
    // "volPointInterpolate(p)" so register it and re-use it. This is the
    // same as the 'cache' functionality from volPointInterpolate but
    // unfortunately that one does not guarantee that the field pointer
    // remain: e.g. some other functionObject might delete the cached version.
    // (volPointInterpolation::interpolate with cache=false deletes any
    //  registered one or if mesh.changing())

    if (!subMeshPtr_.valid())
    {
        const word pointFldName =
            "volPointInterpolate_"
          + type()
          + "("
          + isoField_
          + ')';

        if (fvm.foundObject<pointScalarField>(pointFldName))
        {
            if (debug)
            {
                InfoInFunction
                    << "lookup pointField " << pointFldName << endl;
            }
            const pointScalarField& pfld = fvm.lookupObject<pointScalarField>
            (
                pointFldName
            );

            if (!pfld.upToDate(*volFieldPtr_))
            {
                if (debug)
                {
                    InfoInFunction
                        << "updating pointField "
                        << pointFldName << endl;
                }
                // Update the interpolated value
                volPointInterpolation::New(fvm).interpolate
                (
                    *volFieldPtr_,
                    const_cast<pointScalarField&>(pfld)
                );
            }

            pointFieldPtr_ = &pfld;
        }
        else
        {
            // Not in registry. Interpolate.

            if (debug)
            {
                InfoInFunction
                    << "Checking pointField " << pointFldName
                    << " for same time " << fvm.time().timeName()
                    << endl;
            }

            // Interpolate without cache. Note that we're registering it
            // below so next time round it goes into the condition
            // above.
            tmp<pointScalarField> tpfld
            (
                volPointInterpolation::New(fvm).interpolate
                (
                    *volFieldPtr_,
                    pointFldName,
                    false
                )
            );
            pointFieldPtr_ = tpfld.ptr();
            const_cast<pointScalarField*>(pointFieldPtr_)->store();
        }


        // If averaging redo the volField. Can only be done now since needs the
        // point field.
        if (average_)
        {
            storedVolFieldPtr_.reset
            (
                pointAverage(*pointFieldPtr_).ptr()
            );
            volFieldPtr_ = storedVolFieldPtr_.operator->();
        }


        if (debug)
        {
            InfoInFunction
                << "volField " << volFieldPtr_->name()
                << " min:" << min(*volFieldPtr_).value()
                << " max:" << max(*volFieldPtr_).value() << endl;
            InfoInFunction
                << "pointField " << pointFieldPtr_->name()
                << " min:" << gMin(pointFieldPtr_->primitiveField())
                << " max:" << gMax(pointFieldPtr_->primitiveField()) << endl;
        }
    }
    else
    {
        // Get subMesh variants
        const fvMesh& subFvm = subMeshPtr_().subMesh();

        // Either lookup on the submesh or subset the whole-mesh volField

        if (subFvm.foundObject<volScalarField>(isoField_))
        {
            if (debug)
            {
                InfoInFunction
                    << "Sub-mesh lookup volField "
                    << isoField_ << endl;
            }
            storedVolSubFieldPtr_.clear();
            volSubFieldPtr_ = &subFvm.lookupObject<volScalarField>(isoField_);
        }
        else
        {
            if (debug)
            {
                InfoInFunction
                    << "Sub-setting volField " << isoField_ << endl;
            }
            storedVolSubFieldPtr_.reset
            (
                subMeshPtr_().interpolate
                (
                    *volFieldPtr_
                ).ptr()
            );
            storedVolSubFieldPtr_->checkOut();
            volSubFieldPtr_ = storedVolSubFieldPtr_.operator->();
        }


        // Pointfield on submesh

        word pointFldName =
            "volPointInterpolate("
          + volSubFieldPtr_->name()
          + ')';

        if (subFvm.foundObject<pointScalarField>(pointFldName))
        {
            if (debug)
            {
                InfoInFunction
                    << "Sub-mesh lookup pointField " << pointFldName << endl;
            }
            storedPointSubFieldPtr_.clear();
            pointSubFieldPtr_ = &subFvm.lookupObject<pointScalarField>
            (
                pointFldName
            );
        }
        else
        {
            if (debug)
            {
                InfoInFunction
                    << "Interpolating submesh volField "
                    << volSubFieldPtr_->name()
                    << " to get submesh pointField " << pointFldName << endl;
            }
            storedPointSubFieldPtr_.reset
            (
                volPointInterpolation::New
                (
                    subFvm
                ).interpolate(*volSubFieldPtr_).ptr()
            );
            storedPointSubFieldPtr_->checkOut();
            pointSubFieldPtr_ = storedPointSubFieldPtr_.operator->();
        }


        // If averaging redo the volField. Can only be done now since needs the
        // point field.
        if (average_)
        {
            storedVolSubFieldPtr_.reset
            (
                pointAverage(*pointSubFieldPtr_).ptr()
            );
            volSubFieldPtr_ = storedVolSubFieldPtr_.operator->();
        }


        if (debug)
        {
            InfoInFunction
                << "volSubField "
                << volSubFieldPtr_->name()
                << " min:" << min(*volSubFieldPtr_).value()
                << " max:" << max(*volSubFieldPtr_).value() << endl;
            InfoInFunction
                << "pointSubField "
                << pointSubFieldPtr_->name()
                << " min:" << gMin(pointSubFieldPtr_->primitiveField())
                << " max:" << gMax(pointSubFieldPtr_->primitiveField()) << endl;
        }
    }
}


bool Foam::sampledIsoSurface::updateGeometry() const
{
    const fvMesh& fvm = static_cast<const fvMesh&>(mesh());

    // No update needed
    if (fvm.time().timeIndex() == prevTimeIndex_)
    {
        return false;
    }

    // Get any subMesh
    if (zoneID_.index() != -1 && !subMeshPtr_.valid())
    {
        const polyBoundaryMesh& patches = mesh().boundaryMesh();

        // Patch to put exposed internal faces into
        const label exposedPatchi = patches.findPatchID(exposedPatchName_);

        if (debug)
        {
            Info<< "Allocating subset of size "
                << mesh().cellZones()[zoneID_.index()].size()
                << " with exposed faces into patch "
                << patches[exposedPatchi].name() << endl;
        }

        subMeshPtr_.reset
        (
            new fvMeshSubset(fvm)
        );
        subMeshPtr_().setLargeCellSubset
        (
            labelHashSet(mesh().cellZones()[zoneID_.index()]),
            exposedPatchi
        );
    }


    prevTimeIndex_ = fvm.time().timeIndex();
    getIsoFields();

    // Clear any stored topo
    surfPtr_.clear();
    facesPtr_.clear();

    // Clear derived data
    clearGeom();

    if (subMeshPtr_.valid())
    {
        surfPtr_.reset
        (
            new isoSurface
            (
                *volSubFieldPtr_,
                *pointSubFieldPtr_,
                isoVal_,
                regularise_,
                mergeTol_
            )
        );
    }
    else
    {
        surfPtr_.reset
        (
            new isoSurface
            (
                *volFieldPtr_,
                *pointFieldPtr_,
                isoVal_,
                regularise_,
                mergeTol_
            )
        );
    }


    if (debug)
    {
        Pout<< "sampledIsoSurface::updateGeometry() : constructed iso:"
            << nl
            << "    regularise     : " << regularise_ << nl
            << "    average        : " << average_ << nl
            << "    isoField       : " << isoField_ << nl
            << "    isoValue       : " << isoVal_ << nl;
        if (subMeshPtr_.valid())
        {
            Pout<< "    zone size      : " << subMeshPtr_().subMesh().nCells()
                << nl;
        }
        Pout<< "    points         : " << points().size() << nl
            << "    tris           : " << surface().size() << nl
            << "    cut cells      : " << surface().meshCells().size()
            << endl;
    }

    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sampledIsoSurface::sampledIsoSurface
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    sampledSurface(name, mesh, dict),
    isoField_(dict.lookup("isoField")),
    isoVal_(readScalar(dict.lookup("isoValue"))),
    mergeTol_(dict.lookupOrDefault("mergeTol", 1e-6)),
    regularise_(dict.lookupOrDefault("regularise", true)),
    average_(dict.lookupOrDefault("average", false)),
    zoneID_(dict.lookupOrDefault("zone", word::null), mesh.cellZones()),
    exposedPatchName_(word::null),
    surfPtr_(nullptr),
    facesPtr_(nullptr),
    prevTimeIndex_(-1),
    storedVolFieldPtr_(nullptr),
    volFieldPtr_(nullptr),
    pointFieldPtr_(nullptr)
{
    if (!sampledSurface::interpolate())
    {
        FatalIOErrorInFunction
        (
            dict
        )   << "Non-interpolated iso surface not supported since triangles"
            << " span across cells." << exit(FatalIOError);
    }

    if (zoneID_.index() != -1)
    {
        dict.lookup("exposedPatchName") >> exposedPatchName_;

        if (mesh.boundaryMesh().findPatchID(exposedPatchName_) == -1)
        {
            FatalIOErrorInFunction
            (
                dict
            )   << "Cannot find patch " << exposedPatchName_
                << " in which to put exposed faces." << endl
                << "Valid patches are " << mesh.boundaryMesh().names()
                << exit(FatalIOError);
        }

        if (debug && zoneID_.index() != -1)
        {
            Info<< "Restricting to cellZone " << zoneID_.name()
                << " with exposed internal faces into patch "
                << exposedPatchName_ << endl;
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sampledIsoSurface::~sampledIsoSurface()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::sampledIsoSurface::needsUpdate() const
{
    const fvMesh& fvm = static_cast<const fvMesh&>(mesh());

    return fvm.time().timeIndex() != prevTimeIndex_;
}


bool Foam::sampledIsoSurface::expire()
{
    surfPtr_.clear();
    facesPtr_.clear();
    subMeshPtr_.clear();

    // Clear derived data
    clearGeom();

    // already marked as expired
    if (prevTimeIndex_ == -1)
    {
        return false;
    }

    // force update
    prevTimeIndex_ = -1;
    return true;
}


bool Foam::sampledIsoSurface::update()
{
    return updateGeometry();
}


Foam::tmp<Foam::scalarField> Foam::sampledIsoSurface::sample
(
    const volScalarField& vField
) const
{
    return sampleField(vField);
}


Foam::tmp<Foam::vectorField> Foam::sampledIsoSurface::sample
(
    const volVectorField& vField
) const
{
    return sampleField(vField);
}


Foam::tmp<Foam::sphericalTensorField> Foam::sampledIsoSurface::sample
(
    const volSphericalTensorField& vField
) const
{
    return sampleField(vField);
}


Foam::tmp<Foam::symmTensorField> Foam::sampledIsoSurface::sample
(
    const volSymmTensorField& vField
) const
{
    return sampleField(vField);
}


Foam::tmp<Foam::tensorField> Foam::sampledIsoSurface::sample
(
    const volTensorField& vField
) const
{
    return sampleField(vField);
}


Foam::tmp<Foam::scalarField> Foam::sampledIsoSurface::interpolate
(
    const interpolation<scalar>& interpolator
) const
{
    return interpolateField(interpolator);
}


Foam::tmp<Foam::vectorField> Foam::sampledIsoSurface::interpolate
(
    const interpolation<vector>& interpolator
) const
{
    return interpolateField(interpolator);
}

Foam::tmp<Foam::sphericalTensorField> Foam::sampledIsoSurface::interpolate
(
    const interpolation<sphericalTensor>& interpolator
) const
{
    return interpolateField(interpolator);
}


Foam::tmp<Foam::symmTensorField> Foam::sampledIsoSurface::interpolate
(
    const interpolation<symmTensor>& interpolator
) const
{
    return interpolateField(interpolator);
}


Foam::tmp<Foam::tensorField> Foam::sampledIsoSurface::interpolate
(
    const interpolation<tensor>& interpolator
) const
{
    return interpolateField(interpolator);
}


void Foam::sampledIsoSurface::print(Ostream& os) const
{
    os  << "sampledIsoSurface: " << name() << " :"
        << "  field   :" << isoField_
        << "  value   :" << isoVal_;
}


// ************************************************************************* //
