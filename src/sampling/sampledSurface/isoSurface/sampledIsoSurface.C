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

#include "sampledIsoSurface.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace sampledSurfaces
{
    defineTypeNameAndDebug(isoSurface, 0);
    addToRunTimeSelectionTable(sampledSurface, isoSurface, word);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::sampledSurfaces::isoSurface::updateGeometry() const
{
    const fvMesh& fvm = static_cast<const fvMesh&>(mesh());

    // no update needed
    if (fvm.time().timeIndex() == prevTimeIndex_)
    {
        return false;
    }

    prevTimeIndex_ = fvm.time().timeIndex();

    // Clear derived data
    sampledSurface::clearGeom();

    // Optionally read volScalarField
    autoPtr<volScalarField> readFieldPtr_;

    // 1. see if field in database
    // 2. see if field can be read
    const volScalarField* cellFldPtr = nullptr;
    if (fvm.foundObject<volScalarField>(isoField_))
    {
        if (debug)
        {
            InfoInFunction << "Lookup " << isoField_ << endl;
        }

        cellFldPtr = &fvm.lookupObject<volScalarField>(isoField_);
    }
    else
    {
        // Bit of a hack. Read field and store.

        if (debug)
        {
            InfoInFunction
                << "Reading " << isoField_
                << " from time " <<fvm.time().timeName()
                << endl;
        }

        readFieldPtr_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    isoField_,
                    fvm.time().timeName(),
                    fvm,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                fvm
            )
        );

        cellFldPtr = readFieldPtr_.operator->();
    }
    const volScalarField& cellFld = *cellFldPtr;

    tmp<pointScalarField> pointFld
    (
        volPointInterpolation::New(fvm).interpolate(cellFld)
    );

    PtrList<Foam::isoSurface> isos(isoVals_.size());
    forAll(isos, isoi)
    {
        isos.set
        (
            isoi,
            new Foam::isoSurface
            (
                fvm,
                cellFld.primitiveField(),
                pointFld().primitiveField(),
                isoVals_[isoi],
                filter_
            )
        );
    }

    if (isos.size() == 1)
    {
        // Straight transfer
        const_cast<isoSurface&>
        (
            *this
        ).MeshedSurface<face>::transfer(isos[0]);
        meshCells_ = isos[0].meshCells();
    }
    else
    {
        label nFaces = 0;
        label nPoints = 0;
        forAll(isos, isoi)
        {
            nFaces += isos[isoi].size();
            nPoints += isos[isoi].points().size();
        }

        faceList allFaces(nFaces);
        labelList allMeshCells(nFaces);
        pointField allPoints(nPoints);

        nFaces = 0;
        nPoints = 0;
        forAll(isos, isoi)
        {
            Foam::isoSurface& iso = isos[isoi];

            SubList<face> subAll(allFaces, iso.size(), nFaces);
            subAll = iso;
            // Offset point indices
            if (nPoints > 0)
            {
                forAll(subAll, i)
                {
                    face& f = subAll[i];
                    forAll(f, fp)
                    {
                        f[fp] += nPoints;
                    }
                }
            }
            SubList<label>(allMeshCells, iso.size(), nFaces) = iso.meshCells();
            nFaces += iso.size();

            const pointField& pts = iso.points();
            SubList<point>(allPoints, pts.size(), nPoints) = pts;
            nPoints += pts.size();

            // Clear out asap
            isos.set(isoi, nullptr);
        }

        if (nFaces != allFaces.size() || nPoints != allPoints.size())
        {
            FatalErrorInFunction << "nFaces:" << nFaces
                << " nPoints:" << nPoints << exit(FatalError);
        }


        surfZoneList allZones(1);
        allZones[0] = surfZone
        (
            "allFaces",
            allFaces.size(),    // size
            0,                  // start
            0                   // index
        );

        // Transfer
        const_cast<isoSurface&>
        (
            *this
        ).MeshedSurface<face>::reset
        (
            move(allPoints),
            move(allFaces),
            move(allZones)
        );
        meshCells_.transfer(allMeshCells);
    }
    if (debug)
    {
        Pout<< "sampledSurfaces::isoSurface::updateGeometry() : "
               "constructed iso:"
            << nl
            << "    filtering      : "
            << Foam::isoSurface::filterTypeNames_[filter_] << nl
            << "    isoField       : " << isoField_ << nl;
        if (isoVals_.size() == 1)
        {
            Pout<< "    isoValue       : " << isoVals_[0] << nl;
        }
        else
        {
            Pout<< "    isoValues      : " << isoVals_ << nl;
        }
        Pout<< "    points         : " << points().size() << nl
            << "    faces          : " << faces().size() << nl
            << "    cut cells      : " << meshCells_.size() << endl;
    }

    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sampledSurfaces::isoSurface::isoSurface
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    sampledSurface(name, mesh, dict),
    isoField_(dict.lookup("isoField")),
    isoVals_
    (
        dict.found("isoValues")
      ? scalarField(dict.lookup("isoValues"))
      : scalarField(1, dict.lookup<scalar>("isoValue"))
    ),
    filter_
    (
        dict.found("filtering")
      ? Foam::isoSurface::filterTypeNames_.read(dict.lookup("filtering"))
      : Foam::isoSurface::filterType::full
    ),
    prevTimeIndex_(-1),
    meshCells_(0)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sampledSurfaces::isoSurface::~isoSurface()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::sampledSurfaces::isoSurface::needsUpdate() const
{
    const fvMesh& fvm = static_cast<const fvMesh&>(mesh());

    return fvm.time().timeIndex() != prevTimeIndex_;
}


bool Foam::sampledSurfaces::isoSurface::expire()
{
    // Clear derived data
    sampledSurface::clearGeom();
    MeshedSurface<face>::clearGeom();

    // already marked as expired
    if (prevTimeIndex_ == -1)
    {
        return false;
    }

    // force update
    prevTimeIndex_ = -1;
    return true;
}


bool Foam::sampledSurfaces::isoSurface::update()
{
    return updateGeometry();
}


Foam::tmp<Foam::scalarField>
Foam::sampledSurfaces::isoSurface::sample
(
    const volScalarField& vField
) const
{
    return sampleField(vField);
}


Foam::tmp<Foam::vectorField>
Foam::sampledSurfaces::isoSurface::sample
(
    const volVectorField& vField
) const
{
    return sampleField(vField);
}


Foam::tmp<Foam::sphericalTensorField>
Foam::sampledSurfaces::isoSurface::sample
(
    const volSphericalTensorField& vField
) const
{
    return sampleField(vField);
}


Foam::tmp<Foam::symmTensorField>
Foam::sampledSurfaces::isoSurface::sample
(
    const volSymmTensorField& vField
) const
{
    return sampleField(vField);
}


Foam::tmp<Foam::tensorField>
Foam::sampledSurfaces::isoSurface::sample
(
    const volTensorField& vField
) const
{
    return sampleField(vField);
}


Foam::tmp<Foam::scalarField>
Foam::sampledSurfaces::isoSurface::interpolate
(
    const interpolation<scalar>& interpolator
) const
{
    return interpolateField(interpolator);
}


Foam::tmp<Foam::vectorField>
Foam::sampledSurfaces::isoSurface::interpolate
(
    const interpolation<vector>& interpolator
) const
{
    return interpolateField(interpolator);
}

Foam::tmp<Foam::sphericalTensorField>
Foam::sampledSurfaces::isoSurface::interpolate
(
    const interpolation<sphericalTensor>& interpolator
) const
{
    return interpolateField(interpolator);
}


Foam::tmp<Foam::symmTensorField>
Foam::sampledSurfaces::isoSurface::interpolate
(
    const interpolation<symmTensor>& interpolator
) const
{
    return interpolateField(interpolator);
}


Foam::tmp<Foam::tensorField>
Foam::sampledSurfaces::isoSurface::interpolate
(
    const interpolation<tensor>& interpolator
) const
{
    return interpolateField(interpolator);
}


void Foam::sampledSurfaces::isoSurface::print(Ostream& os) const
{
    os  << "isoSurface: " << name() << " :"
        << "  field:" << isoField_;
    if (isoVals_.size() == 1)
    {
        os << "  value:" << isoVals_[0];
    }
    else
    {
        os << "  values:" << isoVals_;
    }
    os  << "  faces:" << faces().size()
        << "  points:" << points().size();
}


// ************************************************************************* //
