/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

#include "distanceSurface.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace sampledSurfaces
{
    defineTypeNameAndDebug(distanceSurface, 0);
    addToRunTimeSelectionTable(sampledSurface, distanceSurface, word);
}
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::sampledSurfaces::distanceSurface::createGeometry()
{
    if (debug)
    {
        Pout<< "distanceSurface::createGeometry :updating geometry." << endl;
    }

    // Clear any stored topologies
    isoSurfPtr_.clear();

    // Clear derived data
    clearGeom();

    // Get any subMesh
    if (zoneKey_.size() && !subMeshPtr_.valid())
    {
        const polyBoundaryMesh& patches = mesh().boundaryMesh();

        // Patch to put exposed internal faces into
        const label exposedPatchi = patches.findPatchID(exposedPatchName_);

        subMeshPtr_.reset
        (
            new fvMeshSubset(static_cast<const fvMesh&>(mesh()))
        );
        subMeshPtr_().setLargeCellSubset
        (
            labelHashSet(mesh().cellZones().findMatching(zoneKey_).used()),
            exposedPatchi
        );

        DebugInfo
            << "Allocating subset of size " << subMeshPtr_().subMesh().nCells()
            << " with exposed faces into patch "
            << patches[exposedPatchi].name() << endl;
    }

    // Select either the submesh or the underlying mesh
    const fvMesh& mesh =
    (
        subMeshPtr_.valid()
      ? subMeshPtr_().subMesh()
      : static_cast<const fvMesh&>(this->mesh())
    );


    // Distance to cell centres
    // ~~~~~~~~~~~~~~~~~~~~~~~~

    cellDistancePtr_.reset
    (
        new volScalarField
        (
            IOobject
            (
                "cellDistance",
                mesh.time().timeName(),
                mesh.time(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh,
            dimensionedScalar(dimLength, 0)
        )
    );
    volScalarField& cellDistance = cellDistancePtr_();

    // Internal field
    {
        const pointField& cc = mesh.C();
        scalarField& fld = cellDistance.primitiveFieldRef();

        List<pointIndexHit> nearest;
        surfPtr_().findNearest
        (
            cc,
            scalarField(cc.size(), great),
            nearest
        );

        if (signed_)
        {
            List<volumeType> volType;

            surfPtr_().getVolumeType(cc, volType);

            forAll(volType, i)
            {
                volumeType vT = volType[i];

                if (vT == volumeType::outside)
                {
                    fld[i] = Foam::mag(cc[i] - nearest[i].hitPoint());
                }
                else if (vT == volumeType::inside)
                {
                    fld[i] = -Foam::mag(cc[i] - nearest[i].hitPoint());
                }
                else
                {
                    FatalErrorInFunction
                        << "getVolumeType failure, neither INSIDE or OUTSIDE"
                        << exit(FatalError);
                }
            }
        }
        else
        {
            forAll(nearest, i)
            {
                fld[i] = Foam::mag(cc[i] - nearest[i].hitPoint());
            }
        }
    }

    volScalarField::Boundary& cellDistanceBf =
        cellDistance.boundaryFieldRef();

    // Patch fields
    {
        forAll(mesh.C().boundaryField(), patchi)
        {
            const pointField& cc = mesh.C().boundaryField()[patchi];
            fvPatchScalarField& fld = cellDistanceBf[patchi];

            List<pointIndexHit> nearest;
            surfPtr_().findNearest
            (
                cc,
                scalarField(cc.size(), great),
                nearest
            );

            if (signed_)
            {
                List<volumeType> volType;

                surfPtr_().getVolumeType(cc, volType);

                forAll(volType, i)
                {
                    volumeType vT = volType[i];

                    if (vT == volumeType::outside)
                    {
                        fld[i] = Foam::mag(cc[i] - nearest[i].hitPoint());
                    }
                    else if (vT == volumeType::inside)
                    {
                        fld[i] = -Foam::mag(cc[i] - nearest[i].hitPoint());
                    }
                    else
                    {
                        FatalErrorInFunction
                            << "getVolumeType failure, "
                            << "neither INSIDE or OUTSIDE"
                            << exit(FatalError);
                    }
                }
            }
            else
            {
                forAll(nearest, i)
                {
                    fld[i] = Foam::mag(cc[i] - nearest[i].hitPoint());
                }
            }
        }
    }


    // On processor patches the mesh.C() will already be the cell centre
    // on the opposite side so no need to swap cellDistance.


    // Distance to points
    pointDistance_.setSize(mesh.nPoints());
    {
        const pointField& pts = mesh.points();

        List<pointIndexHit> nearest;
        surfPtr_().findNearest
        (
            pts,
            scalarField(pts.size(), great),
            nearest
        );

        if (signed_)
        {
            List<volumeType> volType;

            surfPtr_().getVolumeType(pts, volType);

            forAll(volType, i)
            {
                volumeType vT = volType[i];

                if (vT == volumeType::outside)
                {
                    pointDistance_[i] =
                        Foam::mag(pts[i] - nearest[i].hitPoint());
                }
                else if (vT == volumeType::inside)
                {
                    pointDistance_[i] =
                        -Foam::mag(pts[i] - nearest[i].hitPoint());
                }
                else
                {
                    FatalErrorInFunction
                        << "getVolumeType failure, neither INSIDE or OUTSIDE"
                        << exit(FatalError);
                }
            }
        }
        else
        {
            forAll(nearest, i)
            {
                pointDistance_[i] = Foam::mag(pts[i]-nearest[i].hitPoint());
            }
        }
    }


    if (debug)
    {
        Pout<< "Writing cell distance:" << cellDistance.objectPath() << endl;
        cellDistance.write();
        pointScalarField pDist
        (
            IOobject
            (
                "pointDistance",
                mesh.time().timeName(),
                mesh.time(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            pointMesh::New(mesh),
            dimensionedScalar(dimLength, 0)
        );
        pDist.primitiveFieldRef() = pointDistance_;

        Pout<< "Writing point distance:" << pDist.objectPath() << endl;
        pDist.write();
    }


    //- Direct from cell field and point field.
    isoSurfPtr_.reset
    (
        new isoSurface
        (
            mesh,
            cellDistance,
            pointDistance_,
            distance_,
            filter_
        )
    );

    if (debug)
    {
        print(Pout);
        Pout<< endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sampledSurfaces::distanceSurface::distanceSurface
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    sampledSurface(name, mesh, dict),
    surfPtr_
    (
        searchableSurface::New
        (
            dict.lookup("surfaceType"),
            IOobject
            (
                dict.lookupOrDefault("surfaceName", name),  // name
                mesh.time().constant(),                     // directory
                "triSurface",                               // instance
                mesh.time(),                                // registry
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            dict
        )
    ),
    distance_(dict.lookup<scalar>("distance")),
    signed_(readBool(dict.lookup("signed"))),
    filter_
    (
        dict.found("filtering")
      ? isoSurface::filterTypeNames_.read(dict.lookup("filtering"))
      : isoSurface::filterType::full
    ),
    average_(dict.lookupOrDefault("average", false)),
    zoneKey_(dict.lookupOrDefault("zone", wordRe::null)),
    needsUpdate_(true),
    subMeshPtr_(nullptr),
    isoSurfPtr_(nullptr)
{
    if (zoneKey_.size())
    {
        dict.lookup("exposedPatchName") >> exposedPatchName_;

        if (mesh.boundaryMesh().findPatchID(exposedPatchName_) == -1)
        {
            FatalErrorInFunction
                << "Cannot find patch " << exposedPatchName_
                << " in which to put exposed faces." << endl
                << "Valid patches are " << mesh.boundaryMesh().names()
                << exit(FatalError);
        }

        if (debug)
        {
            Info<< "Restricting to cellZone " << zoneKey_
                << " with exposed internal faces into patch "
                << exposedPatchName_ << endl;
        }

        if (mesh.cellZones().findIndex(zoneKey_) < 0)
        {
            WarningInFunction
                << "cellZone " << zoneKey_
                << " not found - using entire mesh" << endl;
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sampledSurfaces::distanceSurface::~distanceSurface()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::sampledSurfaces::distanceSurface::needsUpdate() const
{
    return needsUpdate_;
}


bool Foam::sampledSurfaces::distanceSurface::expire()
{
    if (debug)
    {
        Pout<< "distanceSurface::expire :"
            << " needsUpdate_:" << needsUpdate_ << endl;
    }

    // Clear derived data
    clearGeom();

    // already marked as expired
    if (needsUpdate_)
    {
        return false;
    }

    needsUpdate_ = true;
    return true;
}


bool Foam::sampledSurfaces::distanceSurface::update()
{
    if (debug)
    {
        Pout<< "distanceSurface::update :"
            << " needsUpdate_:" << needsUpdate_ << endl;
    }

    if (!needsUpdate_)
    {
        return false;
    }

    createGeometry();

    needsUpdate_ = false;
    return true;
}


Foam::tmp<Foam::scalarField>
Foam::sampledSurfaces::distanceSurface::sample
(
    const volScalarField& vField
) const
{
    return sampleField(vField);
}


Foam::tmp<Foam::vectorField>
Foam::sampledSurfaces::distanceSurface::sample
(
    const volVectorField& vField
) const
{
    return sampleField(vField);
}


Foam::tmp<Foam::sphericalTensorField>
Foam::sampledSurfaces::distanceSurface::sample
(
    const volSphericalTensorField& vField
) const
{
    return sampleField(vField);
}


Foam::tmp<Foam::symmTensorField>
Foam::sampledSurfaces::distanceSurface::sample
(
    const volSymmTensorField& vField
) const
{
    return sampleField(vField);
}


Foam::tmp<Foam::tensorField>
Foam::sampledSurfaces::distanceSurface::sample
(
    const volTensorField& vField
) const
{
    return sampleField(vField);
}


Foam::tmp<Foam::scalarField>
Foam::sampledSurfaces::distanceSurface::interpolate
(
    const interpolation<scalar>& interpolator
) const
{
    return interpolateField(interpolator);
}


Foam::tmp<Foam::vectorField>
Foam::sampledSurfaces::distanceSurface::interpolate
(
    const interpolation<vector>& interpolator
) const
{
    return interpolateField(interpolator);
}

Foam::tmp<Foam::sphericalTensorField>
Foam::sampledSurfaces::distanceSurface::interpolate
(
    const interpolation<sphericalTensor>& interpolator
) const
{
    return interpolateField(interpolator);
}


Foam::tmp<Foam::symmTensorField>
Foam::sampledSurfaces::distanceSurface::interpolate
(
    const interpolation<symmTensor>& interpolator
) const
{
    return interpolateField(interpolator);
}


Foam::tmp<Foam::tensorField>
Foam::sampledSurfaces::distanceSurface::interpolate
(
    const interpolation<tensor>& interpolator
) const
{
    return interpolateField(interpolator);
}


void Foam::sampledSurfaces::distanceSurface::print(Ostream& os) const
{
    os  << "distanceSurface: " << name() << " :"
        << "  surface:" << surfPtr_().name()
        << "  distance:" << distance_
        << "  faces:" << faces().size()
        << "  points:" << points().size();
}


// ************************************************************************* //
