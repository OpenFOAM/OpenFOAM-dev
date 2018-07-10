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

#include "distanceSurface.H"
#include "dictionary.H"
#include "volFields.H"
#include "volPointInterpolation.H"
#include "addToRunTimeSelectionTable.H"
#include "fvMesh.H"
#include "volumeType.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(distanceSurface, 0);
    addToRunTimeSelectionTable(sampledSurface, distanceSurface, word);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::distanceSurface::createGeometry()
{
    if (debug)
    {
        Pout<< "distanceSurface::createGeometry :updating geometry." << endl;
    }

    // Clear any stored topologies
    facesPtr_.clear();
    isoSurfCellPtr_.clear();
    isoSurfPtr_.clear();

    // Clear derived data
    clearGeom();

    const fvMesh& fvm = static_cast<const fvMesh&>(mesh());

    // Distance to cell centres
    // ~~~~~~~~~~~~~~~~~~~~~~~~

    cellDistancePtr_.reset
    (
        new volScalarField
        (
            IOobject
            (
                "cellDistance",
                fvm.time().timeName(),
                fvm.time(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            fvm,
            dimensionedScalar("zero", dimLength, 0)
        )
    );
    volScalarField& cellDistance = cellDistancePtr_();

    // Internal field
    {
        const pointField& cc = fvm.C();
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

                if (vT == volumeType::OUTSIDE)
                {
                    fld[i] = Foam::mag(cc[i] - nearest[i].hitPoint());
                }
                else if (vT == volumeType::INSIDE)
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
        forAll(fvm.C().boundaryField(), patchi)
        {
            const pointField& cc = fvm.C().boundaryField()[patchi];
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

                    if (vT == volumeType::OUTSIDE)
                    {
                        fld[i] = Foam::mag(cc[i] - nearest[i].hitPoint());
                    }
                    else if (vT == volumeType::INSIDE)
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
    pointDistance_.setSize(fvm.nPoints());
    {
        const pointField& pts = fvm.points();

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

                if (vT == volumeType::OUTSIDE)
                {
                    pointDistance_[i] =
                        Foam::mag(pts[i] - nearest[i].hitPoint());
                }
                else if (vT == volumeType::INSIDE)
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
                fvm.time().timeName(),
                fvm.time(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            pointMesh::New(fvm),
            dimensionedScalar("zero", dimLength, 0)
        );
        pDist.primitiveFieldRef() = pointDistance_;

        Pout<< "Writing point distance:" << pDist.objectPath() << endl;
        pDist.write();
    }


    //- Direct from cell field and point field.
    if (cell_)
    {
        isoSurfCellPtr_.reset
        (
            new isoSurfaceCell
            (
                fvm,
                cellDistance,
                pointDistance_,
                distance_,
                regularise_
            )
        );
    }
    else
    {
        isoSurfPtr_.reset
        (
            new isoSurface
            (
                cellDistance,
                pointDistance_,
                distance_,
                regularise_
            )
        );
    }

    if (debug)
    {
        print(Pout);
        Pout<< endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::distanceSurface::distanceSurface
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
    distance_(readScalar(dict.lookup("distance"))),
    signed_(readBool(dict.lookup("signed"))),
    cell_(dict.lookupOrDefault("cell", true)),
    regularise_(dict.lookupOrDefault("regularise", true)),
    average_(dict.lookupOrDefault("average", false)),
    zoneKey_(keyType::null),
    needsUpdate_(true),
    isoSurfCellPtr_(nullptr),
    isoSurfPtr_(nullptr),
    facesPtr_(nullptr)
{
//    dict.readIfPresent("zone", zoneKey_);
//
//    if (debug && zoneKey_.size() && mesh.cellZones().findZoneID(zoneKey_) < 0)
//    {
//        Info<< "cellZone " << zoneKey_
//            << " not found - using entire mesh" << endl;
//    }
}



Foam::distanceSurface::distanceSurface
(
    const word& name,
    const polyMesh& mesh,
    const bool interpolate,
    const word& surfaceType,
    const word& surfaceName,
    const scalar distance,
    const bool signedDistance,
    const bool cell,
    const Switch regularise,
    const Switch average
)
:
    sampledSurface(name, mesh, interpolate),
    surfPtr_
    (
        searchableSurface::New
        (
            surfaceType,
            IOobject
            (
                surfaceName,  // name
                mesh.time().constant(),                     // directory
                "triSurface",                               // instance
                mesh.time(),                                // registry
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            dictionary()
        )
    ),
    distance_(distance),
    signed_(signedDistance),
    cell_(cell),
    regularise_(regularise),
    average_(average),
    zoneKey_(keyType::null),
    needsUpdate_(true),
    isoSurfCellPtr_(nullptr),
    isoSurfPtr_(nullptr),
    facesPtr_(nullptr)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::distanceSurface::~distanceSurface()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::distanceSurface::needsUpdate() const
{
    return needsUpdate_;
}


bool Foam::distanceSurface::expire()
{
    if (debug)
    {
        Pout<< "distanceSurface::expire :"
            << " have-facesPtr_:" << facesPtr_.valid()
            << " needsUpdate_:" << needsUpdate_ << endl;
    }

    // Clear any stored topologies
    facesPtr_.clear();

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


bool Foam::distanceSurface::update()
{
    if (debug)
    {
        Pout<< "distanceSurface::update :"
            << " have-facesPtr_:" << facesPtr_.valid()
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


Foam::tmp<Foam::scalarField> Foam::distanceSurface::sample
(
    const volScalarField& vField
) const
{
    return sampleField(vField);
}


Foam::tmp<Foam::vectorField> Foam::distanceSurface::sample
(
    const volVectorField& vField
) const
{
    return sampleField(vField);
}


Foam::tmp<Foam::sphericalTensorField> Foam::distanceSurface::sample
(
    const volSphericalTensorField& vField
) const
{
    return sampleField(vField);
}


Foam::tmp<Foam::symmTensorField> Foam::distanceSurface::sample
(
    const volSymmTensorField& vField
) const
{
    return sampleField(vField);
}


Foam::tmp<Foam::tensorField> Foam::distanceSurface::sample
(
    const volTensorField& vField
) const
{
    return sampleField(vField);
}


Foam::tmp<Foam::scalarField> Foam::distanceSurface::interpolate
(
    const interpolation<scalar>& interpolator
) const
{
    return interpolateField(interpolator);
}


Foam::tmp<Foam::vectorField> Foam::distanceSurface::interpolate
(
    const interpolation<vector>& interpolator
) const
{
    return interpolateField(interpolator);
}

Foam::tmp<Foam::sphericalTensorField> Foam::distanceSurface::interpolate
(
    const interpolation<sphericalTensor>& interpolator
) const
{
    return interpolateField(interpolator);
}


Foam::tmp<Foam::symmTensorField> Foam::distanceSurface::interpolate
(
    const interpolation<symmTensor>& interpolator
) const
{
    return interpolateField(interpolator);
}


Foam::tmp<Foam::tensorField> Foam::distanceSurface::interpolate
(
    const interpolation<tensor>& interpolator
) const
{
    return interpolateField(interpolator);
}


void Foam::distanceSurface::print(Ostream& os) const
{
    os  << "distanceSurface: " << name() << " :"
        << "  surface:" << surfPtr_().name()
        << "  distance:" << distance_
        << "  faces:" << faces().size()
        << "  points:" << points().size();
}


// ************************************************************************* //
