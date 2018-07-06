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

#include "sampledCuttingPlane.H"
#include "dictionary.H"
#include "volFields.H"
#include "volPointInterpolation.H"
#include "addToRunTimeSelectionTable.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(sampledCuttingPlane, 0);
    addNamedToRunTimeSelectionTable
    (
        sampledSurface,
        sampledCuttingPlane,
        word,
        cuttingPlane
    );
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::sampledCuttingPlane::createGeometry()
{
    if (debug)
    {
        Pout<< "sampledCuttingPlane::createGeometry :updating geometry."
            << endl;
    }

    // Clear any stored topologies
    facesPtr_.clear();
    isoSurfPtr_.ptr();
    pointDistance_.clear();
    cellDistancePtr_.clear();

    // Clear derived data
    clearGeom();

    // Get any subMesh
    if (zoneID_.index() != -1 && !subMeshPtr_.valid())
    {
        const polyBoundaryMesh& patches = mesh().boundaryMesh();

        // Patch to put exposed internal faces into
        const label exposedPatchi = patches.findPatchID(exposedPatchName_);

        DebugInfo
            << "Allocating subset of size "
            << mesh().cellZones()[zoneID_.index()].size()
            << " with exposed faces into patch "
            << patches[exposedPatchi].name() << endl;

        subMeshPtr_.reset
        (
            new fvMeshSubset(static_cast<const fvMesh&>(mesh()))
        );
        subMeshPtr_().setLargeCellSubset
        (
            labelHashSet(mesh().cellZones()[zoneID_.index()]),
            exposedPatchi
        );
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
            dimensionedScalar("zero", dimLength, 0)
        )
    );
    volScalarField& cellDistance = cellDistancePtr_();

    // Internal field
    {
        const pointField& cc = mesh.cellCentres();
        scalarField& fld = cellDistance.primitiveFieldRef();

        forAll(cc, i)
        {
            // Signed distance
            fld[i] = (cc[i] - plane_.refPoint()) & plane_.normal();
        }
    }

    volScalarField::Boundary& cellDistanceBf =
        cellDistance.boundaryFieldRef();

    // Patch fields
    {
        forAll(cellDistanceBf, patchi)
        {
            if
            (
                isA<emptyFvPatchScalarField>
                (
                    cellDistanceBf[patchi]
                )
            )
            {
                cellDistanceBf.set
                (
                    patchi,
                    new calculatedFvPatchScalarField
                    (
                        mesh.boundary()[patchi],
                        cellDistance
                    )
                );

                const polyPatch& pp = mesh.boundary()[patchi].patch();
                pointField::subField cc = pp.patchSlice(mesh.faceCentres());

                fvPatchScalarField& fld = cellDistanceBf[patchi];
                fld.setSize(pp.size());
                forAll(fld, i)
                {
                    fld[i] = (cc[i] - plane_.refPoint()) & plane_.normal();
                }
            }
            else
            {
                const pointField& cc = mesh.C().boundaryField()[patchi];
                fvPatchScalarField& fld = cellDistanceBf[patchi];

                forAll(fld, i)
                {
                    fld[i] = (cc[i] - plane_.refPoint()) & plane_.normal();
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

        forAll(pointDistance_, i)
        {
            pointDistance_[i] = (pts[i] - plane_.refPoint()) & plane_.normal();
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
            dimensionedScalar("zero", dimLength, 0)
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
            cellDistance,
            pointDistance_,
            0.0,
            regularise_,
            mergeTol_
        )
        // new isoSurfaceCell
        //(
        //    mesh,
        //    cellDistance,
        //    pointDistance_,
        //    0.0,
        //    regularise_,
        //    mergeTol_
        //)
    );

    if (debug)
    {
        print(Pout);
        Pout<< endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sampledCuttingPlane::sampledCuttingPlane
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    sampledSurface(name, mesh, dict),
    plane_(dict),
    mergeTol_(dict.lookupOrDefault("mergeTol", 1e-6)),
    regularise_(dict.lookupOrDefault("regularise", true)),
    average_(dict.lookupOrDefault("average", false)),
    zoneID_(dict.lookupOrDefault("zone", word::null), mesh.cellZones()),
    exposedPatchName_(word::null),
    needsUpdate_(true),
    subMeshPtr_(nullptr),
    cellDistancePtr_(nullptr),
    isoSurfPtr_(nullptr),
    facesPtr_(nullptr)
{
    if (zoneID_.index() != -1)
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

        if (debug && zoneID_.index() != -1)
        {
            Info<< "Restricting to cellZone " << zoneID_.name()
                << " with exposed internal faces into patch "
                << exposedPatchName_ << endl;
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sampledCuttingPlane::~sampledCuttingPlane()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::sampledCuttingPlane::needsUpdate() const
{
    return needsUpdate_;
}


bool Foam::sampledCuttingPlane::expire()
{
    if (debug)
    {
        Pout<< "sampledCuttingPlane::expire :"
            << " have-facesPtr_:" << facesPtr_.valid()
            << " needsUpdate_:" << needsUpdate_ << endl;
    }

    // Clear any stored topologies
    facesPtr_.clear();

    // Clear derived data
    clearGeom();

    // Already marked as expired
    if (needsUpdate_)
    {
        return false;
    }

    needsUpdate_ = true;
    return true;
}


bool Foam::sampledCuttingPlane::update()
{
    if (debug)
    {
        Pout<< "sampledCuttingPlane::update :"
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


Foam::tmp<Foam::scalarField>
Foam::sampledCuttingPlane::sample
(
    const volScalarField& vField
) const
{
    return sampleField(vField);
}


Foam::tmp<Foam::vectorField>
Foam::sampledCuttingPlane::sample
(
    const volVectorField& vField
) const
{
    return sampleField(vField);
}


Foam::tmp<Foam::sphericalTensorField>
Foam::sampledCuttingPlane::sample
(
    const volSphericalTensorField& vField
) const
{
    return sampleField(vField);
}


Foam::tmp<Foam::symmTensorField>
Foam::sampledCuttingPlane::sample
(
    const volSymmTensorField& vField
) const
{
    return sampleField(vField);
}


Foam::tmp<Foam::tensorField>
Foam::sampledCuttingPlane::sample
(
    const volTensorField& vField
) const
{
    return sampleField(vField);
}


Foam::tmp<Foam::scalarField>
Foam::sampledCuttingPlane::interpolate
(
    const interpolation<scalar>& interpolator
) const
{
    return interpolateField(interpolator);
}


Foam::tmp<Foam::vectorField>
Foam::sampledCuttingPlane::interpolate
(
    const interpolation<vector>& interpolator
) const
{
    return interpolateField(interpolator);
}

Foam::tmp<Foam::sphericalTensorField>
Foam::sampledCuttingPlane::interpolate
(
    const interpolation<sphericalTensor>& interpolator
) const
{
    return interpolateField(interpolator);
}


Foam::tmp<Foam::symmTensorField>
Foam::sampledCuttingPlane::interpolate
(
    const interpolation<symmTensor>& interpolator
) const
{
    return interpolateField(interpolator);
}


Foam::tmp<Foam::tensorField>
Foam::sampledCuttingPlane::interpolate
(
    const interpolation<tensor>& interpolator
) const
{
    return interpolateField(interpolator);
}


void Foam::sampledCuttingPlane::print(Ostream& os) const
{
    os  << "sampledCuttingPlane: " << name() << " :"
        << "  plane:" << plane_
        << "  faces:" << faces().size()
        << "  points:" << points().size();
}


// ************************************************************************* //
