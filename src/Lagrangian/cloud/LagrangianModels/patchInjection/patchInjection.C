/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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

#include "LagrangianSubMesh.H"
#include "patchInjection.H"
#include "LagrangianFields.H"
#include "cloud.H"
#include "cloudFunctionObjectUList.H"
#include "tetIndices.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace Lagrangian
{
    defineTypeNameAndDebug(patchInjection, 0);
    addToRunTimeSelectionTable(LagrangianModel, patchInjection, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::Lagrangian::patchInjection::readCoeffs(const dictionary& modelDict)
{
    patchName_ = modelDict.lookup<word>("patch");

    numberRate_.reset
    (
        Function1<scalar>::New
        (
            "numberRate",
            mesh().time().userUnits(),
            dimRate,
            modelDict
        ).ptr()
    );

    numberDeferred_ = 0;
}


void Foam::Lagrangian::patchInjection::calcSumAreas()
{
    const polyPatch& patch = mesh().mesh().boundaryMesh()[patchName_];

    // Create a point field mid-way through the fraction. We distribute
    // particles uniformly through the time-step, so using the mid-time-step
    // geometry ensures that the particles get distributed evenly even in the
    // presence of motion.
    const tmp<pointField> tPoints
    (
        mesh().mesh().moving()
      ? mesh().mesh().oldPoints()/2 + mesh().mesh().points()/2
      : tmp<pointField>(mesh().mesh().points())
    );
    const pointField& points = tPoints();

    // Count the number of triangles in each face
    labelList patchFaceNTris(patch.size(), label(0));
    forAll(patch, patchFacei)
    {
        patchFaceNTris[patchFacei] = patch[patchFacei].nTriangles();
    }

    // Construct cumulative sums of the areas in the patch faces and patch
    // face-triangles
    patchFaceSumArea_.resize(patch.size());
    patchFaceSumArea_ = 0;
    patchFaceTriSumArea_.resize(patchFaceNTris);
    forAll(patch, patchFacei)
    {
        const label facei = patch.start() + patchFacei;
        const face& f = patch[patchFacei];

        for (label fTrii = 0; fTrii < f.nTriangles(); ++ fTrii)
        {
            patchFaceSumArea_[patchFacei] +=
                tetIndices(mesh().mesh().faceOwner()[facei], facei, fTrii + 1)
               .faceTri(mesh().mesh(), points)
               .mag();

            patchFaceTriSumArea_[patchFacei][fTrii] =
                patchFaceSumArea_[patchFacei];
        }
    }
    for (label patchFacei = 1; patchFacei < patch.size(); ++ patchFacei)
    {
        patchFaceSumArea_[patchFacei] += patchFaceSumArea_[patchFacei - 1];
    }

    // Construct a cumulative sum of the areas across the processes
    procSumArea_ = scalar(-vGreat);
    procSumArea_[Pstream::myProcNo()] = patchFaceSumArea_.last();
    Pstream::listCombineGather(procSumArea_, maxEqOp<scalar>());
    Pstream::listCombineScatter(procSumArea_);
    for (label proci = 1; proci < Pstream::nProcs(); proci ++)
    {
        procSumArea_[proci] += procSumArea_[proci - 1];
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Lagrangian::patchInjection::patchInjection
(
    const word& name,
    const LagrangianMesh& mesh,
    const dictionary& modelDict,
    const dictionary& stateDict
)
:
    LagrangianInjection(name, mesh),
    cloudLagrangianModel(static_cast<const LagrangianModel&>(*this)),
    patchName_(word::null),
    numberRate_(nullptr),
    numberDeferred_(stateDict.lookupOrDefault<scalar>("numberDeferred", 0)),
    globalRndGen_("globalRndGen", stateDict, name, true),
    localRndGen_("localRndGen", stateDict, name, false),
    timeIndex_(-1),
    procSumArea_(Pstream::nProcs()),
    patchFaceSumArea_(),
    patchFaceTriSumArea_()
{
    readCoeffs(modelDict);

    calcSumAreas();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::LagrangianSubMesh Foam::Lagrangian::patchInjection::modify
(
    LagrangianMesh& mesh,
    const LagrangianSubMesh&
) const
{
    const scalar t1 = mesh.time().value();
    const scalar t0 = t1 - mesh.time().deltaT().value();

    // Restart the generators if necessary and set the time index up to date
    localRndGen_.start(timeIndex_ == db().time().timeIndex());
    globalRndGen_.start(timeIndex_ == db().time().timeIndex());
    timeIndex_ = db().time().timeIndex();

    // Calculate the number of particles to inject. Round down to get an
    // integer number. Store the excess to apply at a later time.
    const scalar number = numberRate_->integral(t0, t1) + numberDeferred_;
    const label numberInt = floor(number);
    numberDeferred_ = number - numberInt;

    // Inject at random times throughout the time-step
    scalarField fraction(globalRndGen_.scalar01(numberInt));

    // Binary search for a given volume within a supplied cumulative volume sum.
    // Once the interval has been found, subtract the lower value from that
    // volume so that a subsequent binary search can be done on a sub-set of
    // volumes. Note that the supplied cumulative volume sum does not include
    // the leading zero.
    auto findVolume = []
    (
        const scalarUList& volumes,
        scalar& volume
    )
    {
        label i0 = -1, i1 = volumes.size() - 1;

        while (i0 + 1 != i1)
        {
            const label i = (i0 + i1)/2;

            (volume < volumes[i] ? i1 : i0) = i;
        }

        if (i0 != -1)
        {
            volume -= volumes[i0];
        }

        return i1;
    };

    // Get the mesh index of the first face in the patch
    const polyPatch& patch = mesh.mesh().boundaryMesh()[patchName_];
    const label facei0 = patch.start();

    // Initialise storage for the injection geometry and topology. This is
    // dynamic as we don't know how much will end up on each processor yet.
    DynamicList<barycentric> injectCoordinates(numberInt/Pstream::nProcs());
    DynamicList<label> injectCells(numberInt/Pstream::nProcs());
    DynamicList<label> injectFaces(numberInt/Pstream::nProcs());
    DynamicList<label> injectFaceTris(numberInt/Pstream::nProcs());

    // Create a (global) list of areas at which to inject. Each area is
    // searched for in the cumulative lists to identify a triangle into
    // which to inject.
    scalarField area(globalRndGen_.scalar01(numberInt)*procSumArea_.last());
    forAll(area, areai)
    {
        const label proci = findVolume(procSumArea_, area[areai]);

        if (Pstream::myProcNo() == proci)
        {
            const label patchFacei =
                findVolume(patchFaceSumArea_, area[areai]);
            const label facei = facei0 + patchFacei;
            const label faceTrii =
                findVolume(patchFaceTriSumArea_[patchFacei], area[areai]);

            const barycentric2D r = barycentric2D01(localRndGen_);

            injectCoordinates.append(barycentric(0, r.a(), r.b(), r.c()));
            injectCells.append(mesh.mesh().faceOwner()[facei]);
            injectFaces.append(facei);
            injectFaceTris.append(faceTrii + 1);
        }
    }

    // Create the particles at the identified locations
    LagrangianSubMesh injectionMesh =
        mesh.inject
        (
            *this,
            barycentricField(injectCoordinates),
            labelField(injectCells),
            labelField(injectFaces),
            labelField(injectFaceTris),
            LagrangianMesh::fractionName,
            fraction
        );

    // Execute functions for this sub mesh. This ensures that flux fields and
    // similar include the amount injected through the patch faces.
    cloudFunctionObjectUList(cloud(), true).postCrossFaces
    (
        LagrangianSubMesh
        (
            injectionMesh.mesh(),
            static_cast<LagrangianGroup>
            (
                static_cast<label>(LagrangianGroup::onPatchZero)
              + patch.index()
            ),
            injectionMesh.size(),
            injectionMesh.start()
        ).sub
        (
            mesh.lookupObject<LagrangianScalarInternalDynamicField>
            (
                LagrangianMesh::fractionName
            )
        )
    );

    return injectionMesh;
}


void Foam::Lagrangian::patchInjection::topoChange(const polyTopoChangeMap& map)
{
    calcSumAreas();
}


void Foam::Lagrangian::patchInjection::mapMesh(const polyMeshMap& map)
{
    calcSumAreas();
}


void Foam::Lagrangian::patchInjection::distribute
(
    const polyDistributionMap& map
)
{
    calcSumAreas();
}


void Foam::Lagrangian::patchInjection::correct()
{
    if (mesh().mesh().moving())
    {
        calcSumAreas();
    }
}


bool Foam::Lagrangian::patchInjection::read(const dictionary& modelDict)
{
    if (LagrangianInjection::read(modelDict))
    {
        readCoeffs(modelDict);
        return true;
    }
    else
    {
        return false;
    }
}


void Foam::Lagrangian::patchInjection::writeState(Ostream& os) const
{
    LagrangianInjection::writeState(os);

    writeEntry(os, "numberDeferred", numberDeferred_);
    writeEntry(os, "globalRndGen", globalRndGen_);
}


void Foam::Lagrangian::patchInjection::writeProcessorState(Ostream& os) const
{
    LagrangianInjection::writeProcessorState(os);

    writeEntry(os, "localRndGen", localRndGen_);
}


// ************************************************************************* //
