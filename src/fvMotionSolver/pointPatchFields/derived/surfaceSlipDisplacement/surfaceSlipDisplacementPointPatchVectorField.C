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

#include "surfaceSlipDisplacementPointPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "Time.H"
#include "transformField.H"
#include "dynamicMotionSolverFvMesh.H"
#include "displacementMotionSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<>
const char*
NamedEnum<surfaceSlipDisplacementPointPatchVectorField::projectMode, 3>::
names[] =
{
    "nearest",
    "pointNormal",
    "fixedNormal"
};

const NamedEnum<surfaceSlipDisplacementPointPatchVectorField::projectMode, 3>
    surfaceSlipDisplacementPointPatchVectorField::projectModeNames_;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void surfaceSlipDisplacementPointPatchVectorField::calcProjection
(
    vectorField& displacement
) const
{
    const polyMesh& mesh = patch().boundaryMesh().mesh()();
    const pointField& localPoints = patch().localPoints();
    const labelList& meshPoints = patch().meshPoints();

    // const scalar deltaT = mesh.time().deltaTValue();

    // Construct large enough vector in direction of projectDir so
    // we're guaranteed to hit something.

    //- Per point projection vector:
    const scalar projectLen = mag(mesh.bounds().max()-mesh.bounds().min());

    // For case of fixed projection vector:
    vector projectVec(0, 0, 0);
    if (projectMode_ == FIXEDNORMAL)
    {
        vector n = projectDir_/mag(projectDir_);
        projectVec = projectLen*n;
    }


    // Get fixed points (bit of a hack)
    const pointZone* zonePtr = nullptr;

    if (frozenPointsZone_.size() > 0)
    {
        const meshPointZones& pZones = mesh.pointZones();

        zonePtr = &pZones[frozenPointsZone_];

        Pout<< "surfaceSlipDisplacementPointPatchVectorField : Fixing all "
            << zonePtr->size() << " points in pointZone " << zonePtr->name()
            << endl;
    }

    // Get the motionSolver from the dynamic mesh
    const motionSolver& motion =
        refCast<const dynamicMotionSolverFvMesh>(mesh).motion();

    // Get the starting locations from the motionSolver
    const pointField& points0 =
        refCast<const displacementMotionSolver>(motion).points0();

    pointField start(meshPoints.size());
    forAll(start, i)
    {
        start[i] = points0[meshPoints[i]] + displacement[i];
    }

    label nNotProjected = 0;

    if (projectMode_ == NEAREST)
    {
        List<pointIndexHit> nearest;
        labelList hitSurfaces;
        surfaces().findNearest
        (
            start,
            scalarField(start.size(), sqr(projectLen)),
            hitSurfaces,
            nearest
        );

        forAll(nearest, i)
        {
            if (zonePtr && (zonePtr->whichPoint(meshPoints[i]) >= 0))
            {
                // Fixed point. Reset to point0 location.
                displacement[i] = points0[meshPoints[i]] - localPoints[i];
            }
            else if (nearest[i].hit())
            {
                displacement[i] =
                    nearest[i].hitPoint()
                  - points0[meshPoints[i]];
            }
            else
            {
                nNotProjected++;

                if (debug)
                {
                    Pout<< "    point:" << meshPoints[i]
                        << " coord:" << localPoints[i]
                        << "  did not find any surface within " << projectLen
                        << endl;
                }
            }
        }
    }
    else
    {
        // Do tests on all points. Combine later on.

        // 1. Check if already on surface
        List<pointIndexHit> nearest;
        {
            labelList nearestSurface;
            surfaces().findNearest
            (
                start,
                scalarField(start.size(), sqr(small)),
                nearestSurface,
                nearest
            );
        }

        // 2. intersection. (combined later on with information from nearest
        // above)
        vectorField projectVecs(start.size(), projectVec);

        if (projectMode_ == POINTNORMAL)
        {
            projectVecs = projectLen*patch().pointNormals();
        }

        // Knock out any wedge component
        scalarField offset(start.size(), 0.0);
        if (wedgePlane_ >= 0 && wedgePlane_ < vector::nComponents)
        {
            forAll(offset, i)
            {
                offset[i] = start[i][wedgePlane_];
                start[i][wedgePlane_] = 0;
                projectVecs[i][wedgePlane_] = 0;
            }
        }

        List<pointIndexHit> rightHit;
        {
            labelList rightSurf;
            surfaces().findAnyIntersection
            (
                start,
                start+projectVecs,
                rightSurf,
                rightHit
            );
        }

        List<pointIndexHit> leftHit;
        {
            labelList leftSurf;
            surfaces().findAnyIntersection
            (
                start,
                start-projectVecs,
                leftSurf,
                leftHit
            );
        }

        // 3. Choose either -fixed, nearest, right, left.
        forAll(displacement, i)
        {
            if (zonePtr && (zonePtr->whichPoint(meshPoints[i]) >= 0))
            {
                // Fixed point. Reset to point0 location.
                displacement[i] = points0[meshPoints[i]] - localPoints[i];
            }
            else if (nearest[i].hit())
            {
                // Found nearest.
                displacement[i] =
                    nearest[i].hitPoint()
                  - points0[meshPoints[i]];
            }
            else
            {
                pointIndexHit interPt;

                if (rightHit[i].hit())
                {
                    if (leftHit[i].hit())
                    {
                        if
                        (
                            magSqr(rightHit[i].hitPoint()-start[i])
                          < magSqr(leftHit[i].hitPoint()-start[i])
                        )
                        {
                            interPt = rightHit[i];
                        }
                        else
                        {
                            interPt = leftHit[i];
                        }
                    }
                    else
                    {
                        interPt = rightHit[i];
                    }
                }
                else
                {
                    if (leftHit[i].hit())
                    {
                        interPt = leftHit[i];
                    }
                }


                if (interPt.hit())
                {
                    if (wedgePlane_ >= 0 && wedgePlane_ < vector::nComponents)
                    {
                        interPt.rawPoint()[wedgePlane_] += offset[i];
                    }
                    displacement[i] = interPt.rawPoint()-points0[meshPoints[i]];
                }
                else
                {
                    nNotProjected++;

                    if (debug)
                    {
                        Pout<< "    point:" << meshPoints[i]
                            << " coord:" << localPoints[i]
                            << "  did not find any intersection between"
                            << " ray from " << start[i]-projectVecs[i]
                            << " to " << start[i]+projectVecs[i] << endl;
                    }
                }
            }
        }
    }

    reduce(nNotProjected, sumOp<label>());

    if (nNotProjected > 0)
    {
        Info<< "surfaceSlipDisplacement :"
            << " on patch " << patch().name()
            << " did not project " << nNotProjected
            << " out of " << returnReduce(localPoints.size(), sumOp<label>())
            << " points." << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

surfaceSlipDisplacementPointPatchVectorField::
surfaceSlipDisplacementPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF
)
:
    pointPatchVectorField(p, iF),
    projectMode_(NEAREST),
    projectDir_(Zero),
    wedgePlane_(-1)
{}


surfaceSlipDisplacementPointPatchVectorField::
surfaceSlipDisplacementPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const dictionary& dict
)
:
    pointPatchVectorField(p, iF, dict),
    surfacesDict_(dict.subDict("geometry")),
    projectMode_(projectModeNames_.read(dict.lookup("projectMode"))),
    projectDir_(dict.lookup("projectDirection")),
    wedgePlane_(dict.lookupOrDefault("wedgePlane", -1)),
    frozenPointsZone_(dict.lookupOrDefault("frozenPointsZone", word::null))
{}


surfaceSlipDisplacementPointPatchVectorField::
surfaceSlipDisplacementPointPatchVectorField
(
    const surfaceSlipDisplacementPointPatchVectorField& ppf,
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const pointPatchFieldMapper&
)
:
    pointPatchVectorField(p, iF),
    surfacesDict_(ppf.surfacesDict_),
    projectMode_(ppf.projectMode_),
    projectDir_(ppf.projectDir_),
    wedgePlane_(ppf.wedgePlane_),
    frozenPointsZone_(ppf.frozenPointsZone_)
{}


surfaceSlipDisplacementPointPatchVectorField::
surfaceSlipDisplacementPointPatchVectorField
(
    const surfaceSlipDisplacementPointPatchVectorField& ppf,
    const DimensionedField<vector, pointMesh>& iF
)
:
    pointPatchVectorField(ppf, iF),
    surfacesDict_(ppf.surfacesDict_),
    projectMode_(ppf.projectMode_),
    projectDir_(ppf.projectDir_),
    wedgePlane_(ppf.wedgePlane_),
    frozenPointsZone_(ppf.frozenPointsZone_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const searchableSurfaces&
surfaceSlipDisplacementPointPatchVectorField::surfaces() const
{
    if (surfacesPtr_.empty())
    {
        surfacesPtr_.reset
        (
            new searchableSurfaces
            (
                IOobject
                (
                    "abc",                  // dummy name
                    db().time().constant(),
                    searchableSurface::geometryDir(db().time()),
                    db().time(),
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                surfacesDict_,
                true    // use single region naming shortcut
            )
        );
    }
    return surfacesPtr_();
}


void surfaceSlipDisplacementPointPatchVectorField::evaluate
(
    const Pstream::commsTypes commsType
)
{
    vectorField displacement(this->patchInternalField());

    // Calculate displacement to project points onto surface
    calcProjection(displacement);

    // Get internal field to insert values into
    Field<vector>& iF = const_cast<Field<vector>&>(this->primitiveField());

    // setInInternalField(iF, motionU);
    setInInternalField(iF, displacement);

    pointPatchVectorField::evaluate(commsType);
}


void surfaceSlipDisplacementPointPatchVectorField::write(Ostream& os) const
{
    pointPatchVectorField::write(os);
    writeEntry(os, "geometry", surfacesDict_);
    writeEntry(os, "projectMode", projectModeNames_[projectMode_]);
    writeEntry(os, "projectDirection", projectDir_);
    writeEntry(os, "wedgePlane", wedgePlane_);
    if (frozenPointsZone_ != word::null)
    {
        writeEntry(os, "frozenPointsZone", frozenPointsZone_);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePointPatchTypeField
(
    pointPatchVectorField,
    surfaceSlipDisplacementPointPatchVectorField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
