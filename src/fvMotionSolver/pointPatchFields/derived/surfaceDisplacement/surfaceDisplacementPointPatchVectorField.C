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

#include "surfaceDisplacementPointPatchVectorField.H"
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
NamedEnum<surfaceDisplacementPointPatchVectorField::projectMode, 3>::
names[] =
{
    "nearest",
    "pointNormal",
    "fixedNormal"
};

const NamedEnum<surfaceDisplacementPointPatchVectorField::projectMode, 3>
    surfaceDisplacementPointPatchVectorField::projectModeNames_;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void surfaceDisplacementPointPatchVectorField::calcProjection
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
    vector projectVec(Zero);
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

        Pout<< "surfaceDisplacementPointPatchVectorField : Fixing all "
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
        if (wedgePlane_ >= 0 && wedgePlane_ <= vector::nComponents)
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
                    if (wedgePlane_ >= 0 && wedgePlane_ <= vector::nComponents)
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
        Info<< "surfaceDisplacement :"
            << " on patch " << patch().name()
            << " did not project " << nNotProjected
            << " out of " << returnReduce(localPoints.size(), sumOp<label>())
            << " points." << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

surfaceDisplacementPointPatchVectorField::
surfaceDisplacementPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF
)
:
    fixedValuePointPatchVectorField(p, iF),
    velocity_(Zero),
    projectMode_(NEAREST),
    projectDir_(Zero),
    wedgePlane_(-1)
{}


surfaceDisplacementPointPatchVectorField::
surfaceDisplacementPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const dictionary& dict
)
:
    fixedValuePointPatchVectorField(p, iF, dict),
    velocity_(dict.lookup("velocity")),
    surfacesDict_(dict.subDict("geometry")),
    projectMode_(projectModeNames_.read(dict.lookup("projectMode"))),
    projectDir_(dict.lookup("projectDirection")),
    wedgePlane_(dict.lookupOrDefault("wedgePlane", -1)),
    frozenPointsZone_(dict.lookupOrDefault("frozenPointsZone", word::null))
{
    if (velocity_.x() < 0 || velocity_.y() < 0 || velocity_.z() < 0)
    {
        FatalErrorInFunction
            << "All components of velocity have to be positive : "
            << velocity_ << nl
            << "Set velocity components to a great value if no clipping"
            << " necessary." << exit(FatalError);
    }
}


surfaceDisplacementPointPatchVectorField::
surfaceDisplacementPointPatchVectorField
(
    const surfaceDisplacementPointPatchVectorField& ppf,
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const pointPatchFieldMapper& mapper
)
:
    fixedValuePointPatchVectorField(ppf, p, iF, mapper),
    velocity_(ppf.velocity_),
    surfacesDict_(ppf.surfacesDict_),
    projectMode_(ppf.projectMode_),
    projectDir_(ppf.projectDir_),
    wedgePlane_(ppf.wedgePlane_),
    frozenPointsZone_(ppf.frozenPointsZone_)
{}


surfaceDisplacementPointPatchVectorField::
surfaceDisplacementPointPatchVectorField
(
    const surfaceDisplacementPointPatchVectorField& ppf,
    const DimensionedField<vector, pointMesh>& iF
)
:
    fixedValuePointPatchVectorField(ppf, iF),
    velocity_(ppf.velocity_),
    surfacesDict_(ppf.surfacesDict_),
    projectMode_(ppf.projectMode_),
    projectDir_(ppf.projectDir_),
    wedgePlane_(ppf.wedgePlane_),
    frozenPointsZone_(ppf.frozenPointsZone_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const searchableSurfaces&
surfaceDisplacementPointPatchVectorField::surfaces() const
{
    if (surfacesPtr_.empty())
    {
        surfacesPtr_.reset
        (
            new searchableSurfaces
            (
                IOobject
                (
                    "abc",             // dummy name
                    db().time().constant(),
                    searchableSurface::geometryDir(db().time()),
                    db().time(),
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                surfacesDict_,
                true                // allow single-region shortcut
            )
        );
    }
    return surfacesPtr_();
}


void surfaceDisplacementPointPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    const polyMesh& mesh = patch().boundaryMesh().mesh()();

    vectorField currentDisplacement(this->patchInternalField());

    // Calculate intersections with surface w.r.t points0.
    vectorField displacement(currentDisplacement);
    calcProjection(displacement);

    // offset wrt current displacement
    vectorField offset(displacement-currentDisplacement);

    // Clip offset to maximum displacement possible: velocity*timestep

    const scalar deltaT = mesh.time().deltaTValue();
    const vector clipVelocity = velocity_*deltaT;

    forAll(displacement, i)
    {
        vector& d = offset[i];

        for (direction cmpt = 0; cmpt < vector::nComponents; cmpt++)
        {
            if (d[cmpt] < 0)
            {
                d[cmpt] = max(d[cmpt], -clipVelocity[cmpt]);
            }
            else
            {
                d[cmpt] = min(d[cmpt], clipVelocity[cmpt]);
            }
        }
    }

    this->operator==(currentDisplacement+offset);

    fixedValuePointPatchVectorField::updateCoeffs();
}


void surfaceDisplacementPointPatchVectorField::write(Ostream& os) const
{
    fixedValuePointPatchVectorField::write(os);
    writeEntry(os, "velocity", velocity_);
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
    fixedValuePointPatchVectorField,
    surfaceDisplacementPointPatchVectorField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
