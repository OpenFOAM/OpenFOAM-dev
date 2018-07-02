/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
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

#include "displacementLayeredMotionMotionSolver.H"
#include "addToRunTimeSelectionTable.H"
#include "pointEdgeStructuredWalk.H"
#include "pointFields.H"
#include "PointEdgeWave.H"
#include "syncTools.H"
#include "interpolationTable.H"
#include "pointConstraints.H"
#include "mapPolyMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(displacementLayeredMotionMotionSolver, 0);

    addToRunTimeSelectionTable
    (
        motionSolver,
        displacementLayeredMotionMotionSolver,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::displacementLayeredMotionMotionSolver::calcZoneMask
(
    const label cellZoneI,
    PackedBoolList& isZonePoint,
    PackedBoolList& isZoneEdge
) const
{
    if (cellZoneI == -1)
    {
        isZonePoint.setSize(mesh().nPoints());
        isZonePoint = 1;

        isZoneEdge.setSize(mesh().nEdges());
        isZoneEdge = 1;
    }
    else
    {
        const cellZone& cz = mesh().cellZones()[cellZoneI];

        label nPoints = 0;
        forAll(cz, i)
        {
            const labelList& cPoints = mesh().cellPoints(cz[i]);
            forAll(cPoints, cPointi)
            {
                if (!isZonePoint[cPoints[cPointi]])
                {
                    isZonePoint[cPoints[cPointi]] = 1;
                    nPoints++;
                }
            }
        }
        syncTools::syncPointList
        (
            mesh(),
            isZonePoint,
            orEqOp<unsigned int>(),
            0
        );


        // Mark edge inside cellZone
        label nEdges = 0;
        forAll(cz, i)
        {
            const labelList& cEdges = mesh().cellEdges(cz[i]);
            forAll(cEdges, cEdgeI)
            {
                if (!isZoneEdge[cEdges[cEdgeI]])
                {
                    isZoneEdge[cEdges[cEdgeI]] = 1;
                    nEdges++;
                }
            }
        }
        syncTools::syncEdgeList
        (
            mesh(),
            isZoneEdge,
            orEqOp<unsigned int>(),
            0
        );

        if (debug)
        {
            Info<< "On cellZone " << cz.name()
                << " marked " << returnReduce(nPoints, sumOp<label>())
                << " points and " << returnReduce(nEdges, sumOp<label>())
                << " edges." << endl;
        }
    }
}


// Find distance to starting point
void Foam::displacementLayeredMotionMotionSolver::walkStructured
(
    const label cellZoneI,
    const PackedBoolList& isZonePoint,
    const PackedBoolList& isZoneEdge,
    const labelList& seedPoints,
    const vectorField& seedData,
    scalarField& distance,
    vectorField& data
) const
{
    List<pointEdgeStructuredWalk> seedInfo(seedPoints.size());

    forAll(seedPoints, i)
    {
        seedInfo[i] = pointEdgeStructuredWalk
        (
            points0()[seedPoints[i]],  // location of data
            points0()[seedPoints[i]],  // previous location
            0.0,
            seedData[i]
        );
    }

    // Current info on points
    List<pointEdgeStructuredWalk> allPointInfo(mesh().nPoints());

    // Mark points inside cellZone.
    // Note that we use points0, not mesh.points()
    // so as not to accumulate errors.
    forAll(isZonePoint, pointi)
    {
        if (isZonePoint[pointi])
        {
            allPointInfo[pointi] = pointEdgeStructuredWalk
            (
                points0()[pointi],  // location of data
                vector::max,        // not valid
                0.0,
                Zero        // passive data
            );
        }
    }

    // Current info on edges
    List<pointEdgeStructuredWalk> allEdgeInfo(mesh().nEdges());

    // Mark edges inside cellZone
    forAll(isZoneEdge, edgeI)
    {
        if (isZoneEdge[edgeI])
        {
            allEdgeInfo[edgeI] = pointEdgeStructuredWalk
            (
                mesh().edges()[edgeI].centre(points0()),    // location of data
                vector::max,                                // not valid
                0.0,
                Zero
            );
        }
    }

    // Walk
    PointEdgeWave<pointEdgeStructuredWalk> wallCalc
    (
        mesh(),
        seedPoints,
        seedInfo,

        allPointInfo,
        allEdgeInfo,
        mesh().globalData().nTotalPoints()  // max iterations
    );

    // Extract distance and passive data
    forAll(allPointInfo, pointi)
    {
        if (isZonePoint[pointi])
        {
            distance[pointi] = allPointInfo[pointi].dist();
            data[pointi] = allPointInfo[pointi].data();
        }
    }
}


// Evaluate faceZone patch
Foam::tmp<Foam::vectorField>
Foam::displacementLayeredMotionMotionSolver::faceZoneEvaluate
(
    const faceZone& fz,
    const labelList& meshPoints,
    const dictionary& dict,
    const PtrList<pointVectorField>& patchDisp,
    const label patchi
) const
{
    tmp<vectorField> tfld(new vectorField(meshPoints.size()));
    vectorField& fld = tfld.ref();

    const word type(dict.lookup("type"));

    if (type == "fixedValue")
    {
        fld = vectorField("value", dict, meshPoints.size());
    }
    else if (type == "timeVaryingUniformFixedValue")
    {
        interpolationTable<vector> timeSeries(dict);

        fld = timeSeries(mesh().time().timeOutputValue());
    }
    else if (type == "slip")
    {
        if ((patchi % 2) != 1)
        {
            FatalIOErrorInFunction(*this)
                << "FaceZone:" << fz.name()
                << exit(FatalIOError);
        }
        // Use field set by previous bc
        fld = vectorField(patchDisp[patchi - 1], meshPoints);
    }
    else if (type == "follow")
    {
        // Only on boundary faces - follow boundary conditions
        fld = vectorField(pointDisplacement_, meshPoints);
    }
    else if (type == "uniformFollow")
    {
        // Reads name of name of patch. Then get average point dislacement on
        // patch. That becomes the value of fld.
        const word patchName(dict.lookup("patch"));
        label patchID = mesh().boundaryMesh().findPatchID(patchName);
        pointField pdf
        (
            pointDisplacement_.boundaryField()[patchID].patchInternalField()
        );
        fld = gAverage(pdf);
    }
    else
    {
        FatalIOErrorInFunction(*this)
            << "Unknown faceZonePatch type " << type << " for faceZone "
            << fz.name() << exit(FatalIOError);
    }
    return tfld;
}


void Foam::displacementLayeredMotionMotionSolver::cellZoneSolve
(
    const label cellZoneI,
    const dictionary& zoneDict
)
{
    PackedBoolList isZonePoint(mesh().nPoints());
    PackedBoolList isZoneEdge(mesh().nEdges());
    calcZoneMask(cellZoneI, isZonePoint, isZoneEdge);

    const dictionary& patchesDict = zoneDict.subDict("boundaryField");

    if (patchesDict.size() != 2)
    {
        FatalIOErrorInFunction(*this)
            << "Two faceZones (patches) must be specified per cellZone. "
            << " cellZone:" << cellZoneI
            << " patches:" << patchesDict.toc()
            << exit(FatalIOError);
    }

    PtrList<scalarField> patchDist(patchesDict.size());
    PtrList<pointVectorField> patchDisp(patchesDict.size());

    // Allocate the fields
    label patchi = 0;
    forAllConstIter(dictionary, patchesDict, patchiter)
    {
        const word& faceZoneName = patchiter().keyword();
        label zoneI = mesh().faceZones().findZoneID(faceZoneName);
        if (zoneI == -1)
        {
            FatalIOErrorInFunction(*this)
                << "Cannot find faceZone " << faceZoneName
                << endl << "Valid zones are " << mesh().faceZones().names()
                << exit(FatalIOError);
        }

        // Determine the points of the faceZone within the cellZone
        const faceZone& fz = mesh().faceZones()[zoneI];

        patchDist.set(patchi, new scalarField(mesh().nPoints()));
        patchDisp.set
        (
            patchi,
            new pointVectorField
            (
                IOobject
                (
                    mesh().cellZones()[cellZoneI].name() + "_" + fz.name(),
                    mesh().time().timeName(),
                    mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                pointDisplacement_  // to inherit the boundary conditions
            )
        );

        patchi++;
    }



    // 'correctBoundaryConditions'
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Loops over all the faceZones and walks their boundary values

    // Make sure we can pick up bc values from field
    pointDisplacement_.correctBoundaryConditions();

    patchi = 0;
    forAllConstIter(dictionary, patchesDict, patchiter)
    {
        const word& faceZoneName = patchiter().keyword();
        const dictionary& faceZoneDict = patchiter().dict();

        // Determine the points of the faceZone within the cellZone
        const faceZone& fz = mesh().faceZones()[faceZoneName];
        const labelList& fzMeshPoints = fz().meshPoints();
        DynamicList<label> meshPoints(fzMeshPoints.size());
        forAll(fzMeshPoints, i)
        {
            if (isZonePoint[fzMeshPoints[i]])
            {
                meshPoints.append(fzMeshPoints[i]);
            }
        }

        // Get initial value for all the faceZone points
        tmp<vectorField> tseed = faceZoneEvaluate
        (
            fz,
            meshPoints,
            faceZoneDict,
            patchDisp,
            patchi
        );

        if (debug)
        {
            Info<< "For cellZone:" << cellZoneI
                << " for faceZone:" << fz.name()
                << " nPoints:" << tseed().size()
                << " have patchField:"
                << " max:" << gMax(tseed())
                << " min:" << gMin(tseed())
                << " avg:" << gAverage(tseed())
                << endl;
        }

        // Set distance and transported value
        walkStructured
        (
            cellZoneI,
            isZonePoint,
            isZoneEdge,

            meshPoints,
            tseed,
            patchDist[patchi],
            patchDisp[patchi]
        );

        // Implement real bc.
        patchDisp[patchi].correctBoundaryConditions();

        patchi++;
    }


    // Solve
    // ~~~~~

    if (debug)
    {
        // Normalised distance
        pointScalarField distance
        (
            IOobject
            (
                mesh().cellZones()[cellZoneI].name() + ":distance",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            pointMesh::New(mesh()),
            dimensionedScalar("zero", dimLength, 0.0)
        );

        forAll(distance, pointi)
        {
            scalar d1 = patchDist[0][pointi];
            scalar d2 = patchDist[1][pointi];
            if (d1 + d2 > small)
            {
                scalar s = d1/(d1 + d2);
                distance[pointi] = s;
            }
        }

        Info<< "Writing " << pointScalarField::typeName
            << distance.name() << " to "
            << mesh().time().timeName() << endl;
        distance.write();
    }


    const word interpolationScheme = zoneDict.lookup("interpolationScheme");

    if (interpolationScheme == "oneSided")
    {
        forAll(pointDisplacement_, pointi)
        {
            if (isZonePoint[pointi])
            {
                pointDisplacement_[pointi] = patchDisp[0][pointi];
            }
        }
    }
    else if (interpolationScheme == "linear")
    {
        forAll(pointDisplacement_, pointi)
        {
            if (isZonePoint[pointi])
            {
                scalar d1 = patchDist[0][pointi];
                scalar d2 = patchDist[1][pointi];
                scalar s = d1/(d1 + d2 + vSmall);

                const vector& pd1 = patchDisp[0][pointi];
                const vector& pd2 = patchDisp[1][pointi];

                pointDisplacement_[pointi] = (1 - s)*pd1 + s*pd2;
            }
        }
    }
    else
    {
        FatalErrorInFunction
            << "Invalid interpolationScheme: " << interpolationScheme
            << ". Valid schemes are 'oneSided' and 'linear'"
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::displacementLayeredMotionMotionSolver::
displacementLayeredMotionMotionSolver
(
    const polyMesh& mesh,
    const IOdictionary& dict
)
:
    displacementMotionSolver(mesh, dict, typeName)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::displacementLayeredMotionMotionSolver::
~displacementLayeredMotionMotionSolver()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::pointField>
Foam::displacementLayeredMotionMotionSolver::curPoints() const
{
    tmp<pointField> tcurPoints
    (
        points0() + pointDisplacement_.primitiveField()
    );

    return tcurPoints;
}


void Foam::displacementLayeredMotionMotionSolver::solve()
{
    // The points have moved so before interpolation update the motionSolver
    movePoints(mesh().points());

    // Apply boundary conditions
    pointDisplacement_.boundaryFieldRef().updateCoeffs();

    // Solve motion on all regions (=cellZones)
    const dictionary& regionDicts = coeffDict().subDict("regions");
    forAllConstIter(dictionary, regionDicts, regionIter)
    {
        const word& cellZoneName = regionIter().keyword();
        const dictionary& regionDict = regionIter().dict();

        label zoneI = mesh().cellZones().findZoneID(cellZoneName);

        Info<< "solving for zone: " << cellZoneName << endl;

        if (zoneI == -1)
        {
            FatalIOErrorInFunction(*this)
                << "Cannot find cellZone " << cellZoneName
                << endl << "Valid zones are " << mesh().cellZones().names()
                << exit(FatalIOError);
        }

        cellZoneSolve(zoneI, regionDict);
    }

    // Update pointDisplacement for solved values
    const pointConstraints& pcs =
        pointConstraints::New(pointDisplacement_.mesh());
    pcs.constrainDisplacement(pointDisplacement_, false);
}


void Foam::displacementLayeredMotionMotionSolver::updateMesh
(
    const mapPolyMesh& mpm
)
{
    FatalErrorInFunction
        << "Probably inconsistent with points0MotionSolver" << nl
        << "    Needs to be updated and tested."
        << exit(FatalError);

    displacementMotionSolver::updateMesh(mpm);

    const vectorField displacement(this->newPoints() - points0_);

    forAll(points0_, pointi)
    {
        const label oldPointi = mpm.pointMap()[pointi];

        if (oldPointi >= 0)
        {
            label masterPointi = mpm.reversePointMap()[oldPointi];

            if ((masterPointi != pointi))
            {
                // newly inserted point in this cellZone

                // need to set point0 so that it represents the position that
                // it would have had if it had existed for all time
                points0_[pointi] -= displacement[pointi];
            }
        }
    }
}


// ************************************************************************* //
