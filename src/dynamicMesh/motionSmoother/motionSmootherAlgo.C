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

#include "motionSmootherAlgo.H"
#include "twoDPointCorrector.H"
#include "faceSet.H"
#include "pointSet.H"
#include "fixedValuePointPatchFields.H"
#include "pointConstraints.H"
#include "syncTools.H"
#include "meshTools.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(motionSmootherAlgo, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::motionSmootherAlgo::testSyncPositions
(
    const pointField& fld,
    const scalar maxMag
) const
{
    pointField syncedFld(fld);

    syncTools::syncPointPositions
    (
        mesh_,
        syncedFld,
        minEqOp<point>(),           // combine op
        point(GREAT,GREAT,GREAT)    // null
    );

    forAll(syncedFld, i)
    {
        if (mag(syncedFld[i] - fld[i]) > maxMag)
        {
            FatalErrorIn
            (
                "motionSmootherAlgo::testSyncPositions"
                "("
                    "const pointField&, "
                    "const scalar"
                ")"
            )   << "On point " << i << " point:" << fld[i]
                << " synchronised point:" << syncedFld[i]
                << abort(FatalError);
        }
    }
}


void Foam::motionSmootherAlgo::checkFld(const pointScalarField& fld)
{
    forAll(fld, pointI)
    {
        const scalar val = fld[pointI];

        if ((val > -GREAT) && (val < GREAT))
        {}
        else
        {
            FatalErrorIn
            (
                "motionSmootherAlgo::checkFld"
                "(const pointScalarField&)"
            )   << "Problem : point:" << pointI << " value:" << val
                << abort(FatalError);
        }
    }
}


Foam::labelHashSet Foam::motionSmootherAlgo::getPoints
(
    const labelHashSet& faceLabels
) const
{
    labelHashSet usedPoints(mesh_.nPoints()/100);

    forAllConstIter(labelHashSet, faceLabels, iter)
    {
        const face& f = mesh_.faces()[iter.key()];

        forAll(f, fp)
        {
            usedPoints.insert(f[fp]);
        }
    }

    return usedPoints;
}


Foam::tmp<Foam::scalarField> Foam::motionSmootherAlgo::calcEdgeWeights
(
    const pointField& points
) const
{
    const edgeList& edges = mesh_.edges();

    tmp<scalarField> twght(new scalarField(edges.size()));
    scalarField& wght = twght();

    forAll(edges, edgeI)
    {
        wght[edgeI] = 1.0/(edges[edgeI].mag(points)+SMALL);
    }
    return twght;
}


// Smooth on selected points (usually patch points)
void Foam::motionSmootherAlgo::minSmooth
(
    const scalarField& edgeWeights,
    const PackedBoolList& isAffectedPoint,
    const labelList& meshPoints,
    const pointScalarField& fld,
    pointScalarField& newFld
) const
{
    tmp<pointScalarField> tavgFld = avg
    (
        fld,
        edgeWeights //scalarField(mesh_.nEdges(), 1.0)    // uniform weighting
    );
    const pointScalarField& avgFld = tavgFld();

    forAll(meshPoints, i)
    {
        label pointI = meshPoints[i];
        if (isAffectedPoint.get(pointI) == 1)
        {
            newFld[pointI] = min
            (
                fld[pointI],
                0.5*fld[pointI] + 0.5*avgFld[pointI]
            );
        }
    }

    // Single and multi-patch constraints
    pointConstraints::New(pMesh()).constrain(newFld, false);
}


// Smooth on all internal points
void Foam::motionSmootherAlgo::minSmooth
(
    const scalarField& edgeWeights,
    const PackedBoolList& isAffectedPoint,
    const pointScalarField& fld,
    pointScalarField& newFld
) const
{
    tmp<pointScalarField> tavgFld = avg
    (
        fld,
        edgeWeights //scalarField(mesh_.nEdges(), 1.0)    // uniform weighting
    );
    const pointScalarField& avgFld = tavgFld();

    forAll(fld, pointI)
    {
        if (isAffectedPoint.get(pointI) == 1 && isInternalPoint(pointI))
        {
            newFld[pointI] = min
            (
                fld[pointI],
                0.5*fld[pointI] + 0.5*avgFld[pointI]
            );
        }
    }

   // Single and multi-patch constraints
    pointConstraints::New(pMesh()).constrain(newFld, false);

}


// Scale on all internal points
void Foam::motionSmootherAlgo::scaleField
(
    const labelHashSet& pointLabels,
    const scalar scale,
    pointScalarField& fld
) const
{
    forAllConstIter(labelHashSet, pointLabels, iter)
    {
        if (isInternalPoint(iter.key()))
        {
            fld[iter.key()] *= scale;
        }
    }

    // Single and multi-patch constraints
    pointConstraints::New(pMesh()).constrain(fld, false);
}


// Scale on selected points (usually patch points)
void Foam::motionSmootherAlgo::scaleField
(
    const labelList& meshPoints,
    const labelHashSet& pointLabels,
    const scalar scale,
    pointScalarField& fld
) const
{
    forAll(meshPoints, i)
    {
        label pointI = meshPoints[i];

        if (pointLabels.found(pointI))
        {
            fld[pointI] *= scale;
        }
    }
}


// Lower on internal points
void Foam::motionSmootherAlgo::subtractField
(
    const labelHashSet& pointLabels,
    const scalar f,
    pointScalarField& fld
) const
{
    forAllConstIter(labelHashSet, pointLabels, iter)
    {
        if (isInternalPoint(iter.key()))
        {
            fld[iter.key()] = max(0.0, fld[iter.key()]-f);
        }
    }

    // Single and multi-patch constraints
    pointConstraints::New(pMesh()).constrain(fld);
}


// Scale on selected points (usually patch points)
void Foam::motionSmootherAlgo::subtractField
(
    const labelList& meshPoints,
    const labelHashSet& pointLabels,
    const scalar f,
    pointScalarField& fld
) const
{
    forAll(meshPoints, i)
    {
        label pointI = meshPoints[i];

        if (pointLabels.found(pointI))
        {
            fld[pointI] = max(0.0, fld[pointI]-f);
        }
    }
}


bool Foam::motionSmootherAlgo::isInternalPoint(const label pointI) const
{
    return isInternalPoint_.get(pointI) == 1;
}


void Foam::motionSmootherAlgo::getAffectedFacesAndPoints
(
    const label nPointIter,
    const faceSet& wrongFaces,

    labelList& affectedFaces,
    PackedBoolList& isAffectedPoint
) const
{
    isAffectedPoint.setSize(mesh_.nPoints());
    isAffectedPoint = 0;

    faceSet nbrFaces(mesh_, "checkFaces", wrongFaces);

    // Find possible points influenced by nPointIter iterations of
    // scaling and smoothing by doing pointCellpoint walk.
    // Also update faces-to-be-checked to extend one layer beyond the points
    // that will get updated.

    for (label i = 0; i < nPointIter; i++)
    {
        pointSet nbrPoints(mesh_, "grownPoints", getPoints(nbrFaces.toc()));

        forAllConstIter(pointSet, nbrPoints, iter)
        {
            const labelList& pCells = mesh_.pointCells(iter.key());

            forAll(pCells, pCellI)
            {
                const cell& cFaces = mesh_.cells()[pCells[pCellI]];

                forAll(cFaces, cFaceI)
                {
                    nbrFaces.insert(cFaces[cFaceI]);
                }
            }
        }
        nbrFaces.sync(mesh_);

        if (i == nPointIter - 2)
        {
            forAllConstIter(faceSet, nbrFaces, iter)
            {
                const face& f = mesh_.faces()[iter.key()];
                forAll(f, fp)
                {
                    isAffectedPoint.set(f[fp], 1);
                }
            }
        }
    }

    affectedFaces = nbrFaces.toc();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::motionSmootherAlgo::motionSmootherAlgo
(
    polyMesh& mesh,
    pointMesh& pMesh,
    indirectPrimitivePatch& pp,
    pointVectorField& displacement,
    pointScalarField& scale,
    pointField& oldPoints,
    const labelList& adaptPatchIDs,
    const dictionary& paramDict
)
:
    mesh_(mesh),
    pMesh_(pMesh),
    pp_(pp),
    displacement_(displacement),
    scale_(scale),
    oldPoints_(oldPoints),
    adaptPatchIDs_(adaptPatchIDs),
    paramDict_(paramDict),
    isInternalPoint_(mesh_.nPoints(), 1)
{
    updateMesh();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::motionSmootherAlgo::~motionSmootherAlgo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::polyMesh& Foam::motionSmootherAlgo::mesh() const
{
    return mesh_;
}


const Foam::pointMesh& Foam::motionSmootherAlgo::pMesh() const
{
    return pMesh_;
}


const Foam::indirectPrimitivePatch& Foam::motionSmootherAlgo::patch() const
{
    return pp_;
}


const Foam::labelList& Foam::motionSmootherAlgo::adaptPatchIDs() const
{
    return adaptPatchIDs_;
}


const Foam::dictionary& Foam::motionSmootherAlgo::paramDict() const
{
    return paramDict_;
}


void Foam::motionSmootherAlgo::correct()
{
    oldPoints_ = mesh_.points();

    scale_ = 1.0;

    // No need to update twoDmotion corrector since only holds edge labels
    // which will remain the same as before. So unless the mesh was distorted
    // severely outside of motionSmootherAlgo there will be no need.
}


void Foam::motionSmootherAlgo::setDisplacementPatchFields
(
    const labelList& patchIDs,
    pointVectorField& displacement
)
{
    // Adapt the fixedValue bc's (i.e. copy internal point data to
    // boundaryField for all affected patches)
    forAll(patchIDs, i)
    {
        label patchI = patchIDs[i];

        displacement.boundaryField()[patchI] ==
            displacement.boundaryField()[patchI].patchInternalField();
    }

    // Make consistent with non-adapted bc's by evaluating those now and
    // resetting the displacement from the values.
    // Note that we're just doing a correctBoundaryConditions with
    // fixedValue bc's first.
    labelHashSet adaptPatchSet(patchIDs);

    const lduSchedule& patchSchedule = displacement.mesh().globalData().
        patchSchedule();

    forAll(patchSchedule, patchEvalI)
    {
        label patchI = patchSchedule[patchEvalI].patch;

        if (!adaptPatchSet.found(patchI))
        {
            if (patchSchedule[patchEvalI].init)
            {
                displacement.boundaryField()[patchI]
                    .initEvaluate(Pstream::scheduled);
            }
            else
            {
                displacement.boundaryField()[patchI]
                    .evaluate(Pstream::scheduled);
            }
        }
    }

    // Multi-patch constraints
    pointConstraints::New(displacement.mesh()).constrainCorners(displacement);

    // Adapt the fixedValue bc's (i.e. copy internal point data to
    // boundaryField for all affected patches) to take the changes caused
    // by multi-corner constraints into account.
    forAll(patchIDs, i)
    {
        label patchI = patchIDs[i];

        displacement.boundaryField()[patchI] ==
            displacement.boundaryField()[patchI].patchInternalField();
    }
}


void Foam::motionSmootherAlgo::setDisplacementPatchFields()
{
    setDisplacementPatchFields(adaptPatchIDs_, displacement_);
}


void Foam::motionSmootherAlgo::setDisplacement
(
    const labelList& patchIDs,
    const indirectPrimitivePatch& pp,
    pointField& patchDisp,
    pointVectorField& displacement
)
{
    const polyMesh& mesh = displacement.mesh()();

    // See comment in .H file about shared points.
    // We want to disallow effect of loose coupled points - we only
    // want to see effect of proper fixedValue boundary conditions

    const labelList& cppMeshPoints =
        mesh.globalData().coupledPatch().meshPoints();

    forAll(cppMeshPoints, i)
    {
        displacement[cppMeshPoints[i]] = vector::zero;
    }

    const labelList& ppMeshPoints = pp.meshPoints();

    // Set internal point data from displacement on combined patch points.
    forAll(ppMeshPoints, patchPointI)
    {
        displacement[ppMeshPoints[patchPointI]] = patchDisp[patchPointI];
    }


    // Combine any coupled points
    syncTools::syncPointList
    (
        mesh,
        displacement,
        maxMagEqOp(),           // combine op
        vector::zero            // null value
    );


    // Adapt the fixedValue bc's (i.e. copy internal point data to
    // boundaryField for all affected patches)
    setDisplacementPatchFields(patchIDs, displacement);


    if (debug)
    {
        OFstream str(mesh.db().path()/"changedPoints.obj");
        label nVerts = 0;
        forAll(ppMeshPoints, patchPointI)
        {
            const vector& newDisp = displacement[ppMeshPoints[patchPointI]];

            if (mag(newDisp-patchDisp[patchPointI]) > SMALL)
            {
                const point& pt = mesh.points()[ppMeshPoints[patchPointI]];

                meshTools::writeOBJ(str, pt);
                nVerts++;
                //Pout<< "Point:" << pt
                //    << " oldDisp:" << patchDisp[patchPointI]
                //    << " newDisp:" << newDisp << endl;
            }
        }
        Pout<< "Written " << nVerts << " points that are changed to file "
            << str.name() << endl;
    }

    // Now reset input displacement
    forAll(ppMeshPoints, patchPointI)
    {
        patchDisp[patchPointI] = displacement[ppMeshPoints[patchPointI]];
    }
}


void Foam::motionSmootherAlgo::setDisplacement(pointField& patchDisp)
{
    setDisplacement(adaptPatchIDs_, pp_, patchDisp, displacement_);
}


// correctBoundaryConditions with fixedValue bc's first.
void Foam::motionSmootherAlgo::correctBoundaryConditions
(
    pointVectorField& displacement
) const
{
    labelHashSet adaptPatchSet(adaptPatchIDs_);

    const lduSchedule& patchSchedule = mesh_.globalData().patchSchedule();

    // 1. evaluate on adaptPatches
    forAll(patchSchedule, patchEvalI)
    {
        label patchI = patchSchedule[patchEvalI].patch;

        if (adaptPatchSet.found(patchI))
        {
            if (patchSchedule[patchEvalI].init)
            {
                displacement.boundaryField()[patchI]
                    .initEvaluate(Pstream::blocking);
            }
            else
            {
                displacement.boundaryField()[patchI]
                    .evaluate(Pstream::blocking);
            }
        }
    }


    // 2. evaluate on non-AdaptPatches
    forAll(patchSchedule, patchEvalI)
    {
        label patchI = patchSchedule[patchEvalI].patch;

        if (!adaptPatchSet.found(patchI))
        {
            if (patchSchedule[patchEvalI].init)
            {
                displacement.boundaryField()[patchI]
                    .initEvaluate(Pstream::blocking);
            }
            else
            {
                displacement.boundaryField()[patchI]
                    .evaluate(Pstream::blocking);
            }
        }
    }

    // Multi-patch constraints
    pointConstraints::New(displacement.mesh()).constrainCorners(displacement);

    // Correct for problems introduced by corner constraints
    syncTools::syncPointList
    (
        mesh_,
        displacement,
        maxMagEqOp(),           // combine op
        vector::zero            // null value
    );
}



void Foam::motionSmootherAlgo::modifyMotionPoints(pointField& newPoints) const
{
    // Correct for 2-D motion
    const twoDPointCorrector& twoDCorrector = twoDPointCorrector::New(mesh_);

    if (twoDCorrector.required())
    {
        Info<< "Correcting 2-D mesh motion";

        if (mesh_.globalData().parallel())
        {
            WarningIn("motionSmootherAlgo::modifyMotionPoints(pointField&)")
                << "2D mesh-motion probably not correct in parallel" << endl;
        }

        // We do not want to move 3D planes so project all points onto those
        const pointField& oldPoints = mesh_.points();
        const edgeList& edges = mesh_.edges();

        const labelList& neIndices = twoDCorrector.normalEdgeIndices();
        const vector& pn = twoDCorrector.planeNormal();

        forAll(neIndices, i)
        {
            const edge& e = edges[neIndices[i]];

            point& pStart = newPoints[e.start()];

            pStart += pn*(pn & (oldPoints[e.start()] - pStart));

            point& pEnd = newPoints[e.end()];

            pEnd += pn*(pn & (oldPoints[e.end()] - pEnd));
        }

        // Correct tangentially
        twoDCorrector.correctPoints(newPoints);
        Info<< " ...done" << endl;
    }

    if (debug)
    {
        Pout<< "motionSmootherAlgo::modifyMotionPoints :"
            << " testing sync of newPoints."
            << endl;
        testSyncPositions(newPoints, 1e-6*mesh_.bounds().mag());
    }
}


void Foam::motionSmootherAlgo::movePoints()
{
    // Make sure to clear out tetPtIs since used in checks (sometimes, should
    // really check)
    mesh_.clearAdditionalGeom();
    pp_.movePoints(mesh_.points());
}


Foam::scalar Foam::motionSmootherAlgo::setErrorReduction
(
    const scalar errorReduction
)
{
    scalar oldErrorReduction = readScalar(paramDict_.lookup("errorReduction"));

    paramDict_.remove("errorReduction");
    paramDict_.add("errorReduction", errorReduction);

    return oldErrorReduction;
}


bool Foam::motionSmootherAlgo::scaleMesh
(
    labelList& checkFaces,
    const bool smoothMesh,
    const label nAllowableErrors
)
{
    List<labelPair> emptyBaffles;
    return scaleMesh
    (
        checkFaces,
        emptyBaffles,
        smoothMesh,
        nAllowableErrors
    );
}


bool Foam::motionSmootherAlgo::scaleMesh
(
    labelList& checkFaces,
    const List<labelPair>& baffles,
    const bool smoothMesh,
    const label nAllowableErrors
)
{
    return scaleMesh
    (
        checkFaces,
        baffles,
        paramDict_,
        paramDict_,
        smoothMesh,
        nAllowableErrors
    );
}


Foam::tmp<Foam::pointField> Foam::motionSmootherAlgo::curPoints() const
{
    // Set newPoints as old + scale*displacement

    // Create overall displacement with same b.c.s as displacement_
    wordList actualPatchTypes;
    {
        const pointBoundaryMesh& pbm = displacement_.mesh().boundary();
        actualPatchTypes.setSize(pbm.size());
        forAll(pbm, patchI)
        {
            actualPatchTypes[patchI] = pbm[patchI].type();
        }
    }

    wordList actualPatchFieldTypes;
    {
        const pointVectorField::GeometricBoundaryField& pfld =
            displacement_.boundaryField();
        actualPatchFieldTypes.setSize(pfld.size());
        forAll(pfld, patchI)
        {
            if (isA<fixedValuePointPatchField<vector> >(pfld[patchI]))
            {
                // Get rid of funny
                actualPatchFieldTypes[patchI] =
                    fixedValuePointPatchField<vector>::typeName;
            }
            else
            {
                actualPatchFieldTypes[patchI] = pfld[patchI].type();
            }
        }
    }

    pointVectorField totalDisplacement
    (
        IOobject
        (
            "totalDisplacement",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        scale_*displacement_,
        actualPatchFieldTypes,
        actualPatchTypes
    );
    correctBoundaryConditions(totalDisplacement);

    if (debug)
    {
        Pout<< "scaleMesh : testing sync of totalDisplacement" << endl;
        testSyncField
        (
            totalDisplacement,
            maxMagEqOp(),
            vector::zero,   // null value
            1e-6*mesh_.bounds().mag()
        );
    }

    tmp<pointField> tnewPoints(oldPoints_ + totalDisplacement.internalField());

    // Correct for 2-D motion
    modifyMotionPoints(tnewPoints());

    return tnewPoints;
}


bool Foam::motionSmootherAlgo::scaleMesh
(
    labelList& checkFaces,
    const List<labelPair>& baffles,
    const dictionary& paramDict,
    const dictionary& meshQualityDict,
    const bool smoothMesh,
    const label nAllowableErrors
)
{
    if (!smoothMesh && adaptPatchIDs_.empty())
    {
        FatalErrorIn
        (
            "motionSmootherAlgo::scaleMesh"
            "("
                "labelList&, "
                "const List<labelPair>&, "
                "const dictionary&, "
                "const dictionary&, "
                "const bool, "
                "const label"
            ")"
        )
            << "You specified both no movement on the internal mesh points"
            << " (smoothMesh = false)" << nl
            << "and no movement on the patch (adaptPatchIDs is empty)" << nl
            << "Hence nothing to adapt."
            << exit(FatalError);
    }

    if (debug)
    {
        // Had a problem with patches moved non-synced. Check transformations.
        const polyBoundaryMesh& patches = mesh_.boundaryMesh();

        Pout<< "Entering scaleMesh : coupled patches:" << endl;
        forAll(patches, patchI)
        {
            if (patches[patchI].coupled())
            {
                const coupledPolyPatch& pp =
                    refCast<const coupledPolyPatch>(patches[patchI]);

                Pout<< '\t' << patchI << '\t' << pp.name()
                    << " parallel:" << pp.parallel()
                    << " separated:" << pp.separated()
                    << " forwardT:" << pp.forwardT().size()
                    << endl;
            }
        }
    }

    const scalar errorReduction =
        readScalar(paramDict.lookup("errorReduction"));
    const label nSmoothScale =
        readLabel(paramDict.lookup("nSmoothScale"));


    // Note: displacement_ should already be synced already from setDisplacement
    // but just to make sure.
    syncTools::syncPointList
    (
        mesh_,
        displacement_,
        maxMagEqOp(),
        vector::zero    // null value
    );

    Info<< "Moving mesh using displacement scaling :"
        << " min:" << gMin(scale_.internalField())
        << "  max:" << gMax(scale_.internalField())
        << endl;

    // Get points using current displacement and scale. Optionally 2D corrected.
    pointField newPoints(curPoints());

    // Move. No need to do 2D correction since curPoints already done that.
    mesh_.movePoints(newPoints);
    movePoints();


    // Check. Returns parallel number of incorrect faces.
    faceSet wrongFaces(mesh_, "wrongFaces", mesh_.nFaces()/100+100);
    checkMesh(false, mesh_, meshQualityDict, checkFaces, baffles, wrongFaces);

    if (returnReduce(wrongFaces.size(), sumOp<label>()) <= nAllowableErrors)
    {
        return true;
    }
    else
    {
        // Sync across coupled faces by extending the set.
        wrongFaces.sync(mesh_);

        // Special case:
        // if errorReduction is set to zero, extend wrongFaces
        // to face-Cell-faces to ensure quick return to previously valid mesh

        if (mag(errorReduction) < SMALL)
        {
            labelHashSet newWrongFaces(wrongFaces);
            forAllConstIter(labelHashSet, wrongFaces, iter)
            {
                label own = mesh_.faceOwner()[iter.key()];
                const cell& ownFaces = mesh_.cells()[own];

                forAll(ownFaces, cfI)
                {
                    newWrongFaces.insert(ownFaces[cfI]);
                }

                if (iter.key() < mesh_.nInternalFaces())
                {
                    label nei = mesh_.faceNeighbour()[iter.key()];
                    const cell& neiFaces = mesh_.cells()[nei];

                    forAll(neiFaces, cfI)
                    {
                        newWrongFaces.insert(neiFaces[cfI]);
                    }
                }
            }
            wrongFaces.transfer(newWrongFaces);
            wrongFaces.sync(mesh_);
        }


        // Find out points used by wrong faces and scale displacement.
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        pointSet usedPoints(mesh_, "usedPoints", getPoints(wrongFaces));
        usedPoints.sync(mesh_);



        // Grow a few layers to determine
        // - points to be smoothed
        // - faces to be checked in next iteration
        PackedBoolList isAffectedPoint(mesh_.nPoints());
        getAffectedFacesAndPoints
        (
            nSmoothScale,       // smoothing iterations
            wrongFaces,         // error faces
            checkFaces,
            isAffectedPoint
        );

        if (debug)
        {
            Pout<< "Faces in error:" << wrongFaces.size()
                << "  with points:" << usedPoints.size()
                << endl;
        }

        if (adaptPatchIDs_.size())
        {
            // Scale conflicting patch points
            scaleField(pp_.meshPoints(), usedPoints, errorReduction, scale_);
            //subtractField(pp_.meshPoints(), usedPoints, 0.1, scale_);
        }
        if (smoothMesh)
        {
            // Scale conflicting internal points
            scaleField(usedPoints, errorReduction, scale_);
            //subtractField(usedPoints, 0.1, scale_);
        }

        scalarField eWeights(calcEdgeWeights(oldPoints_));

        for (label i = 0; i < nSmoothScale; i++)
        {
            if (adaptPatchIDs_.size())
            {
                // Smooth patch values
                pointScalarField oldScale(scale_);
                minSmooth
                (
                    eWeights,
                    isAffectedPoint,
                    pp_.meshPoints(),
                    oldScale,
                    scale_
                );
                checkFld(scale_);
            }
            if (smoothMesh)
            {
                // Smooth internal values
                pointScalarField oldScale(scale_);
                minSmooth(eWeights, isAffectedPoint, oldScale, scale_);
                checkFld(scale_);
            }
        }

        syncTools::syncPointList
        (
            mesh_,
            scale_,
            maxEqOp<scalar>(),
            -GREAT              // null value
        );


        if (debug)
        {
            Pout<< "scale_ after smoothing :"
                << " min:" << Foam::gMin(scale_)
                << " max:" << Foam::gMax(scale_)
                << endl;
        }

        return false;
    }
}


void Foam::motionSmootherAlgo::updateMesh()
{
    const pointBoundaryMesh& patches = pMesh_.boundary();

    // Check whether displacement has fixed value b.c. on adaptPatchID
    forAll(adaptPatchIDs_, i)
    {
        label patchI = adaptPatchIDs_[i];

        if
        (
           !isA<fixedValuePointPatchVectorField>
            (
                displacement_.boundaryField()[patchI]
            )
        )
        {
            FatalErrorIn
            (
                "motionSmootherAlgo::updateMesh"
            )   << "Patch " << patches[patchI].name()
                << " has wrong boundary condition "
                << displacement_.boundaryField()[patchI].type()
                << " on field " << displacement_.name() << nl
                << "Only type allowed is "
                << fixedValuePointPatchVectorField::typeName
                << exit(FatalError);
        }
    }


    // Determine internal points. Note that for twoD there are no internal
    // points so we use the points of adaptPatchIDs instead

    const labelList& meshPoints = pp_.meshPoints();

    forAll(meshPoints, i)
    {
        isInternalPoint_.unset(meshPoints[i]);
    }

    // Calculate master edge addressing
    isMasterEdge_ = syncTools::getMasterEdges(mesh_);
}


// ************************************************************************* //
