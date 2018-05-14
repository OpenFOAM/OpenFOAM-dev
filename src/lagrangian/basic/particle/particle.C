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

#include "particle.H"
#include "transform.H"
#include "treeDataCell.H"
#include "cubicEqn.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::scalar Foam::particle::negativeSpaceDisplacementFactor = 1.01;

Foam::label Foam::particle::particleCount_ = 0;

namespace Foam
{
    defineTypeNameAndDebug(particle, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::particle::stationaryTetReverseTransform
(
    vector& centre,
    scalar& detA,
    barycentricTensor& T
) const
{
    barycentricTensor A = stationaryTetTransform();

    const vector ab = A.b() - A.a();
    const vector ac = A.c() - A.a();
    const vector ad = A.d() - A.a();
    const vector bc = A.c() - A.b();
    const vector bd = A.d() - A.b();

    centre = A.a();

    detA = ab & (ac ^ ad);

    T = barycentricTensor
    (
        bd ^ bc,
        ac ^ ad,
        ad ^ ab,
        ab ^ ac
    );
}


void Foam::particle::movingTetReverseTransform
(
    const scalar fraction,
    Pair<vector>& centre,
    FixedList<scalar, 4>& detA,
    FixedList<barycentricTensor, 3>& T
) const
{
    Pair<barycentricTensor> A = movingTetTransform(fraction);

    const Pair<vector> ab(A[0].b() - A[0].a(), A[1].b() - A[1].a());
    const Pair<vector> ac(A[0].c() - A[0].a(), A[1].c() - A[1].a());
    const Pair<vector> ad(A[0].d() - A[0].a(), A[1].d() - A[1].a());
    const Pair<vector> bc(A[0].c() - A[0].b(), A[1].c() - A[1].b());
    const Pair<vector> bd(A[0].d() - A[0].b(), A[1].d() - A[1].b());

    centre[0] = A[0].a();
    centre[1] = A[1].a();

    detA[0] = ab[0] & (ac[0] ^ ad[0]);
    detA[1] =
        (ab[1] & (ac[0] ^ ad[0]))
      + (ab[0] & (ac[1] ^ ad[0]))
      + (ab[0] & (ac[0] ^ ad[1]));
    detA[2] =
        (ab[0] & (ac[1] ^ ad[1]))
      + (ab[1] & (ac[0] ^ ad[1]))
      + (ab[1] & (ac[1] ^ ad[0]));
    detA[3] = ab[1] & (ac[1] ^ ad[1]);

    T[0] = barycentricTensor
    (
        bd[0] ^ bc[0],
        ac[0] ^ ad[0],
        ad[0] ^ ab[0],
        ab[0] ^ ac[0]
    );
    T[1] = barycentricTensor
    (
        (bd[0] ^ bc[1]) + (bd[1] ^ bc[0]),
        (ac[0] ^ ad[1]) + (ac[1] ^ ad[0]),
        (ad[0] ^ ab[1]) + (ad[1] ^ ab[0]),
        (ab[0] ^ ac[1]) + (ab[1] ^ ac[0])
    );
    T[2] = barycentricTensor
    (
        bd[1] ^ bc[1],
        ac[1] ^ ad[1],
        ad[1] ^ ab[1],
        ab[1] ^ ac[1]
    );
}


void Foam::particle::reflect()
{
    Swap(coordinates_.c(), coordinates_.d());
}


void Foam::particle::rotate(const bool reverse)
{
    if (!reverse)
    {
        scalar temp = coordinates_.b();
        coordinates_.b() = coordinates_.c();
        coordinates_.c() = coordinates_.d();
        coordinates_.d() = temp;
    }
    else
    {
        scalar temp = coordinates_.d();
        coordinates_.d() = coordinates_.c();
        coordinates_.c() = coordinates_.b();
        coordinates_.b() = temp;
    }
}


void Foam::particle::changeTet(const label tetTriI)
{
    const bool isOwner = mesh_.faceOwner()[tetFacei_] == celli_;

    const label firstTetPtI = 1;
    const label lastTetPtI = mesh_.faces()[tetFacei_].size() - 2;

    if (tetTriI == 1)
    {
        changeFace(tetTriI);
    }
    else if (tetTriI == 2)
    {
        if (isOwner)
        {
            if (tetPti_ == lastTetPtI)
            {
                changeFace(tetTriI);
            }
            else
            {
                reflect();
                tetPti_ += 1;
            }
        }
        else
        {
            if (tetPti_ == firstTetPtI)
            {
                changeFace(tetTriI);
            }
            else
            {
                reflect();
                tetPti_ -= 1;
            }
        }
    }
    else if (tetTriI == 3)
    {
        if (isOwner)
        {
            if (tetPti_ == firstTetPtI)
            {
                changeFace(tetTriI);
            }
            else
            {
                reflect();
                tetPti_ -= 1;
            }
        }
        else
        {
            if (tetPti_ == lastTetPtI)
            {
                changeFace(tetTriI);
            }
            else
            {
                reflect();
                tetPti_ += 1;
            }
        }
    }
    else
    {
        FatalErrorInFunction
            << "Changing tet without changing cell should only happen when the "
            << "track is on triangle 1, 2 or 3."
            << exit(FatalError);
    }
}


void Foam::particle::changeFace(const label tetTriI)
{
    // Get the old topology
    const triFace triOldIs(currentTetIndices().faceTriIs(mesh_));

    // Get the shared edge and the pre-rotation
    edge sharedEdge;
    if (tetTriI == 1)
    {
        sharedEdge = edge(triOldIs[1], triOldIs[2]);
    }
    else if (tetTriI == 2)
    {
        sharedEdge = edge(triOldIs[2], triOldIs[0]);
    }
    else if (tetTriI == 3)
    {
        sharedEdge = edge(triOldIs[0], triOldIs[1]);
    }
    else
    {
        FatalErrorInFunction
            << "Changing face without changing cell should only happen when the"
            << " track is on triangle 1, 2 or 3."
            << exit(FatalError);

        sharedEdge = edge(-1, -1);
    }

    // Find the face in the same cell that shares the edge, and the
    // corresponding tetrahedra point
    tetPti_ = -1;
    forAll(mesh_.cells()[celli_], cellFaceI)
    {
        const label newFaceI = mesh_.cells()[celli_][cellFaceI];
        const class face& newFace = mesh_.faces()[newFaceI];
        const label newOwner = mesh_.faceOwner()[newFaceI];

        // Exclude the current face
        if (tetFacei_ == newFaceI)
        {
            continue;
        }

        // Loop over the edges, looking for the shared one. Note that we have to
        // match the direction of the edge as well as the end points in order to
        // avoid false positives when dealing with coincident ACMI faces.
        const label edgeComp = newOwner == celli_ ? -1 : +1;
        label edgeI = 0;
        for
        (
            ;
            edgeI < newFace.size()
         && edge::compare(sharedEdge, newFace.faceEdge(edgeI)) != edgeComp;
            ++ edgeI
        );

        // If the face does not contain the edge, then move on to the next face
        if (edgeI >= newFace.size())
        {
            continue;
        }

        // Make the edge index relative to the base point
        const label newBaseI = max(0, mesh_.tetBasePtIs()[newFaceI]);
        edgeI = (edgeI - newBaseI + newFace.size()) % newFace.size();

        // If the edge is next the base point (i.e., the index is 0 or n - 1),
        // then we swap it for the adjacent edge. This new edge is opposite the
        // base point, and defines the tet with the original edge in it.
        edgeI = min(max(1, edgeI), newFace.size() - 2);

        // Set the new face and tet point
        tetFacei_ = newFaceI;
        tetPti_ = edgeI;

        // Exit the loop now that the tet point has been found
        break;
    }

    if (tetPti_ == -1)
    {
        FatalErrorInFunction
            << "The search for an edge-connected face and tet-point failed."
            << exit(FatalError);
    }

    // Pre-rotation puts the shared edge opposite the base of the tetrahedron
    if (sharedEdge.otherVertex(triOldIs[1]) == -1)
    {
        rotate(false);
    }
    else if (sharedEdge.otherVertex(triOldIs[2]) == -1)
    {
        rotate(true);
    }

    // Get the new topology
    const triFace triNewIs = currentTetIndices().faceTriIs(mesh_);

    // Reflect to account for the change of triangle orientation on the new face
    reflect();

    // Post rotation puts the shared edge back in the correct location
    if (sharedEdge.otherVertex(triNewIs[1]) == -1)
    {
        rotate(true);
    }
    else if (sharedEdge.otherVertex(triNewIs[2]) == -1)
    {
        rotate(false);
    }
}


void Foam::particle::changeCell()
{
    // Set the cell to be the one on the other side of the face
    const label ownerCellI = mesh_.faceOwner()[tetFacei_];
    const bool isOwner = celli_ == ownerCellI;
    celli_ = isOwner ? mesh_.faceNeighbour()[tetFacei_] : ownerCellI;

    // Reflect to account for the change of triangle orientation in the new cell
    reflect();
}


void Foam::particle::changeToMasterPatch()
{
    label thisPatch = patch();

    forAll(mesh_.cells()[celli_], cellFaceI)
    {
        // Skip the current face and any internal faces
        const label otherFaceI = mesh_.cells()[celli_][cellFaceI];
        if (facei_ == otherFaceI || mesh_.isInternalFace(otherFaceI))
        {
            continue;
        }

        // Compare the two faces. If they are the same, chose the one with the
        // lower patch index. In the case of an ACMI-wall pair, this will be
        // the ACMI, which is what we want.
        const class face& thisFace = mesh_.faces()[facei_];
        const class face& otherFace = mesh_.faces()[otherFaceI];
        if (face::compare(thisFace, otherFace) != 0)
        {
            const label otherPatch =
                mesh_.boundaryMesh().whichPatch(otherFaceI);
            if (thisPatch > otherPatch)
            {
                facei_ = otherFaceI;
                thisPatch = otherPatch;
            }
        }
    }

    tetFacei_ = facei_;
}


void Foam::particle::locate
(
    const vector& position,
    const vector* direction,
    const label celli,
    const bool boundaryFail,
    const string boundaryMsg
)
{
    celli_ = celli;

    // Find the cell, if it has not been given
    if (celli_ < 0)
    {
        celli_ = mesh_.cellTree().findInside(position);
    }
    if (celli_ < 0)
    {
        FatalErrorInFunction
            << "Cell not found for particle position " << position << "."
            << exit(FatalError);
    }

    // Put the particle at the cell centre and in a random tet
    coordinates_ = barycentric(1, 0, 0, 0);
    tetFacei_ = mesh_.cells()[celli_][0];
    tetPti_ = 1;
    facei_ = -1;

    // Track to the injection point
    track(position - mesh_.cellCentres()[celli_], 0);
    if (!onFace())
    {
        return;
    }

    // We hit a boundary ...
    if (boundaryFail)
    {
        FatalErrorInFunction << boundaryMsg << exit(FatalError);
    }
    else
    {
        // Re-do the track, but this time do the bit tangential to the
        // direction/patch first. This gets us as close as possible to the
        // original path/position.

        if (direction == nullptr)
        {
            const polyPatch& p = mesh_.boundaryMesh()[patch()];
            direction = &p.faceNormals()[p.whichFace(facei_)];
        }

        const vector n = *direction/mag(*direction);
        const vector s = position - mesh_.cellCentres()[celli_];
        const vector sN = (s & n)*n;
        const vector sT = s - sN;

        coordinates_ = barycentric(1, 0, 0, 0);
        tetFacei_ = mesh_.cells()[celli_][0];
        tetPti_ = 1;
        facei_ = -1;

        track(sT, 0);
        track(sN, 0);

        static label nWarnings = 0;
        static const label maxNWarnings = 100;
        if (nWarnings < maxNWarnings)
        {
            WarningInFunction << boundaryMsg << endl;
            ++ nWarnings;
        }
        if (nWarnings == maxNWarnings)
        {
            WarningInFunction
                << "Suppressing any further warnings about particles being "
                << "located outside of the mesh." << endl;
            ++ nWarnings;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::particle::particle
(
    const polyMesh& mesh,
    const barycentric& coordinates,
    const label celli,
    const label tetFacei,
    const label tetPti
)
:
    mesh_(mesh),
    coordinates_(coordinates),
    celli_(celli),
    tetFacei_(tetFacei),
    tetPti_(tetPti),
    facei_(-1),
    stepFraction_(0.0),
    origProc_(Pstream::myProcNo()),
    origId_(getNewParticleID())
{}


Foam::particle::particle
(
    const polyMesh& mesh,
    const vector& position,
    const label celli
)
:
    mesh_(mesh),
    coordinates_(- vGreat, - vGreat, - vGreat, - vGreat),
    celli_(celli),
    tetFacei_(-1),
    tetPti_(-1),
    facei_(-1),
    stepFraction_(0.0),
    origProc_(Pstream::myProcNo()),
    origId_(getNewParticleID())
{
    locate
    (
        position,
        nullptr,
        celli,
        false,
        "Particle initialised with a location outside of the mesh."
    );
}


Foam::particle::particle(const particle& p)
:
    mesh_(p.mesh_),
    coordinates_(p.coordinates_),
    celli_(p.celli_),
    tetFacei_(p.tetFacei_),
    tetPti_(p.tetPti_),
    facei_(p.facei_),
    stepFraction_(p.stepFraction_),
    origProc_(p.origProc_),
    origId_(p.origId_)
{}


Foam::particle::particle(const particle& p, const polyMesh& mesh)
:
    mesh_(mesh),
    coordinates_(p.coordinates_),
    celli_(p.celli_),
    tetFacei_(p.tetFacei_),
    tetPti_(p.tetPti_),
    facei_(p.facei_),
    stepFraction_(p.stepFraction_),
    origProc_(p.origProc_),
    origId_(p.origId_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::particle::track
(
    const vector& displacement,
    const scalar fraction
)
{
    scalar f = trackToFace(displacement, fraction);

    while (onInternalFace())
    {
        changeCell();

        f *= trackToFace(f*displacement, f*fraction);
    }

    return f;
}


Foam::scalar Foam::particle::trackToCell
(
    const vector& displacement,
    const scalar fraction
)
{
    const scalar f = trackToFace(displacement, fraction);

    if (onInternalFace())
    {
        changeCell();
    }

    return f;
}


Foam::scalar Foam::particle::trackToFace
(
    const vector& displacement,
    const scalar fraction
)
{
    scalar f = 1;

    label tetTriI = onFace() ? 0 : -1;

    facei_ = -1;

    while (true)
    {
        f *= trackToTri(f*displacement, f*fraction, tetTriI);

        if (tetTriI == -1)
        {
            // The track has completed within the current tet
            return 0;
        }
        else if (tetTriI == 0)
        {
            // The track has hit a face, so set the current face and return
            facei_ = tetFacei_;
            return f;
        }
        else
        {
            // Move to the next tet and continue the track
            changeTet(tetTriI);
        }
    }
}


Foam::scalar Foam::particle::trackToStationaryTri
(
    const vector& displacement,
    const scalar fraction,
    label& tetTriI
)
{
    const vector x0 = position();
    const vector x1 = displacement;
    const barycentric y0 = coordinates_;

    if (debug)
    {
        Info<< "Particle " << origId() << endl << "Tracking from " << x0
            << " along " << x1 << " to " << x0 + x1 << endl;
    }

    // Get the tet geometry
    vector centre;
    scalar detA;
    barycentricTensor T;
    stationaryTetReverseTransform(centre, detA, T);

    if (debug)
    {
        vector o, b, v1, v2;
        stationaryTetGeometry(o, b, v1, v2);
        Info<< "Tet points o=" << o << ", b=" << b
            << ", v1=" << v1 << ", v2=" << v2 << endl
            << "Tet determinant = " << detA << endl
            << "Start local coordinates = " << y0 << endl;
    }

    // Get the factor by which the displacement is increased
    const scalar f = detA >= 0 ? 1 : negativeSpaceDisplacementFactor;

    // Calculate the local tracking displacement
    barycentric Tx1(f*x1 & T);

    if (debug)
    {
        Info<< "Local displacement = " << Tx1 << "/" << detA << endl;
    }

    // Calculate the hit fraction
    label iH = -1;
    scalar muH = std::isnormal(detA) && detA <= 0 ? vGreat : 1/detA;
    for (label i = 0; i < 4; ++ i)
    {
        if (std::isnormal(Tx1[i]) && Tx1[i] < 0)
        {
            scalar mu = - y0[i]/Tx1[i];

            if (debug)
            {
                Info<< "Hit on tet face " << i << " at local coordinate "
                    << y0 + mu*Tx1 << ", " << mu*detA*100 << "% of the "
                    << "way along the track" << endl;
            }

            if (0 <= mu && mu < muH)
            {
                iH = i;
                muH = mu;
            }
        }
    }

    // Set the new coordinates
    barycentric yH = y0 + muH*Tx1;

    // Clamp to zero any negative coordinates generated by round-off error
    for (label i = 0; i < 4; ++ i)
    {
        yH.replace(i, i == iH ? 0 : max(0, yH[i]));
    }

    // Re-normalise if within the tet
    if (iH == -1)
    {
        yH /= cmptSum(yH);
    }

    // Set the new position and hit index
    coordinates_ = yH;
    tetTriI = iH;

    if (debug)
    {
        if (iH != -1)
        {
            Info<< "Track hit tet face " << iH << " first" << endl;
        }
        else
        {
            Info<< "Track hit no tet faces" << endl;
        }
        Info<< "End local coordinates = " << yH << endl
            << "End global coordinates = " << position() << endl
            << "Tracking displacement = " << position() - x0 << endl
            << muH*detA*100 << "% of the step from " << stepFraction_ << " to "
            << stepFraction_ + fraction << " completed" << endl << endl;
    }

    // Set the proportion of the track that has been completed
    stepFraction_ += fraction*muH*detA;

    return iH != -1 ? 1 - muH*detA : 0;
}


Foam::scalar Foam::particle::trackToMovingTri
(
    const vector& displacement,
    const scalar fraction,
    label& tetTriI
)
{
    const vector x0 = position();
    const vector x1 = displacement;
    const barycentric y0 = coordinates_;

    if (debug)
    {
        Info<< "Particle " << origId() << endl << "Tracking from " << x0
            << " along " << x1 << " to " << x0 + x1 << endl;
    }

    // Get the tet geometry
    Pair<vector> centre;
    FixedList<scalar, 4> detA;
    FixedList<barycentricTensor, 3> T;
    movingTetReverseTransform(fraction, centre, detA, T);

    // Get the factor by which the displacement is increased
    const scalar f = detA[0] >= 0 ? 1 : negativeSpaceDisplacementFactor;

    // Get the relative global position
    const vector x0Rel = x0 - centre[0];
    const vector x1Rel = f*x1 - centre[1];

    // Form the determinant and hit equations
    cubicEqn detAEqn(sqr(detA[0])*detA[3], detA[0]*detA[2], detA[1], 1);
    const barycentric yC(1, 0, 0, 0);
    const barycentric hitEqnA =
        ((x1Rel & T[2])                  + detA[3]*yC)*sqr(detA[0]);
    const barycentric hitEqnB =
        ((x1Rel & T[1]) + (x0Rel & T[2]) + detA[2]*yC)*detA[0];
    const barycentric hitEqnC =
        ((x1Rel & T[0]) + (x0Rel & T[1]) + detA[1]*yC);
    const barycentric hitEqnD = y0;
    FixedList<cubicEqn, 4> hitEqn;
    forAll(hitEqn, i)
    {
        hitEqn[i] = cubicEqn(hitEqnA[i], hitEqnB[i], hitEqnC[i], hitEqnD[i]);
    }

    // Calculate the hit fraction
    label iH = -1;
    scalar muH = std::isnormal(detA[0]) && detA[0] <= 0 ? vGreat : 1/detA[0];
    for (label i = 0; i < 4; ++ i)
    {
        const Roots<3> mu = hitEqn[i].roots();

        for (label j = 0; j < 3; ++ j)
        {
            if (mu.type(j) == roots::real && hitEqn[i].derivative(mu[j]) < 0)
            {
                if (0 <= mu[j] && mu[j] < muH)
                {
                    iH = i;
                    muH = mu[j];
                }
            }
        }
    }

    // Set the new coordinates
    barycentric yH
    (
        hitEqn[0].value(muH),
        hitEqn[1].value(muH),
        hitEqn[2].value(muH),
        hitEqn[3].value(muH)
    );
    // !!! <-- This fails if the tet collapses onto the particle, as detA tends
    // to zero at the hit. In this instance, we can differentiate the hit and
    // detA polynomials to find a limiting location, but this will not be on a
    // triangle. We will then need to do a second track through the degenerate
    // tet to find the final hit position. This second track is over zero
    // distance and therefore can be of the static mesh type. This has not yet
    // been implemented.
    const scalar detAH = detAEqn.value(muH);
    if (!std::isnormal(detAH))
    {
        FatalErrorInFunction
            << "A moving tet collapsed onto a particle. This is not supported. "
            << "The mesh is too poor, or the motion too severe, for particle "
            << "tracking to function." << exit(FatalError);
    }
    yH /= detAH;

    // Clamp to zero any negative coordinates generated by round-off error
    for (label i = 0; i < 4; ++ i)
    {
        yH.replace(i, i == iH ? 0 : max(0, yH[i]));
    }

    // Re-normalise if within the tet
    if (iH == -1)
    {
        yH /= cmptSum(yH);
    }

    // Set the new position and hit index
    coordinates_ = yH;
    tetTriI = iH;

    // Set the proportion of the track that has been completed
    stepFraction_ += fraction*muH*detA[0];

    if (debug)
    {
        if (iH != -1)
        {
            Info<< "Track hit tet face " << iH << " first" << endl;
        }
        else
        {
            Info<< "Track hit no tet faces" << endl;
        }
        Info<< "End local coordinates = " << yH << endl
            << "End global coordinates = " << position() << endl;
    }

    return iH != -1 ? 1 - muH*detA[0] : 0;
}


Foam::scalar Foam::particle::trackToTri
(
    const vector& displacement,
    const scalar fraction,
    label& tetTriI
)
{
    if (mesh_.moving())
    {
        return trackToMovingTri(displacement, fraction, tetTriI);
    }
    else
    {
        return trackToStationaryTri(displacement, fraction, tetTriI);
    }
}


Foam::vector Foam::particle::deviationFromMeshCentre() const
{
    if (cmptMin(mesh_.geometricD()) == -1)
    {
        vector pos = position(), posC = pos;
        meshTools::constrainToMeshCentre(mesh_, posC);
        return pos - posC;
    }
    else
    {
        return vector::zero;
    }
}


void Foam::particle::transformProperties(const tensor&)
{}


void Foam::particle::transformProperties(const vector&)
{}


void Foam::particle::prepareForParallelTransfer()
{
    // Convert the face index to be local to the processor patch
    facei_ = mesh_.boundaryMesh()[patch()].whichFace(facei_);
}


void Foam::particle::correctAfterParallelTransfer
(
    const label patchi,
    trackingData& td
)
{
    const coupledPolyPatch& ppp =
        refCast<const coupledPolyPatch>(mesh_.boundaryMesh()[patchi]);

    if (!ppp.parallel())
    {
        const tensor& T =
        (
            ppp.forwardT().size() == 1
          ? ppp.forwardT()[0]
          : ppp.forwardT()[facei_]
        );
        transformProperties(T);
    }
    else if (ppp.separated())
    {
        const vector& s =
        (
            (ppp.separation().size() == 1)
          ? ppp.separation()[0]
          : ppp.separation()[facei_]
        );
        transformProperties(-s);
    }

    // Set the topology
    celli_ = ppp.faceCells()[facei_];
    facei_ += ppp.start();
    tetFacei_ = facei_;
    // Faces either side of a coupled patch are numbered in opposite directions
    // as their normals both point away from their connected cells. The tet
    // point therefore counts in the opposite direction from the base point.
    tetPti_ = mesh_.faces()[tetFacei_].size() - 1 - tetPti_;

    // Reflect to account for the change of triangle orientation in the new cell
    reflect();

    // Note that the position does not need transforming explicitly. The face-
    // triangle on the receive patch is the transformation of the one on the
    // send patch, so whilst the barycentric coordinates remain the same, the
    // change of triangle implicitly transforms the position.
}


void Foam::particle::prepareForInteractionListReferral
(
    const vectorTensorTransform& transform
)
{
    // Get the transformed position
    const vector pos = transform.invTransformPosition(position());

    // Break the topology
    celli_ = -1;
    tetFacei_ = -1;
    tetPti_ = -1;
    facei_ = -1;

    // Store the position in the barycentric data
    coordinates_ = barycentric(1 - cmptSum(pos), pos.x(), pos.y(), pos.z());

    // Transform the properties
    transformProperties(- transform.t());
    if (transform.hasR())
    {
        transformProperties(transform.R().T());
    }
}


void Foam::particle::correctAfterInteractionListReferral(const label celli)
{
    // Get the position from the barycentric data
    const vector pos(coordinates_.b(), coordinates_.c(), coordinates_.d());

    // Create some arbitrary topology for the supplied cell
    celli_ = celli;
    tetFacei_ = mesh_.cells()[celli_][0];
    tetPti_ = 1;
    facei_ = -1;

    // Get the reverse transform and directly set the coordinates from the
    // position. This isn't likely to be correct; the particle is probably not
    // in this tet. It will, however, generate the correct vector when the
    // position method is called. A referred particle should never be tracked,
    // so this approximate topology is good enough. By using the nearby cell we
    // minimize the error associated with the incorrect topology.
    coordinates_ = barycentric(1, 0, 0, 0);
    if (mesh_.moving())
    {
        Pair<vector> centre;
        FixedList<scalar, 4> detA;
        FixedList<barycentricTensor, 3> T;
        movingTetReverseTransform(0, centre, detA, T);
        coordinates_ += (pos - centre[0]) & T[0]/detA[0];
    }
    else
    {
        vector centre;
        scalar detA;
        barycentricTensor T;
        stationaryTetReverseTransform(centre, detA, T);
        coordinates_ += (pos - centre) & T/detA;
    }
}


Foam::label Foam::particle::procTetPt
(
    const polyMesh& procMesh,
    const label procCell,
    const label procTetFace
) const
{
    // The tet point on the procMesh differs from the current tet point if the
    // mesh and procMesh faces are of differing orientation. The change is the
    // same as in particle::correctAfterParallelTransfer.

    if
    (
        (mesh_.faceOwner()[tetFacei_] == celli_)
     == (procMesh.faceOwner()[procTetFace] == procCell)
    )
    {
        return tetPti_;
    }
    else
    {
        return procMesh.faces()[procTetFace].size() - 1 - tetPti_;
    }
}


void Foam::particle::autoMap
(
    const vector& position,
    const mapPolyMesh& mapper
)
{
    locate
    (
        position,
        nullptr,
        mapper.reverseCellMap()[celli_],
        true,
        "Particle mapped to a location outside of the mesh."
    );
}


// * * * * * * * * * * * * * * Friend Operators * * * * * * * * * * * * * * //

bool Foam::operator==(const particle& pA, const particle& pB)
{
    return (pA.origProc() == pB.origProc() && pA.origId() == pB.origId());
}


bool Foam::operator!=(const particle& pA, const particle& pB)
{
    return !(pA == pB);
}


// ************************************************************************* //
