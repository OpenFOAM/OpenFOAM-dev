/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024 OpenFOAM Foundation
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

#include "tracking.H"
#include "cubicEqn.H"
#include "indexedOctree.H"
#include "treeDataCell.H"

#include "emptyPolyPatch.H"
#include "wedgePolyPatch.H"
#include "cyclicPolyPatch.H"
#include "nonConformalCyclicPolyPatch.H"
#include "processorPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace tracking
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Tetrahedra

    //- Get the reverse transform associated with the current tet. The
    //  conversion is detA*y = (x - centre) & T. The variables x, y and
    //  centre have the same meaning as for the forward transform. T is
    //  the transposed inverse of the forward transform tensor, A,
    //  multiplied by its determinant, detA. This separation allows
    //  the barycentric tracking algorithm to function on inverted or
    //  degenerate tetrahedra.
    void stationaryTetReverseTransform
    (
        const polyMesh& mesh,
        const label celli,
        const label facei,
        const label faceTrii,
        vector& centre,
        scalar& detA,
        barycentricTensor& T
    );

    //- Get the reverse transformation associated with the current,
    //  moving, tet. This is of the same form as for the static case. As
    //  with the moving geometry, a function of the tracking fraction is
    //  returned for each component. The functions are higher order than
    //  for the forward transform; the determinant is cubic, and the
    //  tensor is quadratic.
    void movingTetReverseTransform
    (
        const polyMesh& mesh,
        const label celli,
        const label facei,
        const label faceTrii,
        const scalar startStepFraction,
        const scalar endStepFraction,
        Pair<vector>& centre,
        FixedList<scalar, 4>& detA,
        FixedList<barycentricTensor, 3>& T
    );


// Tracking

    //- The counter nTracksBehind is the number of tracks carried out that
    //  ended in a step fraction less than the maximum reached so far. Once
    //  this reaches maxNTracksBehind, tracking is abandoned for the current
    //  step.
    //
    //  This is needed because when tetrahedra are inverted a straight
    //  trajectory can form a closed loop through regions of overlapping
    //  positive and negative space. Without this break clause, such loops can
    //  result in a valid track which never ends.
    //
    //  Because the test is susceptible to round off error, a track of zero
    //  length will also increment the counter. As such, it is important that
    //  maxNTracksBehind is set large enough so that valid small tracks do not
    //  result in the track being abandoned. The largest number of valid small
    //  tracks that are likely to be performed sequentially is equal to the
    //  number of tetrahedra that can meet at a point. An estimate of this
    //  number is therefore used to set maxNTracksBehind.
    static const label maxNTracksBehind = 48;

    //- See toTri. For a stationary mesh.
    Tuple2<label, scalar> toStationaryTri
    (
        const polyMesh& mesh,
        const vector& displacement,
        const scalar fraction,
        barycentric& coordinates,
        label& celli,
        label& facei,
        label& faceTrii,
        scalar& stepFraction,
        scalar& stepFractionBehind,
        label& nTracksBehind,
        const string& debugPrefix = NullObjectRef<string>()
    );

    //- See toTri. For a moving mesh.
    Tuple2<label, scalar> toMovingTri
    (
        const polyMesh& mesh,
        const vector& displacement,
        const scalar fraction,
        barycentric& coordinates,
        label& celli,
        label& facei,
        label& faceTrii,
        scalar& stepFraction,
        scalar& stepFractionBehind,
        label& nTracksBehind,
        const string& debugPrefix = NullObjectRef<string>()
    );

    //- Track along the displacement for a given fraction of the overall
    //  time-step. End when the track is complete or when a tet triangle is
    //  hit. Return the index of the tet triangle that was hit, or -1 if the
    //  end position was reached. Also return the proportion of the
    //  displacement still to be completed.
    Tuple2<label, scalar> toTri
    (
        const polyMesh& mesh,
        const vector& displacement,
        const scalar fraction,
        barycentric& coordinates,
        label& celli,
        label& facei,
        label& faceTrii,
        scalar& stepFraction,
        scalar& stepFractionBehind,
        label& nTracksBehind,
        const string& debugPrefix = NullObjectRef<string>()
    );


// Transformations

    //- Reflection transform. Corrects the coordinates when the track moves
    //  between two tets which share a base vertex, but for which the other two
    //  non cell-centre vertices are reversed. All hits which retain the same
    //  face behave this way, as do face hits.
    void reflect(barycentric& coordinates);

    //- Rotation transform. Corrects the coordinates when the track moves
    //  between two tets with different base vertices, but are otherwise
    //  similarly oriented. Hits which change the face within the cell make use
    //  of both this and the reflect transform.
    void rotate(const bool reverse, barycentric& coordinates);


// Topology Changes

    //- Change face-triangle within a cell. Called after a tet-triangle is hit.
    void changeFaceTri
    (
        const polyMesh& mesh,
        const label tetTrii,
        barycentric& coordinates,
        const label celli,
        label& facei,
        label& faceTrii
    );

    //- Change face within a cell. Called (if necessary) by changeFaceTri.
    void changeFace
    (
        const polyMesh& mesh,
        const label tetTrii,
        barycentric& coordinates,
        const label celli,
        label& facei,
        label& faceTrii
    );

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace tracking
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::tracking::stationaryTetReverseTransform
(
    const polyMesh& mesh,
    const label celli,
    const label facei,
    const label faceTrii,
    vector& centre,
    scalar& detA,
    barycentricTensor& T
)
{
    barycentricTensor A = stationaryTetTransform(mesh, celli, facei, faceTrii);

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


void Foam::tracking::movingTetReverseTransform
(
    const polyMesh& mesh,
    const label celli,
    const label facei,
    const label faceTrii,
    const scalar startStepFraction,
    const scalar endStepFraction,
    Pair<vector>& centre,
    FixedList<scalar, 4>& detA,
    FixedList<barycentricTensor, 3>& T
)
{
    Pair<barycentricTensor> A =
        movingTetTransform
        (
            mesh,
            celli,
            facei,
            faceTrii,
            startStepFraction,
            endStepFraction
        );

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


Foam::Tuple2<Foam::label, Foam::scalar> Foam::tracking::toStationaryTri
(
    const polyMesh& mesh,
    const vector& displacement,
    const scalar fraction,
    barycentric& coordinates,
    label& celli,
    label& facei,
    label& faceTrii,
    scalar& stepFraction,
    scalar& stepFractionBehind,
    label& nTracksBehind,
    const string& debugPrefix_
)
{
    const bool debug = notNull(debugPrefix_);
    #define debugPrefix debugPrefix_.c_str() << ": "
    #define debugIndent string(debugPrefix_.size(), ' ').c_str() << ": "

    const vector x0 =
        position(mesh, coordinates, celli, facei, faceTrii, stepFraction);
    const vector x1 = displacement;
    const barycentric y0 = coordinates;

    DebugInfo
        << debugPrefix << "Tracking from " << x0
        << " along " << x1 << " to " << x0 + x1 << nl;

    // Get the tet geometry
    vector centre;
    scalar detA;
    barycentricTensor T;
    stationaryTetReverseTransform
    (
        mesh,
        celli,
        facei,
        faceTrii,
        centre,
        detA,
        T
    );

    if (debug)
    {
        vector o, b, v1, v2;
        stationaryTetGeometry
        (
            mesh,
            celli,
            facei,
            faceTrii,
            o,
            b,
            v1,
            v2
        );

        Info<< debugIndent << "Tet points o=" << o << ", b=" << b
            << ", v1=" << v1 << ", v2=" << v2 << nl
            << debugIndent << "Tet determinant = " << detA << nl
            << debugIndent << "Start local coordinates = " << y0 << nl;
    }

    // Calculate the local tracking displacement
    barycentric Tx1(x1 & T);

    DebugInfo
        << debugIndent << "Local displacement = " << Tx1 << "/" << detA << nl;

    // Calculate the hit fraction
    label iH = -1;
    scalar muH = detA > vSmall ? 1/detA : vGreat;
    for (label i = 0; i < 4; ++ i)
    {
        if (Tx1[i] < - vSmall && Tx1[i] < - mag(detA)*small)
        {
            scalar mu = - y0[i]/Tx1[i];

            DebugInfo
                << debugIndent << "Hit on tet face " << i
                << " at local coordinate " << y0 + mu*Tx1 << ", "
                << mu*detA*100 << "% of the " << "way along the track"
                << nl;

            if (0 <= mu && mu < muH)
            {
                iH = i;
                muH = mu;
            }
        }
    }

    // If there has been no hit on a degenerate or inverted tet then the
    // displacement must be within the round off error. Advance the step
    // fraction without moving and return.
    if (iH == -1 && muH == vGreat)
    {
        stepFraction += fraction;
        return Tuple2<label, scalar>(-1, 0);
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

    // Set the new position
    coordinates = yH;

    // Set the proportion of the track that has been completed
    stepFraction += fraction*muH*detA;

    if (debug)
    {
        if (iH != -1)
        {
            Info<< debugIndent << "Track hit tet face " << iH << " first" << nl;
        }
        else
        {
            Info<< debugIndent << "Track hit no tet faces" << nl;
        }

        const vector xH =
            position(mesh, coordinates, celli, facei, faceTrii, stepFraction);

        Info<< debugIndent << "End local coordinates = " << yH << nl
            << debugIndent << "End global coordinates = " << xH << nl
            << debugIndent << "Tracking displacement = " << xH - x0 << nl
            << debugIndent << muH*detA*100 << "% of the step from "
            << stepFraction - fraction*muH*detA << " to "
            << stepFraction - fraction*muH*detA + fraction
            << " completed" << nl << endl;
    }

    // Accumulate fraction behind
    if (muH*detA < small || nTracksBehind > 0)
    {
        stepFractionBehind += (fraction != 0 ? fraction : 1)*muH*detA;

        if (stepFractionBehind > rootSmall)
        {
            stepFractionBehind = 0;
            nTracksBehind = 0;
        }
        else
        {
            ++ nTracksBehind;
        }
    }

    #undef debugPrefix
    #undef debugIndent

    return Tuple2<label, scalar>(iH, iH != -1 ? 1 - muH*detA : 0);
}


Foam::Tuple2<Foam::label, Foam::scalar> Foam::tracking::toMovingTri
(
    const polyMesh& mesh,
    const vector& displacement,
    const scalar fraction,
    barycentric& coordinates,
    label& celli,
    label& facei,
    label& faceTrii,
    scalar& stepFraction,
    scalar& stepFractionBehind,
    label& nTracksBehind,
    const string& debugPrefix_
)
{
    const bool debug = notNull(debugPrefix_);
    #define debugPrefix debugPrefix_.c_str() << ": "
    #define debugIndent string(debugPrefix_.size(), ' ').c_str() << ": "

    const vector x0 =
        position(mesh, coordinates, celli, facei, faceTrii, stepFraction);
    const vector x1 = displacement;
    const barycentric y0 = coordinates;

    DebugInfo
        << debugPrefix << "Tracking from " << x0
        << " along " << x1 << " to " << x0 + x1 << nl;

    // Get the tet geometry
    Pair<vector> centre;
    FixedList<scalar, 4> detA;
    FixedList<barycentricTensor, 3> T;
    movingTetReverseTransform
    (
        mesh,
        celli,
        facei,
        faceTrii,
        stepFraction,
        fraction,
        centre,
        detA,
        T
    );

    if (debug)
    {
        Pair<vector> o, b, v1, v2;
        movingTetGeometry
        (
            mesh,
            celli,
            facei,
            faceTrii,
            stepFraction,
            fraction,
            o,
            b,
            v1,
            v2
        );

        Info<< debugIndent << "Tet points o=" << o[0] << ", b=" << b[0]
            << ", v1=" << v1[0] << ", v2=" << v2[0] << nl
            << debugIndent << "Tet determinant = " << detA[0] << nl
            << debugIndent << "Start local coordinates = " << y0[0] << nl;
    }

    // Get the relative global position
    const vector x0Rel = x0 - centre[0];
    const vector x1Rel = x1 - centre[1];

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

    if (debug)
    {
        for (label i = 0; i < 4; ++ i)
        {
            Info<< debugPrefix << (i ? "             " : "Hit equation ")
                << i << " = " << hitEqn[i] << nl;
        }
        Info<< debugPrefix << " DetA equation = " << detA << nl;
    }

    // Calculate the hit fraction
    label iH = -1;
    scalar muH = detA[0] > vSmall ? 1/detA[0] : vGreat;
    for (label i = 0; i < 4; ++ i)
    {
        const Roots<3> mu = hitEqn[i].roots();

        for (label j = 0; j < 3; ++ j)
        {
            if
            (
                mu.type(j) == rootType::real
             && hitEqn[i].derivative(mu[j]) < - vSmall
             && hitEqn[i].derivative(mu[j]) < - mag(detA[0])*small
            )
            {
                if (debug)
                {
                    const barycentric yH
                    (
                        hitEqn[0].value(mu[j]),
                        hitEqn[1].value(mu[j]),
                        hitEqn[2].value(mu[j]),
                        hitEqn[3].value(mu[j])
                    );
                    const scalar detAH = detAEqn.value(mu[j]);

                    Info<< debugPrefix << "Hit on tet face " << i
                        << " at local coordinate "
                        << (mag(detAH) > vSmall ? name(yH/detAH) : "???")
                        << ", " << mu[j]*detA[0]*100 << "% of the "
                        << "way along the track" << nl;
                }

                if (0 <= mu[j] && mu[j] < muH)
                {
                    iH = i;
                    muH = mu[j];
                }
            }
        }
    }

    // If there has been no hit on a degenerate or inverted tet then the
    // displacement must be within the round off error. Advance the step
    // fraction without moving and return.
    if (iH == -1 && muH == vGreat)
    {
        stepFraction += fraction;
        return Tuple2<label, scalar>(-1, 0);
    }

    // Set the new coordinates
    barycentric yH
    (
        hitEqn[0].value(muH),
        hitEqn[1].value(muH),
        hitEqn[2].value(muH),
        hitEqn[3].value(muH)
    );
    // !!! <-- This fails if the tet collapses onto the track, as detA tends
    // to zero at the hit. In this instance, we can differentiate the hit and
    // detA polynomials to find a limiting location, but this will not be on a
    // triangle. We will then need to do a second track through the degenerate
    // tet to find the final hit position. This second track is over zero
    // distance and therefore can be of the static mesh type. This has not yet
    // been implemented.
    const scalar detAH = detAEqn.value(muH);
    if (mag(detAH) < vSmall)
    {
        FatalErrorInFunction
            << "A moving tet collapsed onto a track. This is not supported. "
            << "The mesh is too poor, or the motion too severe, for tracking "
            << "to function." << exit(FatalError);
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
    coordinates = yH;

    // Set the proportion of the track that has been completed
    stepFraction += fraction*muH*detA[0];

    if (debug)
    {
        if (iH != -1)
        {
            Info<< debugPrefix << "Track hit tet face " << iH << " first" << nl;
        }
        else
        {
            Info<< debugPrefix << "Track hit no tet faces" << nl;
        }

        const vector xH =
            position(mesh, coordinates, celli, facei, faceTrii, stepFraction);

        Info<< debugPrefix << "End local coordinates = " << yH << nl
            << debugPrefix << "End global coordinates = " << xH << nl
            << debugPrefix << "Tracking displacement = " << xH - x0 << nl
            << debugPrefix << muH*detA[0]*100 << "% of the step from "
            << stepFraction - fraction*muH*detA[0] << " to "
            << stepFraction - fraction*muH*detA[0] + fraction
            << " completed" << nl << endl;
    }

    // Accumulate fraction behind
    if (muH*detA[0] < small || nTracksBehind > 0)
    {
        stepFractionBehind += (fraction != 0 ? fraction : 1)*muH*detA[0];

        if (stepFractionBehind > rootSmall)
        {
            stepFractionBehind = 0;
            nTracksBehind = 0;
        }
        else
        {
            ++ nTracksBehind;
        }
    }

    #undef debugPrefix
    #undef debugIndent

    return Tuple2<label, scalar>(iH, iH != -1 ? 1 - muH*detA[0] : 0);
}


Foam::Tuple2<Foam::label, Foam::scalar> Foam::tracking::toTri
(
    const polyMesh& mesh,
    const vector& displacement,
    const scalar fraction,
    barycentric& coordinates,
    label& celli,
    label& facei,
    label& faceTrii,
    scalar& stepFraction,
    scalar& stepFractionBehind,
    label& nTracksBehind,
    const string& debugPrefix
)
{
    return
        mesh.moving() && (stepFraction != 1 || fraction != 0)
      ? toMovingTri
        (
            mesh, displacement, fraction,
            coordinates, celli, facei, faceTrii, stepFraction,
            stepFractionBehind, nTracksBehind,
            debugPrefix
        )
      : toStationaryTri
        (
            mesh, displacement, fraction,
            coordinates, celli, facei, faceTrii, stepFraction,
            stepFractionBehind, nTracksBehind,
            debugPrefix
        );
}


void Foam::tracking::reflect(barycentric& coordinates)
{
    Swap(coordinates.c(), coordinates.d());
}


void Foam::tracking::rotate(const bool reverse, barycentric& coordinates)
{
    if (!reverse)
    {
        scalar temp = coordinates.b();
        coordinates.b() = coordinates.c();
        coordinates.c() = coordinates.d();
        coordinates.d() = temp;
    }
    else
    {
        scalar temp = coordinates.d();
        coordinates.d() = coordinates.c();
        coordinates.c() = coordinates.b();
        coordinates.b() = temp;
    }
}


void Foam::tracking::changeFaceTri
(
    const polyMesh& mesh,
    const label tetTrii,
    barycentric& coordinates,
    const label celli,
    label& facei,
    label& faceTrii
)
{
    const bool isOwner = mesh.faceOwner()[facei] == celli;

    const label firstTetPti = 1;
    const label lastTetPti = mesh.faces()[facei].size() - 2;

    if (tetTrii == 1)
    {
        changeFace(mesh, tetTrii, coordinates, celli, facei, faceTrii);
    }
    else if (tetTrii == 2)
    {
        if (isOwner)
        {
            if (faceTrii == lastTetPti)
            {
                changeFace(mesh, tetTrii, coordinates, celli, facei, faceTrii);
            }
            else
            {
                reflect(coordinates);
                faceTrii += 1;
            }
        }
        else
        {
            if (faceTrii == firstTetPti)
            {
                changeFace(mesh, tetTrii, coordinates, celli, facei, faceTrii);
            }
            else
            {
                reflect(coordinates);
                faceTrii -= 1;
            }
        }
    }
    else if (tetTrii == 3)
    {
        if (isOwner)
        {
            if (faceTrii == firstTetPti)
            {
                changeFace(mesh, tetTrii, coordinates, celli, facei, faceTrii);
            }
            else
            {
                reflect(coordinates);
                faceTrii -= 1;
            }
        }
        else
        {
            if (faceTrii == lastTetPti)
            {
                changeFace(mesh, tetTrii, coordinates, celli, facei, faceTrii);
            }
            else
            {
                reflect(coordinates);
                faceTrii += 1;
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


void Foam::tracking::changeFace
(
    const polyMesh& mesh,
    const label tetTrii,
    barycentric& coordinates,
    const label celli,
    label& facei,
    label& faceTrii
)
{
    // Get the old topology
    const triFace triOldIs(tetIndices(celli, facei, faceTrii).faceTriIs(mesh));

    // Get the shared edge and the pre-rotation
    edge sharedEdge;
    if (tetTrii == 1)
    {
        sharedEdge = edge(triOldIs[1], triOldIs[2]);
    }
    else if (tetTrii == 2)
    {
        sharedEdge = edge(triOldIs[2], triOldIs[0]);
    }
    else if (tetTrii == 3)
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
    faceTrii = -1;
    forAll(mesh.cells()[celli], cellFacei)
    {
        const label newFacei = mesh.cells()[celli][cellFacei];
        const class face& newFace = mesh.faces()[newFacei];
        const label newOwner = mesh.faceOwner()[newFacei];

        // Exclude the current face
        if (facei == newFacei)
        {
            continue;
        }

        // Loop over the edges, looking for the shared one
        const label edgeComp = newOwner == celli ? -1 : +1;
        label edgei = 0;
        for
        (
            ;
            edgei < newFace.size()
         && edge::compare(sharedEdge, newFace.faceEdge(edgei)) != edgeComp;
            ++ edgei
        );

        // If the face does not contain the edge, then move on to the next face
        if (edgei >= newFace.size())
        {
            continue;
        }

        // Make the edge index relative to the base point
        const label newBasei = max(0, mesh.tetBasePtIs()[newFacei]);
        edgei = (edgei - newBasei + newFace.size()) % newFace.size();

        // If the edge is next the base point (i.e., the index is 0 or n - 1),
        // then we swap it for the adjacent edge. This new edge is opposite the
        // base point, and defines the tet with the original edge in it.
        edgei = min(max(1, edgei), newFace.size() - 2);

        // Set the new face and tet point
        facei = newFacei;
        faceTrii = edgei;

        // Exit the loop now that the tet point has been found
        break;
    }

    if (faceTrii == -1)
    {
        FatalErrorInFunction
            << "The search for an edge-connected face and tet-point failed."
            << exit(FatalError);
    }

    // Pre-rotation puts the shared edge opposite the base of the tetrahedron
    if (sharedEdge.otherVertex(triOldIs[1]) == -1)
    {
        rotate(false, coordinates);
    }
    else if (sharedEdge.otherVertex(triOldIs[2]) == -1)
    {
        rotate(true, coordinates);
    }

    // Get the new topology
    const triFace triNewIs(tetIndices(celli, facei, faceTrii).faceTriIs(mesh));

    // Reflect to account for the change of triangle orientation on the new face
    reflect(coordinates);

    // Post rotation puts the shared edge back in the correct location
    if (sharedEdge.otherVertex(triNewIs[1]) == -1)
    {
        rotate(true, coordinates);
    }
    else if (sharedEdge.otherVertex(triNewIs[2]) == -1)
    {
        rotate(false, coordinates);
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::barycentric Foam::tracking::coordinates
(
    const polyMesh& mesh,
    const point& position,
    const label celli,
    const label facei,
    const label faceTrii,
    const scalar stepFraction
)
{
    static const barycentric o(1, 0, 0, 0);

    if (mesh.moving() && stepFraction != 1)
    {
        Pair<vector> centre;
        FixedList<scalar, 4> detA;
        FixedList<barycentricTensor, 3> T;
        movingTetReverseTransform
        (
            mesh,
            celli, facei, faceTrii, stepFraction, 0,
            centre, detA, T
        );

        return o + ((position - centre[0]) & T[0]/detA[0]);
    }
    else
    {
        vector centre;
        scalar detA;
        barycentricTensor T;
        stationaryTetReverseTransform
        (
            mesh,
            celli, facei, faceTrii,
            centre, detA, T
        );

        return o + ((position - centre) & T/detA);
    }
}


Foam::Pair<Foam::vector> Foam::tracking::faceNormalAndDisplacement
(
    const polyMesh& mesh,
    const barycentric& coordinates,
    const label celli,
    const label facei,
    const label faceTrii,
    const scalar stepFraction
)
{
    if (mesh.moving() && stepFraction != 1)
    {
        Pair<vector> centre, base, vertex1, vertex2;
        movingTetGeometry
        (
            mesh,
            celli,
            facei,
            faceTrii,
            stepFraction,
            1,
            centre,
            base,
            vertex1,
            vertex2
        );

        // Use the coordinates to interpolate the motion of the three face
        // vertices to the current coordinates
        return
            Pair<vector>
            (
                triPointRef(base[0], vertex1[0], vertex2[0]).normal(),
                coordinates.b()*base[1]
              + coordinates.c()*vertex1[1]
              + coordinates.d()*vertex2[1]
            );
    }
    else
    {
        vector centre, base, vertex1, vertex2;
        stationaryTetGeometry
        (
            mesh,
            celli,
            facei,
            faceTrii,
            centre,
            base,
            vertex1,
            vertex2
        );

        return
            Pair<vector>
            (
                triPointRef(base, vertex1, vertex2).normal(),
                Zero
            );
    }
}


Foam::Tuple2<bool, Foam::scalar> Foam::tracking::toFace
(
    const polyMesh& mesh,
    const vector& displacement,
    const scalar fraction,
    barycentric& coordinates,
    label& celli,
    label& facei,
    label& faceTrii,
    scalar& stepFraction,
    scalar& stepFractionBehind,
    label& nTracksBehind,
    const string& debugPrefix
)
{
    scalar f = 1;

    // Loop the tets in the current cell until the track ends or a face is hit
    while (nTracksBehind < maxNTracksBehind)
    {
        const Tuple2<label, scalar> tetTriiAndF =
            toTri
            (
                mesh, f*displacement, f*fraction,
                coordinates, celli, facei, faceTrii, stepFraction,
                stepFractionBehind, nTracksBehind,
                debugPrefix
            );

        const label tetTrii = tetTriiAndF.first();

        f *= tetTriiAndF.second();

        if (tetTrii == -1)
        {
            // The track has completed within the current tet
            return Tuple2<bool, scalar>(false, 0);
        }
        else if (tetTrii == 0)
        {
            // The track has hit a face
            return Tuple2<bool, scalar>(true, f);
        }
        else
        {
            // Move to the next tet and continue the track
            changeFaceTri
            (
                mesh, tetTrii,
                coordinates, celli, facei, faceTrii
            );
        }
    }

    // Warn if stuck, and incorrectly advance the step fraction to completion
    WarningInFunction
        << "Track got stuck at "
        << position(mesh, coordinates, celli, facei, faceTrii, stepFraction)
        << endl;

    stepFraction += f*fraction;

    stepFractionBehind = 0;
    nTracksBehind = 0;

    return Tuple2<bool, scalar>(false, 0);
}


Foam::Tuple2<bool, Foam::scalar> Foam::tracking::toCell
(
    const polyMesh& mesh,
    const vector& displacement,
    const scalar fraction,
    barycentric& coordinates,
    label& celli,
    label& facei,
    label& faceTrii,
    scalar& stepFraction,
    scalar& stepFractionBehind,
    label& nTracksBehind,
    const string& debugPrefix
)
{
    const Tuple2<bool, scalar> onFaceAndF =
        toFace
        (
            mesh, displacement, fraction,
            coordinates, celli, facei, faceTrii, stepFraction,
            stepFractionBehind, nTracksBehind,
            debugPrefix
        );

    const bool onInternalFace =
        onFaceAndF.first() && mesh.isInternalFace(facei);

    if (onInternalFace)
    {
        crossInternalFace(mesh, coordinates, celli, facei, faceTrii);
    }

    return Tuple2<bool, scalar>(onInternalFace, onFaceAndF.second());
}


Foam::Tuple2<bool, Foam::scalar> Foam::tracking::toBoundary
(
    const polyMesh& mesh,
    const vector& displacement, const scalar fraction,
    barycentric& coordinates,
    label& celli,
    label& facei,
    label& faceTrii,
    scalar& stepFraction,
    scalar& stepFractionBehind,
    label& nTracksBehind,
    const string& debugPrefix
)
{
    scalar f = 1;

    // Loop multiple cells until the track ends or a boundary face is hit
    while (true)
    {
        const Tuple2<bool, scalar> onFaceAndF =
            toFace
            (
                mesh, f*displacement, f*fraction,
                coordinates, celli, facei, faceTrii, stepFraction,
                stepFractionBehind, nTracksBehind,
                debugPrefix
            );

        f *= onFaceAndF.second();

        const bool onInternalFace =
            onFaceAndF.first() && mesh.isInternalFace(facei);

        if (onInternalFace)
        {
            crossInternalFace(mesh, coordinates, celli, facei, faceTrii);
        }
        else
        {
            const bool onBoundaryFace =
                onFaceAndF.first() && !mesh.isInternalFace(facei);

            return Tuple2<bool, scalar>(onBoundaryFace, onBoundaryFace ? f : 0);
        }
    }
}


bool Foam::tracking::locate
(
    const polyMesh& mesh,
    const point& position,
    barycentric& coordinates,
    label& celli,
    label& facei,
    label& faceTrii,
    const scalar stepFraction,
    const string& debugPrefix
)
{
    // Find the cell, if it has not been given
    if (celli < 0)
    {
        celli = mesh.cellTree().findInside(position);
    }
    if (celli < 0)
    {
        FatalErrorInFunction
            << "Cell not found for position " << position << "."
            << exit(FatalError);
    }

    // Track from the centre of the cell to the desired position
    const vector displacement = position - mesh.cellCentres()[celli];

    // Loop all cell tets to find the one containing the position. Track
    // through each tet from the cell centre. If a tet contains the position
    // then the track will end with a single toTri.
    const class cell& c = mesh.cells()[celli];
    scalar minF = vGreat;
    label minFacei = -1, minFaceTrii = -1;
    forAll(c, cellFacei)
    {
        const class face& f = mesh.faces()[c[cellFacei]];
        for (label tetPti = 1; tetPti < f.size() - 1; ++ tetPti)
        {
            coordinates = barycentric(1, 0, 0, 0);
            facei = c[cellFacei];
            faceTrii = tetPti;
            scalar stepFractionCopy = stepFraction, stepFractionBehind = 0;
            label nTracksBehind = 0;

            const Tuple2<label, scalar> tetTriiAndF =
                toTri
                (
                    mesh, displacement, 0,
                    coordinates, celli, facei, faceTrii, stepFractionCopy,
                    stepFractionBehind, nTracksBehind,
                    debugPrefix
                );

            if (tetTriiAndF.first() == -1)
            {
                return true;
            }

            if (tetTriiAndF.second() < minF)
            {
                minF = tetTriiAndF.second();
                minFacei = facei;
                minFaceTrii = faceTrii;
            }
        }
    }

    // The track must be (hopefully only slightly) outside the cell. Track into
    // the tet which got the furthest.
    coordinates = barycentric(1, 0, 0, 0);
    facei = minFacei;
    faceTrii = minFaceTrii;
    scalar stepFractionCopy = stepFraction, stepFractionBehind = 0;
    label nTracksBehind = 0;
    const Tuple2<bool, scalar> onBoundaryAndF =
        toBoundary
        (
            mesh, displacement, 0,
            coordinates, celli, facei, faceTrii, stepFractionCopy,
            stepFractionBehind, nTracksBehind,
            debugPrefix
        );

    // Return successful if in a cell
    return !onBoundaryAndF.first();
}


void Foam::tracking::crossInternalFace
(
    const polyMesh& mesh,
    barycentric& coordinates,
    label& celli,
    label& facei,
    label& faceTrii
)
{
    // Set the cell to be the one on the other side of the face
    const label ownerCelli = mesh.faceOwner()[facei];
    const bool isOwner = celli == ownerCelli;
    celli = isOwner ? mesh.faceNeighbour()[facei] : ownerCelli;

    // Reflect to account for the change of triangle orientation in the new cell
    reflect(coordinates);
}


void Foam::tracking::crossWedge
(
    const wedgePolyPatch& inPatch,
    barycentric& coordinates,
    label& celli,
    label& facei,
    label& faceTrii,
    const scalar stepFraction
)
{
    const polyMesh& mesh = inPatch.boundaryMesh().mesh();

    // Move to the other side of the mesh. Note that we track here because the
    // tet indices aren't guaranteed to be consistent between pairs of wedge
    // patches, so we'd have to search for the opposite position to jump there,
    // and its quicker and more robust to track across the cell instead. We can
    // do a simplistic track here because we know we are just moving through a
    // cell; no other patch hits or transfers need to be considered.
    scalar stepFractionCopy = stepFraction, stepFractionBehind = 0;
    label nTracksBehind = 0;
    toBoundary
    (
        mesh,
      - mesh.bounds().mag()*inPatch.centreNormal(),
        0,
        coordinates,
        celli,
        facei,
        faceTrii,
        stepFractionCopy,
        stepFractionBehind,
        nTracksBehind
    );
}


void Foam::tracking::crossCyclic
(
    const cyclicPolyPatch& inPatch,
    barycentric& coordinates,
    label& celli,
    label& facei,
    label& faceTrii
)
{
    const polyMesh& mesh = inPatch.boundaryMesh().mesh();

    // Set the topology to the neighbouring face and cell
    facei = facei - inPatch.start() + inPatch.nbrPatch().start();
    celli = mesh.faceOwner()[facei];

    // Faces either side of a conformal coupled patch are numbered in opposite
    // directions as their normals both point away from their connected
    // cells. The tet point therefore counts in the opposite direction from
    // the base point.
    faceTrii = mesh.faces()[facei].size() - 1 - faceTrii;

    // Reflect to account for the change of tri orientation in the new cell
    reflect(coordinates);
}


void Foam::tracking::inProcessor
(
    const processorPolyPatch& inPatch,
    label& celli,
    label& facei
)
{
    // Break the topology. Convert the mesh face label to a patch-local face
    // label, as this is the same on both sides of the processor interface.
    // Invalidate the cell index.
    celli = -1;
    facei -= inPatch.start();
}


void Foam::tracking::outProcessor
(
    const processorPolyPatch& outPatch,
    barycentric& coordinates,
    label& celli,
    label& facei,
    label& faceTrii
)
{
    const polyMesh& mesh = outPatch.boundaryMesh().mesh();

    // Restore the topology. Convert the patch-local face label to a mesh face
    // label. Set the cell as the face's owner.
    celli = outPatch.faceCells()[facei];
    facei += outPatch.start();

    // See tracking::crossCyclic
    faceTrii = mesh.faces()[facei].size() - 1 - faceTrii;

    // See tracking::crossCyclic
    reflect(coordinates);
}


// ************************************************************************* //
