/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "polyMesh.H"
#include "Time.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

inline void Foam::particle::findTris
(
    const vector& position,
    DynamicList<label>& faceList,
    const tetPointRef& tet,
    const FixedList<vector, 4>& tetAreas,
    const FixedList<label, 4>& tetPlaneBasePtIs,
    const scalar tol
) const
{
    faceList.clear();

    const point Ct = tet.centre();

    for (label i = 0; i < 4; i++)
    {
        scalar lambda = tetLambda
        (
            Ct,
            position,
            i,
            tetAreas[i],
            tetPlaneBasePtIs[i],
            cellI_,
            tetFaceI_,
            tetPtI_,
            tol
        );

        if ((lambda > 0.0) && (lambda < 1.0))
        {
            faceList.append(i);
        }
    }
}


inline Foam::scalar Foam::particle::tetLambda
(
    const vector& from,
    const vector& to,
    const label triI,
    const vector& n,
    const label tetPlaneBasePtI,
    const label cellI,
    const label tetFaceI,
    const label tetPtI,
    const scalar tol
) const
{
    const pointField& pPts = mesh_.points();

    if (mesh_.moving())
    {
        return movingTetLambda
        (
            from,
            to,
            triI,
            n,
            tetPlaneBasePtI,
            cellI,
            tetFaceI,
            tetPtI,
            tol
        );
    }

    const point& base = pPts[tetPlaneBasePtI];

    scalar lambdaNumerator = (base - from) & n;
    scalar lambdaDenominator = (to - from) & n;

    // n carries the area of the tet faces, so the dot product with a
    // delta-length has the units of volume.  Comparing the component of each
    // delta-length in the direction of n times the face area to a fraction of
    // the cell volume.

    if (mag(lambdaDenominator) < tol)
    {
        if (mag(lambdaNumerator) < tol)
        {
            // Track starts on the face, and is potentially
            // parallel to it.  +-tol/+-tol is not a good
            // comparison, return 0.0, in anticipation of tet
            // centre correction.

            return 0.0;
        }
        else
        {
            if (mag((to - from)) < tol/mag(n))
            {
                // 'Zero' length track (compared to the tolerance, which is
                // based on the cell volume, divided by the tet face area), not
                // along the face, face cannot be crossed.

                return GREAT;
            }
            else
            {
                // Trajectory is non-zero and parallel to face

                lambdaDenominator = sign(lambdaDenominator)*SMALL;
            }
        }
    }

    return lambdaNumerator/lambdaDenominator;
}


inline Foam::scalar Foam::particle::movingTetLambda
(
    const vector& from,
    const vector& to,
    const label triI,
    const vector& n,
    const label tetPlaneBasePtI,
    const label cellI,
    const label tetFaceI,
    const label tetPtI,
    const scalar tol
) const
{
    const pointField& pPts = mesh_.points();
    const pointField& oldPPts = mesh_.oldPoints();

    // Base point of plane at end of motion
    const point& b = pPts[tetPlaneBasePtI];

    // n: Normal of plane at end of motion

    // Base point of plane at start of timestep
    const point& b00 = oldPPts[tetPlaneBasePtI];

    // Base point of plane at start of tracking portion (cast forward by
    // stepFraction)
    point b0 = b00 + stepFraction_*(b - b00);

    // Normal of plane at start of tracking portion
    vector n0 = vector::zero;

    {
        tetIndices tetIs(cellI, tetFaceI, tetPtI, mesh_);

        // tet at timestep start
        tetPointRef tet00 = tetIs.oldTet(mesh_);

        // tet at timestep end
        tetPointRef tet = tetIs.tet(mesh_);

        point tet0PtA = tet00.a() + stepFraction_*(tet.a() - tet00.a());
        point tet0PtB = tet00.b() + stepFraction_*(tet.b() - tet00.b());
        point tet0PtC = tet00.c() + stepFraction_*(tet.c() - tet00.c());
        point tet0PtD = tet00.d() + stepFraction_*(tet.d() - tet00.d());

        // Tracking portion start tet (cast forward by stepFraction)
        tetPointRef tet0(tet0PtA, tet0PtB, tet0PtC, tet0PtD);

        switch (triI)
        {
            case 0:
            {
                n0 = tet0.Sa();
                break;
            }
            case 1:
            {
                n0 = tet0.Sb();
                break;
            }
            case 2:
            {
                n0 = tet0.Sc();
                break;
            }
            case 3:
            {
                n0 = tet0.Sd();
                break;
            }
            default:
            {
                break;
            }
        }
    }

    if (mag(n0) < SMALL)
    {
        // If the old normal is zero (for example in layer addition)
        // then use the current normal;

        n0 = n;
    }

    scalar lambdaNumerator = 0;
    scalar lambdaDenominator = 0;

    vector dP = to - from;
    vector dN = n - n0;
    vector dB = b - b0;
    vector dS = from - b0;

    if (mag(dN) > SMALL)
    {
        scalar a = (dP - dB) & dN;
        scalar b = ((dP - dB) & n0) + (dS & dN);
        scalar c = dS & n0;

        if (mag(a) > SMALL)
        {

            // Solve quadratic for lambda
            scalar discriminant = sqr(b) - 4.0*a*c;

            if (discriminant < 0)
            {
                // Imaginary roots only - face not crossed
                return GREAT;
            }
            else
            {
                // Numerical Recipes in C,
                // Second Edition (1992),
                // Section 5.6.
                // q = -0.5*(b + sgn(b)*sqrt(b^2 - 4ac))
                // x1 = q/a
                // x2 = c/q

                scalar q = -0.5*(b + sign(b)*Foam::sqrt(discriminant));

                if (mag(q) < VSMALL)
                {
                    // If q is zero, then l1 = q/a is the required
                    // value of lambda, and is zero.

                    return 0.0;
                }

                scalar l1 = q/a;
                scalar l2 = c/q;

                // There will be two roots, a big one and a little
                // one, choose the little one.

                if (mag(l1) < mag(l2))
                {
                    return l1;
                }
                else
                {
                    return l2;
                }
            }
        }
        {
            // When a is zero, solve the first order polynomial

            lambdaNumerator = -c;
            lambdaDenominator = b;
        }
    }
    else
    {
        // when n = n0 is zero, there is no plane rotation, solve the
        // first order polynomial

        lambdaNumerator = -(dS & n0);
        lambdaDenominator = ((dP - dB) & n0);

    }

    if (mag(lambdaDenominator) < tol)
    {
        if (mag(lambdaNumerator) < tol)
        {
            // Track starts on the face, and is potentially
            // parallel to it.  +-tol)/+-tol is not a good
            // comparison, return 0.0, in anticipation of tet
            // centre correction.

            return 0.0;
        }
        else
        {
            if (mag((to - from)) < tol/mag(n))
            {
                // Zero length track, not along the face, face
                // cannot be crossed.

                return GREAT;
            }
            else
            {
                // Trajectory is non-zero and parallel to face

                lambdaDenominator = sign(lambdaDenominator)*SMALL;
            }
        }
    }

    return lambdaNumerator/lambdaDenominator;
}



inline void Foam::particle::tetNeighbour(label triI)
{
    const labelList& pOwner = mesh_.faceOwner();
    const faceList& pFaces = mesh_.faces();

    bool own = (pOwner[tetFaceI_] == cellI_);

    const Foam::face& f = pFaces[tetFaceI_];

    label tetBasePtI = mesh_.tetBasePtIs()[tetFaceI_];

    if (tetBasePtI == -1)
    {
        FatalErrorIn
        (
            "inline void Foam::particle::tetNeighbour(label triI)"
        )
            << "No base point for face " << tetFaceI_ << ", " << f
            << ", produces a valid tet decomposition."
            << abort(FatalError);
    }

    label facePtI = (tetPtI_ + tetBasePtI) % f.size();
    label otherFacePtI = f.fcIndex(facePtI);

    switch (triI)
    {
        case 0:
        {
            // Crossing this triangle changes tet to that in the
            // neighbour cell over tetFaceI

            // Modification of cellI_ will happen by other indexing,
            // tetFaceI_ and tetPtI don't change.

            break;
        }
        case 1:
        {
            crossEdgeConnectedFace
            (
                cellI_,
                tetFaceI_,
                tetPtI_,
                Foam::edge(f[facePtI], f[otherFacePtI])
            );

            break;
        }
        case 2:
        {
            if (own)
            {
                if (tetPtI_ < f.size() - 2)
                {
                    tetPtI_ = f.fcIndex(tetPtI_);
                }
                else
                {
                    crossEdgeConnectedFace
                    (
                        cellI_,
                        tetFaceI_,
                        tetPtI_,
                        Foam::edge(f[tetBasePtI], f[otherFacePtI])
                    );
                }
            }
            else
            {
                if (tetPtI_ > 1)
                {
                    tetPtI_ = f.rcIndex(tetPtI_);
                }
                else
                {
                    crossEdgeConnectedFace
                    (
                        cellI_,
                        tetFaceI_,
                        tetPtI_,
                        Foam::edge(f[tetBasePtI], f[facePtI])
                    );
                }
            }

            break;
        }
        case 3:
        {
            if (own)
            {
                if (tetPtI_ > 1)
                {
                    tetPtI_ = f.rcIndex(tetPtI_);
                }
                else
                {
                    crossEdgeConnectedFace
                    (
                        cellI_,
                        tetFaceI_,
                        tetPtI_,
                        Foam::edge(f[tetBasePtI], f[facePtI])
                    );
                }
            }
            else
            {
                if (tetPtI_ < f.size() - 2)
                {
                    tetPtI_ = f.fcIndex(tetPtI_);
                }
                else
                {
                    crossEdgeConnectedFace
                    (
                        cellI_,
                        tetFaceI_,
                        tetPtI_,
                        Foam::edge(f[tetBasePtI], f[otherFacePtI])
                    );
                }
            }

            break;
        }
        default:
        {
            FatalErrorIn
            (
                "inline void "
                "Foam::particle::tetNeighbour(label triI)"
            )
                << "Tet tri face index error, can only be 0..3, supplied "
                << triI << abort(FatalError);

            break;
        }
    }
}


inline void Foam::particle::crossEdgeConnectedFace
(
    const label& cellI,
    label& tetFaceI,
    label& tetPtI,
    const edge& e
)
{
    const faceList& pFaces = mesh_.faces();
    const cellList& pCells = mesh_.cells();

    const Foam::face& f = pFaces[tetFaceI];

    const Foam::cell& thisCell = pCells[cellI];

    forAll(thisCell, cFI)
    {
        // Loop over all other faces of this cell and
        // find the one that shares this edge

        label fI = thisCell[cFI];

        if (tetFaceI == fI)
        {
            continue;
        }

        const Foam::face& otherFace = pFaces[fI];

        label edDir = otherFace.edgeDirection(e);

        if (edDir == 0)
        {
            continue;
        }
        else if (f == pFaces[fI])
        {
            // This is a necessary condition if using duplicate baffles
            // (so coincident faces). We need to make sure we don't cross into
            // the face with the same vertices since we might enter a tracking
            // loop where it never exits. This test should be cheap
            // for most meshes so can be left in for 'normal' meshes.

            continue;
        }
        else
        {
            //Found edge on other face
            tetFaceI = fI;

            label eIndex = -1;

            if (edDir == 1)
            {
                // Edge is in the forward circulation of this face, so
                // work with the start point of the edge

                eIndex = findIndex(otherFace, e.start());
            }
            else
            {
                // edDir == -1, so the edge is in the reverse
                // circulation of this face, so work with the end
                // point of the edge

                eIndex = findIndex(otherFace, e.end());
            }

            label tetBasePtI = mesh_.tetBasePtIs()[fI];

            if (tetBasePtI == -1)
            {
                FatalErrorIn
                (
                    "inline void "
                    "Foam::particle::crossEdgeConnectedFace"
                    "("
                        "const label& cellI,"
                        "label& tetFaceI,"
                        "label& tetPtI,"
                        "const edge& e"
                    ")"
                )
                    << "No base point for face " << fI << ", " << f
                    << ", produces a decomposition that has a minimum "
                    << "volume greater than tolerance."
                    << abort(FatalError);
            }

            // Find eIndex relative to the base point on new face
            eIndex -= tetBasePtI;

            if (neg(eIndex))
            {
                eIndex = (eIndex + otherFace.size()) % otherFace.size();
            }

            if (eIndex == 0)
            {
                // The point is the base point, so this is first tet
                // in the face circulation

                tetPtI = 1;
            }
            else if (eIndex == otherFace.size() - 1)
            {
                // The point is the last before the base point, so
                // this is the last tet in the face circulation

                tetPtI = otherFace.size() - 2;
            }
            else
            {
                tetPtI = eIndex;
            }

            break;
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::label Foam::particle::getNewParticleID() const
{
    label id = particleCount_++;

    if (id == labelMax)
    {
        WarningIn("particle::getNewParticleID() const")
            << "Particle counter has overflowed. This might cause problems"
            << " when reconstructing particle tracks." << endl;
    }
    return id;
}


inline const Foam::polyMesh& Foam::particle::mesh() const
{
    return mesh_;
}


inline const Foam::vector& Foam::particle::position() const
{
    return position_;
}


inline Foam::vector& Foam::particle::position()
{
    return position_;
}


inline Foam::label Foam::particle::cell() const
{
    return cellI_;
}


inline Foam::label& Foam::particle::cell()
{
    return cellI_;
}


inline Foam::label Foam::particle::tetFace() const
{
    return tetFaceI_;
}


inline Foam::label& Foam::particle::tetFace()
{
    return tetFaceI_;
}


inline Foam::label Foam::particle::tetPt() const
{
    return tetPtI_;
}


inline Foam::label& Foam::particle::tetPt()
{
    return tetPtI_;
}


inline Foam::tetIndices Foam::particle::currentTetIndices() const
{
    return tetIndices(cellI_, tetFaceI_, tetPtI_, mesh_);
}


inline Foam::tetPointRef Foam::particle::currentTet() const
{
    return currentTetIndices().tet(mesh_);
}


inline Foam::vector Foam::particle::normal() const
{
    return currentTetIndices().faceTri(mesh_).normal();
}


inline Foam::vector Foam::particle::oldNormal() const
{
    return currentTetIndices().oldFaceTri(mesh_).normal();
}


inline Foam::label Foam::particle::face() const
{
    return faceI_;
}


inline Foam::label& Foam::particle::face()
{
    return faceI_;
}


inline void Foam::particle::initCellFacePt()
{
    if (cellI_ == -1)
    {
        mesh_.findCellFacePt
        (
            position_,
            cellI_,
            tetFaceI_,
            tetPtI_
        );

        if (cellI_ == -1)
        {
            FatalErrorIn("void Foam::particle::initCellFacePt()")
                << "cell, tetFace and tetPt search failure at position "
                << position_ << abort(FatalError);
        }
    }
    else
    {
        mesh_.findTetFacePt(cellI_, position_, tetFaceI_, tetPtI_);

        if (tetFaceI_ == -1 || tetPtI_ == -1)
        {
            label oldCellI = cellI_;

            mesh_.findCellFacePt
            (
                position_,
                cellI_,
                tetFaceI_,
                tetPtI_
            );

            if (cellI_ == -1 || tetFaceI_ == -1 || tetPtI_ == -1)
            {
                // The particle has entered this function with a cell
                // number, but hasn't been able to find a cell to
                // occupy.

                if (!mesh_.pointInCellBB(position_, oldCellI, 0.1))
                {
                    // If the position is not inside the (slightly
                    // extended) bound-box of the cell that it thought
                    // it should be in, then this is considered an
                    // error.

                    FatalErrorIn("void Foam::particle::initCellFacePt()")
                        << "cell, tetFace and tetPt search failure at position "
                        << position_ << nl << "for requested cell " << oldCellI
                        << abort(FatalError);
                }

                // The position is in the (slightly extended)
                // bound-box of the cell.  This situation may arise
                // because the face decomposition of the cell is not
                // the same as when the particle acquired the cell
                // index.  For example, it has been read into a mesh
                // that has made a different face base-point decision
                // for a boundary face and now this particle is in a
                // position that is not in the mesh.  Gradually move
                // the particle towards the centre of the cell that it
                // thought that it was in.

                cellI_ = oldCellI;

                point newPosition = position_;

                const point& cC = mesh_.cellCentres()[cellI_];

                label trap(1.0/trackingCorrectionTol + 1);

                label iterNo = 0;

                do
                {
                    newPosition += trackingCorrectionTol*(cC - position_);

                    mesh_.findTetFacePt
                    (
                        cellI_,
                        newPosition,
                        tetFaceI_,
                        tetPtI_
                    );

                    iterNo++;

                } while (tetFaceI_ < 0  && iterNo <= trap);

                if (tetFaceI_ == -1)
                {
                    FatalErrorIn("void Foam::particle::initCellFacePt()")
                        << "cell, tetFace and tetPt search failure at position "
                        << position_ << abort(FatalError);
                }

                if (debug)
                {
                    WarningIn("void Foam::particle::initCellFacePt()")
                        << "Particle moved from " << position_
                        << " to " << newPosition
                        << " in cell " << cellI_
                        << " tetFace " << tetFaceI_
                        << " tetPt " << tetPtI_ << nl
                        << "    (A fraction of "
                        << 1.0 - mag(cC - newPosition)/mag(cC - position_)
                        << " of the distance to the cell centre)"
                        << " because a decomposition tetFace and tetPt "
                        << "could not be found."
                        << endl;
                }

                position_ = newPosition;
            }

            if (debug && cellI_ != oldCellI)
            {
                WarningIn("void Foam::particle::initCellFacePt()")
                    << "Particle at position " << position_
                    << " searched for a cell, tetFace and tetPt." << nl
                    << "    Found"
                    << " cell " << cellI_
                    << " tetFace " << tetFaceI_
                    << " tetPt " << tetPtI_ << nl
                    << "    This is a different cell to that which was supplied"
                    << " (" << oldCellI << ")." << nl
                    << endl;
            }
        }
    }
}


inline bool Foam::particle::onBoundary() const
{
    return faceI_ != -1 && faceI_ >= mesh_.nInternalFaces();
}


inline Foam::scalar& Foam::particle::stepFraction()
{
    return stepFraction_;
}


inline Foam::scalar Foam::particle::stepFraction() const
{
    return stepFraction_;
}


inline Foam::label Foam::particle::origProc() const
{
    return origProc_;
}


inline Foam::label& Foam::particle::origProc()
{
    return origProc_;
}


inline Foam::label Foam::particle::origId() const
{
    return origId_;
}


inline Foam::label& Foam::particle::origId()
{
    return origId_;
}


inline bool Foam::particle::softImpact() const
{
    return false;
}


inline Foam::scalar Foam::particle::currentTime() const
{
    return
        mesh_.time().value()
      + stepFraction_*mesh_.time().deltaTValue();
}


inline bool Foam::particle::internalFace(const label faceI) const
{
    return mesh_.isInternalFace(faceI);
}


bool Foam::particle::boundaryFace(const label faceI) const
{
    return !internalFace(faceI);
}


inline Foam::label Foam::particle::patch(const label faceI) const
{
    return mesh_.boundaryMesh().whichPatch(faceI);
}


inline Foam::label Foam::particle::patchFace
(
    const label patchI,
    const label faceI
) const
{
    return mesh_.boundaryMesh()[patchI].whichFace(faceI);
}


inline Foam::label Foam::particle::faceInterpolation() const
{
    return faceI_;
}


// ************************************************************************* //
