/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2025 OpenFOAM Foundation
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

#include "cellPoint_interpolation.H"
#include "polyMeshTetDecomposition.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class Type>
inline void Foam::interpolations::cellPointBase<Type>::tetTransform
(
    const tetIndices& tetIs,
    vector& x0,
    scalar& detA,
    barycentricTensor& T
) const
{
    const triFace triIs(tetIs.faceTriIs(this->mesh_));

    const vectorField& ccs = this->mesh_.cellCentres();
    const pointField& ps = this->mesh_.points();

    const barycentricTensor A
    (
        ccs[tetIs.cell()],
        ps[triIs[0]],
        ps[triIs[1]],
        ps[triIs[2]]
    );

    const vector ba = A.a() - A.b();
    const vector bc = A.c() - A.b();
    const vector bd = A.d() - A.b();
    const vector ac = A.c() - A.a();
    const vector ad = A.d() - A.a();

    x0 = A.b();

    detA = ba & (bc ^ bd);

    T = barycentricTensor(bc ^ bd, ad ^ ac, bd ^ ba, ba ^ bc);
}


template<class Type>
inline void Foam::interpolations::cellPointBase<Type>::tetFaceTriTransform
(
    const tetIndices& tetIs,
    vector& x0,
    scalar& detA,
    barycentricTensor& T
) const
{
    const triFace triIs(tetIs.faceTriIs(this->mesh_));

    const vectorField& ccs = this->mesh_.cellCentres();
    const pointField& ps = this->mesh_.points();

    const barycentricTensor A
    (
        ccs[tetIs.cell()],
        ps[triIs[0]],
        ps[triIs[1]],
        ps[triIs[2]]
    );

    const vector bc = A.c() - A.b();
    const vector bd = A.d() - A.b();
    const vector dc = A.c() - A.d();

    x0 = A.b();

    const vector bc_bd = bc ^ bd;

    detA = mag(bc_bd);

    const vector n = detA != 0 ? bc_bd/detA : vector::zero;

    T = barycentricTensor(vector::zero, dc ^ n, bd ^ n, n ^ bc);
}


template<class Type>
inline Foam::barycentric Foam::interpolations::cellPointBase<Type>::y0() const
{
    return barycentric(0, 1, 0, 0);
}


template<class Type>
inline Foam::tetIndices Foam::interpolations::cellPointBase<Type>::findTet
(
    const point& x,
    const label celli,
    const label facei,
    vector& x0,
    scalar& detA,
    barycentricTensor& T
) const
{
    const List<tetIndices> tetIss =
        facei == -1
      ? polyMeshTetDecomposition::cellTetIndices
        (
            this->mesh_,
            celli
        )
      : polyMeshTetDecomposition::faceTetIndices
        (
            this->mesh_,
            facei,
            this->mesh_.faceOwner()[facei]
        );

    // Find the containing tetrahedron
    forAll(tetIss, teti)
    {
        if (facei == -1)
        {
            tetTransform(tetIss[teti], x0, detA, T);
        }
        else
        {
            tetFaceTriTransform(tetIss[teti], x0, detA, T);
        }

        if (detA <= 0) continue;

        const barycentric detAY = detA*y0() + ((x - x0) & T);

        if (cmptMin(detAY) + detA*small > 0)
        {
            return tetIss[teti];
        }
    }

    // A containing tet/tri was not found, so look instead for the nearest
    scalar minDist = vGreat;
    label minTeti = -1;
    forAll(tetIss, teti)
    {
        const scalar dist =
            facei == -1
          ? tetIss[teti].tet(this->mesh_).nearestPoint(x).distance()
          : tetIss[teti].faceTri(this->mesh_).nearestPoint(x).distance();

        if (dist < minDist)
        {
            minDist = dist;
            minTeti = teti;
        }
    }

    // Get the transform for the nearest tetrahedron
    if (facei == -1)
    {
        tetTransform(tetIss[minTeti], x0, detA, T);
    }
    else
    {
        tetFaceTriTransform(tetIss[minTeti], x0, detA, T);
    }

    return tetIss[minTeti];
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
inline Type Foam::interpolations::cellPointBase<Type>::interpolate
(
    const vector& position,
    const label celli,
    const label facei
) const
{
    vector x0;
    scalar detA;
    barycentricTensor T;

    // Find the tetrahedron
    const tetIndices tetIs = this->findTet(position, celli, facei, x0, detA, T);

    // Determine the barycentric coordinates. If the tetrahedron is degenerate
    // then just take the centre.
    const barycentric y =
        detA == 0
      ? barycentric::uniform(0.25)
      : y0() + ((position - x0) & T)/detA;

    // Defer to the barycentric interpolation overload. Note that the face
    // index is not needed.
    return cellPointBase<Type>::interpolate(y, tetIs, -1);
}


template<class Type>
inline Type Foam::interpolations::cellPointBase<Type>::interpolate
(
    const barycentric& coordinates,
    const tetIndices& tetIs,
    const label facei
) const
{
    // Ignore the face index. If the location is on the face then the first
    // component of the barycentric coordinates should be zero, meaning we get
    // interpolation in the face anyway.

    const triFace triIs = tetIs.faceTriIs(this->mesh_);

    return
        this->psi_[tetIs.cell()]*coordinates[0]
      + this->psip_[triIs[0]]*coordinates[1]
      + this->psip_[triIs[1]]*coordinates[2]
      + this->psip_[triIs[2]]*coordinates[3];
}


template<class Type>
inline Foam::interpolationGradType<Type>
Foam::interpolations::cellPointGradBase<Type>::interpolateGrad
(
    const tetIndices& tetIs,
    const scalar detA,
    const barycentricTensor T
) const
{
    const triFace triIs = tetIs.faceTriIs(this->mesh_);

    const Type& psia = this->psi_[tetIs.cell()];
    const Type& psib = this->psip_[triIs[0]];
    const Type& psic = this->psip_[triIs[1]];
    const Type& psid = this->psip_[triIs[2]];

    return
        interpolationGradType<Type>
        (
            T.x().a()*psia + T.x().b()*psib + T.x().c()*psic + T.x().d()*psid,
            T.y().a()*psia + T.y().b()*psib + T.y().c()*psic + T.y().d()*psid,
            T.z().a()*psia + T.z().b()*psib + T.z().c()*psic + T.z().d()*psid
        )/detA;
}


template<class Type>
inline Foam::interpolationGradType<Type>
Foam::interpolations::cellPointGradBase<Type>::interpolateGrad
(
    const vector& position,
    const label celli,
    const label facei
) const
{
    vector x0;
    scalar detA;
    barycentricTensor T;

    const tetIndices tetIs = this->findTet(position, celli, facei, x0, detA, T);

    return interpolateGrad(tetIs, detA, T);
}


template<class Type>
inline Foam::interpolationGradType<Type>
Foam::interpolations::cellPointGradBase<Type>::interpolateGrad
(
    const barycentric&,
    const tetIndices& tetIs,
    const label
) const
{
    vector x0;
    scalar detA;
    barycentricTensor T;

    this->tetTransform(tetIs, x0, detA, T);

    return interpolateGrad(tetIs, detA, T);
}


// ************************************************************************* //
