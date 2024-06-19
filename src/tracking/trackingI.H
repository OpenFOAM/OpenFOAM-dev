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
#include "tetIndices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace tracking
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Tetrahedra

    //- Get the vertices of the current tet
    inline void stationaryTetGeometry
    (
        const polyMesh& mesh,
        const label celli,
        const label facei,
        const label faceTrii,
        vector& centre,
        vector& base,
        vector& vertex1,
        vector& vertex2
    );

    //- Get the transformation associated with the current tet. This
    //  will convert a barycentric position within the tet to a
    //  cartesian position in the global coordinate system. The
    //  conversion is x = A & y, where x is the cartesian position, y is
    //  the barycentric position and A is the transformation tensor.
    inline barycentricTensor stationaryTetTransform
    (
        const polyMesh& mesh,
        const label celli,
        const label facei,
        const label faceTrii
    );

    //- Get the vertices of the current moving tet. Two values are
    //  returned for each vertex. The first is a constant, and the
    //  second is a linear coefficient of the track fraction.
    inline void movingTetGeometry
    (
        const polyMesh& mesh,
        const label celli,
        const label facei,
        const label faceTrii,
        const scalar startStepFraction,
        const scalar endStepFraction,
        Pair<vector>& centre,
        Pair<vector>& base,
        Pair<vector>& vertex1,
        Pair<vector>& vertex2
    );

    //- Get the transformation associated with the current, moving, tet.
    //  This is of the same form as for the static case. As with the
    //  moving geometry, a linear function of the tracking fraction is
    //  returned for each component.
    inline Pair<barycentricTensor> movingTetTransform
    (
        const polyMesh& mesh,
        const label celli,
        const label facei,
        const label faceTrii,
        const scalar startStepFraction,
        const scalar endStepFraction
    );


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace tracking
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

inline void Foam::tracking::stationaryTetGeometry
(
    const polyMesh& mesh,
    const label celli,
    const label facei,
    const label faceTrii,
    vector& centre,
    vector& base,
    vector& vertex1,
    vector& vertex2
)
{
    const triFace triIs(tetIndices(celli, facei, faceTrii).faceTriIs(mesh));
    const vectorField& ccs = mesh.cellCentres();
    const pointField& pts = mesh.points();

    centre = ccs[celli];
    base = pts[triIs[0]];
    vertex1 = pts[triIs[1]];
    vertex2 = pts[triIs[2]];
}


inline Foam::barycentricTensor Foam::tracking::stationaryTetTransform
(
    const polyMesh& mesh,
    const label celli,
    const label facei,
    const label faceTrii
)
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

    return barycentricTensor(centre, base, vertex1, vertex2);
}


inline void Foam::tracking::movingTetGeometry
(
    const polyMesh& mesh,
    const label celli,
    const label facei,
    const label faceTrii,
    const scalar startStepFraction,
    const scalar endStepFraction,
    Pair<vector>& centre,
    Pair<vector>& base,
    Pair<vector>& vertex1,
    Pair<vector>& vertex2
)
{
    const triFace triIs(tetIndices(celli, facei, faceTrii).faceTriIs(mesh));
    const pointField& ptsOld = mesh.oldPoints();
    const pointField& ptsNew = mesh.points();
    const vector ccOld = mesh.oldCellCentres()[celli];
    const vector ccNew = mesh.cellCentres()[celli];

    const scalar f0 = startStepFraction, f1 = endStepFraction;

    centre[0] = ccOld + f0*(ccNew - ccOld);
    base[0] = ptsOld[triIs[0]] + f0*(ptsNew[triIs[0]] - ptsOld[triIs[0]]);
    vertex1[0] = ptsOld[triIs[1]] + f0*(ptsNew[triIs[1]] - ptsOld[triIs[1]]);
    vertex2[0] = ptsOld[triIs[2]] + f0*(ptsNew[triIs[2]] - ptsOld[triIs[2]]);

    centre[1] = f1*(ccNew - ccOld);
    base[1] = f1*(ptsNew[triIs[0]] - ptsOld[triIs[0]]);
    vertex1[1] = f1*(ptsNew[triIs[1]] - ptsOld[triIs[1]]);
    vertex2[1] = f1*(ptsNew[triIs[2]] - ptsOld[triIs[2]]);
}


inline Foam::Pair<Foam::barycentricTensor> Foam::tracking::movingTetTransform
(
    const polyMesh& mesh,
    const label celli,
    const label facei,
    const label faceTrii,
    const scalar startStepFraction,
    const scalar endStepFraction
)
{
    Pair<vector> centre, base, vertex1, vertex2;
    movingTetGeometry
    (
        mesh,
        celli,
        facei,
        faceTrii,
        startStepFraction,
        endStepFraction,
        centre,
        base,
        vertex1,
        vertex2
    );

    return
        Pair<barycentricTensor>
        (
            barycentricTensor(centre[0], base[0], vertex1[0], vertex2[0]),
            barycentricTensor(centre[1], base[1], vertex1[1], vertex2[1])
        );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

inline Foam::point Foam::tracking::position
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
        return
            movingTetTransform(mesh, celli, facei, faceTrii, stepFraction, 1)[0]
          & coordinates;
    }
    else
    {
        return
            stationaryTetTransform(mesh, celli, facei, faceTrii)
          & coordinates;
    }
}


// ************************************************************************* //
