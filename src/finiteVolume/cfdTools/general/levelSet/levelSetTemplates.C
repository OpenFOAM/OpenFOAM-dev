/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenFOAM Foundation
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

#include "levelSet.H"
#include "cut.H"
#include "polyMeshTetDecomposition.H"
#include "tetIndices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::DimensionedField<Type, Foam::volMesh>> Foam::levelSetAverage
(
    const fvMesh& mesh,
    const scalarField& levelC,
    const scalarField& levelP,
    const DimensionedField<Type, volMesh>& positiveC,
    const DimensionedField<Type, pointMesh>& positiveP,
    const DimensionedField<Type, volMesh>& negativeC,
    const DimensionedField<Type, pointMesh>& negativeP
)
{
    DimensionedField<Type, volMesh> sum
    (
        IOobject
        (
            "levelSetIntegral",
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensioned<Type>
        (
            "0",
            positiveC.dimensions()*dimVolume,
            pTraits<Type>::zero
        )
    );

    forAll(sum, cI)
    {
        const List<tetIndices> cellTetIs =
            polyMeshTetDecomposition::cellTetIndices(mesh, cI);

        forAll(cellTetIs, cellTetI)
        {
            const tetIndices& tetIs = cellTetIs[cellTetI];
            const face& f = mesh.faces()[tetIs.face()];

            const label pI0 = f[tetIs.faceBasePt()];
            const label pIA = f[tetIs.facePtA()];
            const label pIB = f[tetIs.facePtB()];

            const FixedList<point, 4>
                tet =
                {
                    mesh.cellCentres()[cI],
                    mesh.points()[pI0],
                    mesh.points()[pIA],
                    mesh.points()[pIB]
                };
            const FixedList<scalar, 4>
                level =
                {
                    levelC[cI],
                    levelP[pI0],
                    levelP[pIA],
                    levelP[pIA]
                };
            const cut::volumeIntegrateOp<Type>
                positive = FixedList<Type, 4>
                ({
                    positiveC[cI],
                    positiveP[pI0],
                    positiveP[pIA],
                    positiveP[pIB]
                });
            const cut::volumeIntegrateOp<Type>
                negative = FixedList<Type, 4>
                ({
                    negativeC[cI],
                    negativeP[pI0],
                    negativeP[pIA],
                    negativeP[pIB]
                });

            sum[cI] += tetCut(tet, level, positive, negative);
        }
    }

    return sum/mesh.V();
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::levelSetAverage
(
    const fvPatch& patch,
    const scalarField& levelF,
    const scalarField& levelP,
    const Field<Type>& positiveF,
    const Field<Type>& positiveP,
    const Field<Type>& negativeF,
    const Field<Type>& negativeP
)
{
    typedef typename outerProduct<Type, vector>::type sumType;

    Field<sumType> sum(patch.size(), pTraits<sumType>::zero);

    forAll(sum, fI)
    {
        const face& f = patch.patch().localFaces()[fI];

        for(label eI = 0; eI < f.size(); ++ eI)
        {
            const edge e = f.faceEdge(eI);

            const FixedList<point, 3>
                tri =
                {
                    patch.patch().faceCentres()[fI],
                    patch.patch().localPoints()[e[0]],
                    patch.patch().localPoints()[e[1]]
                };
            const FixedList<scalar, 3>
                level =
                {
                    levelF[fI],
                    levelP[e[0]],
                    levelP[e[1]]
                };
            const cut::areaIntegrateOp<Type>
                positive = FixedList<Type, 3>
                ({
                    positiveF[fI],
                    positiveP[e[0]],
                    positiveP[e[1]]
                });
            const cut::areaIntegrateOp<Type>
                negative = FixedList<Type, 3>
                ({
                    negativeF[fI],
                    negativeP[e[0]],
                    negativeP[e[1]]
                });

            sum[fI] += triCut(tri, level, positive, negative);
        }
    }

    return sum & patch.Sf()/sqr(patch.magSf());
}


// ************************************************************************* //
