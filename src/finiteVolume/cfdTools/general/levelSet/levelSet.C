/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2022 OpenFOAM Foundation
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
#include "cutTriTet.H"
#include "polyMeshTetDecomposition.H"
#include "tetIndices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField> Foam::levelSetFraction
(
    const fvMesh& mesh,
    const scalarField& levelC,
    const scalarField& levelP,
    const bool above
)
{
    typedef cutTriTet::noOp noOp;
    typedef cutTriTet::volumeOp volumeOp;

    tmp<scalarField> tResult(new scalarField(mesh.nCells(), Zero));
    scalarField& result = tResult.ref();

    forAll(result, cI)
    {
        const List<tetIndices> cellTetIs =
            polyMeshTetDecomposition::cellTetIndices(mesh, cI);

        scalar v = 0, r = 0;

        forAll(cellTetIs, cellTetI)
        {
            const triFace triIs = cellTetIs[cellTetI].faceTriIs(mesh);

            const FixedList<point, 4>
                tet =
                {
                    mesh.cellCentres()[cI],
                    mesh.points()[triIs[0]],
                    mesh.points()[triIs[1]],
                    mesh.points()[triIs[2]]
                };
            const FixedList<scalar, 4>
                level =
                {
                    levelC[cI],
                    levelP[triIs[0]],
                    levelP[triIs[1]],
                    levelP[triIs[2]]
                };

            v += volumeOp()(tet);

            if (above)
            {
                r += tetCut(tet, level, volumeOp(), noOp());
            }
            else
            {
                r += tetCut(tet, level, noOp(), volumeOp());
            }
        }

        result[cI] = r/v;
    }

    return tResult;
}


Foam::tmp<Foam::scalarField> Foam::levelSetFraction
(
    const fvPatch& patch,
    const scalarField& levelF,
    const scalarField& levelP,
    const bool above
)
{
    typedef cutTriTet::noOp noOp;
    typedef cutTriTet::areaMagOp areaMagOp;

    tmp<scalarField> tResult(new scalarField(patch.size(), 0));
    scalarField& result = tResult.ref();

    forAll(result, fI)
    {
        const face& f = patch.patch().localFaces()[fI];

        scalar a = 0, r = 0;

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

            a += areaMagOp()(tri);

            if (above)
            {
                r += triCut(tri, level, areaMagOp(), noOp());
            }
            else
            {
                r += triCut(tri, level, noOp(), areaMagOp());
            }
        }

        result[fI] = r/a;
    }

    return tResult;
}


Foam::tmp<Foam::volScalarField> Foam::levelSetFraction
(
    const volScalarField& levelC,
    const pointScalarField& levelP,
    const bool above
)
{
    const fvMesh& mesh = levelC.mesh();

    tmp<volScalarField> tResult
    (
        volScalarField::New
        (
            "levelSetFraction",
            mesh,
            dimensionedScalar(dimless, 0)
        )
    );
    volScalarField& result = tResult.ref();

    result.primitiveFieldRef() =
        levelSetFraction
        (
            mesh,
            levelC.primitiveField(),
            levelP.primitiveField(),
            above
        );

    forAll(mesh.boundary(), patchi)
    {
        result.boundaryFieldRef()[patchi] =
            levelSetFraction
            (
                mesh.boundary()[patchi],
                levelC.boundaryField()[patchi],
                levelP.boundaryField()[patchi].patchInternalField()(),
                above
            );
    }

    return tResult;
}


// ************************************************************************* //
