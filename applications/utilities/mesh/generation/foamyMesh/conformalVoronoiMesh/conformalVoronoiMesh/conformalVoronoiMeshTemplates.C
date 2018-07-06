/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
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

#include "conformalVoronoiMesh.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Triangulation>
bool Foam::conformalVoronoiMesh::distributeBackground(const Triangulation& mesh)
{
    if (!Pstream::parRun())
    {
        return false;
    }

    Info<< nl << "Redistributing points" << endl;

    timeCheck("Before distribute");

    label iteration = 0;

    scalar previousLoadUnbalance = 0;

    while (true)
    {
        scalar maxLoadUnbalance = mesh.calculateLoadUnbalance();

        if
        (
            maxLoadUnbalance <= foamyHexMeshControls().maxLoadUnbalance()
         || maxLoadUnbalance <= previousLoadUnbalance
        )
        {
            // If this is the first iteration, return false, if it was a
            // subsequent one, return true;
            return iteration != 0;
        }

        previousLoadUnbalance = maxLoadUnbalance;

        Info<< "    Total number of vertices before redistribution "
            << returnReduce(label(mesh.number_of_vertices()), sumOp<label>())
            << endl;

        const fvMesh& bMesh = decomposition_().mesh();

        volScalarField cellWeights
        (
            IOobject
            (
                "cellWeights",
                bMesh.time().timeName(),
                bMesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            bMesh,
            dimensionedScalar("weight", dimless, 1e-2),
            zeroGradientFvPatchScalarField::typeName
        );

        meshSearch cellSearch(bMesh, polyMesh::FACE_PLANES);

        labelList cellVertices(bMesh.nCells(), label(0));

        for
        (
            typename Triangulation::Finite_vertices_iterator vit
                = mesh.finite_vertices_begin();
            vit != mesh.finite_vertices_end();
            ++vit
        )
        {
            // Only store real vertices that are not feature vertices
            if (vit->real() && !vit->featurePoint())
            {
                pointFromPoint v = topoint(vit->point());

                label celli = cellSearch.findCell(v);

                if (celli == -1)
                {
//                     Pout<< "findCell conformalVoronoiMesh::distribute "
//                         << "findCell "
//                         << vit->type() << " "
//                         << vit->index() << " "
//                         << v << " "
//                         << celli
//                         << " find nearest celli ";

                    celli = cellSearch.findNearestCell(v);
                }

                cellVertices[celli]++;
            }
        }

        scalarField& cwi = cellWeights.primitiveFieldRef();

        forAll(cellVertices, cI)
        {
            // Give a small but finite weight for empty cells.  Some
            // decomposition methods have difficulty with integer overflows in
            // the sum of the normalised weight field.
            cwi[cI] = max(cellVertices[cI], 1e-2);
        }

        autoPtr<mapDistributePolyMesh> mapDist = decomposition_().distribute
        (
            cellWeights
        );

        cellShapeControl_.shapeControlMesh().distribute(decomposition_);

        distribute();

        timeCheck("After distribute");

        iteration++;
    }

    return true;
}


// ************************************************************************* //
