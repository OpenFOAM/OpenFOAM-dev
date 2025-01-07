/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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

#include "cellPointLagrangianAccumulator.H"
#include "cellPointLagrangianAddressor.H"
#include "LagrangianFields.H"
#include "LagrangianSubFields.H"
#include "tetIndices.H"
#include "syncTools.H"

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Type>
void Foam::cellPointLagrangianAccumulator::accumulate
(
    const LagrangianSubSubField<Type>& lPsi,
    Field<Type>& cPsi
) const
{
    const LagrangianSubMesh& lSubMesh = lPsi.mesh();
    const LagrangianMesh& lMesh = lSubMesh.mesh();

    // Do the simple cell-cell contributions
    forAll(lPsi.mesh(), subi)
    {
        const label i = lSubMesh.start() + subi;

        const barycentric& coordinates = lMesh.coordinates()[i];
        const label celli = lMesh.celli()[i];

        cPsi[celli] += coordinates.a()*lPsi[subi];
    }

    // Do the more complicated cell-point-cell contributions...

    // Sum contributions into the point workspace
    DynamicList<Type>& accumulatingPointValues =
        this->accumulatingPointValues<Type>();
    forAll(lPsi.mesh(), subi)
    {
        const label i = lSubMesh.start() + subi;

        const barycentric& coordinates = lMesh.coordinates()[i];
        const label celli = lMesh.celli()[i];
        const label facei = lMesh.facei()[i];
        const label faceTrii = lMesh.faceTrii()[i];

        const triFace triPoints =
            tetIndices(celli, facei, faceTrii).faceTriIs(mesh());

        forAll(triPoints, triPointi)
        {
            const label pointi = triPoints[triPointi];

            if (pointAccumulatingPoint_[pointi] == -1)
            {
                pointAccumulatingPoint_[pointi] =
                    accumulatingPointPoint_.size();

                accumulatingPointPoint_.append(pointi);
                accumulatingPointValues.append(pTraits<Type>::zero);
            }

            accumulatingPointValues[pointAccumulatingPoint_[pointi]] +=
                coordinates[triPointi + 1]*lPsi[subi];
        }
    }

    // Expand to include any points that are being accumulated to remotely
    cellPointLagrangianAddressor::New(mesh()).sync
    (
        pointAccumulatingPoint_,
        accumulatingPointPoint_
    );
    accumulatingPointValues.resize
    (
        accumulatingPointPoint_.size(),
        pTraits<Type>::zero
    );

    // Synchronise
    syncTools::syncPointList
    (
        mesh(),
        accumulatingPointPoint_,
        accumulatingPointValues,
        plusEqOp<Type>(),
        pTraits<Type>::zero
    );

    // Accumulate back into the cells
    const labelListList& pointCells = mesh().pointCells();
    forAll(accumulatingPointPoint_, accumulatingPointi)
    {
        const label pointi = accumulatingPointPoint_[accumulatingPointi];

        forAll(pointCells[pointi], pointCelli)
        {
            const label celli = pointCells[pointi][pointCelli];

            cPsi[celli] +=
                pointCellWeights_[pointi][pointCelli]
               *accumulatingPointValues[accumulatingPointi];
        }
    }

    // Reset the workspace
    forAll(accumulatingPointPoint_, accumulatingPointi)
    {
        const label pointi = accumulatingPointPoint_[accumulatingPointi];

        pointAccumulatingPoint_[pointi] = -1;
    }

    accumulatingPointPoint_.clear();
    accumulatingPointValues.clear();
}


// ************************************************************************* //
