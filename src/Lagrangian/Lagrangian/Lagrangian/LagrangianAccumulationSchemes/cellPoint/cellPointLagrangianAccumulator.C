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
#include "syncTools.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cellPointLagrangianAccumulator, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::cellPointLagrangianAccumulator::calcPointCellWeights()
{
    const labelListList& pointCells = mesh().pointCells();

    List<scalar> pointWeights(mesh().nPoints(), scalar(0));

    pointCellWeights_.setSize(pointCells);

    forAll(pointCells, pointi)
    {
        forAll(pointCells[pointi], pointCelli)
        {
            const label celli = pointCells[pointi][pointCelli];

            const scalar w =
                1/mag(mesh().points()[pointi] - mesh().cellCentres()[celli]);

            pointWeights[pointi] += w;

            pointCellWeights_[pointi][pointCelli] = w;
        }
    }

    syncTools::syncPointList
    (
        mesh(),
        pointWeights,
        plusEqOp<scalar>(),
        scalar(0)
    );

    forAll(pointCells, pointi)
    {
        forAll(pointCells[pointi], pointCelli)
        {
            pointCellWeights_[pointi][pointCelli] /= pointWeights[pointi];
        }
    }
}


#define ACCESS_ACCUMULATING_POINT_TYPES(Type, nullArg)                         \
namespace Foam                                                                 \
{                                                                              \
    template<>                                                                 \
    DynamicList<Type>&                                                         \
    cellPointLagrangianAccumulator::accumulatingPointValues() const            \
    {                                                                          \
        autoPtr<DynamicList<Type>>& ptr =                                      \
            CAT3(accumulatingPoint, CAPITALIZE(Type), ValuesPtr_);             \
                                                                               \
        if (!ptr.valid())                                                      \
        {                                                                      \
            ptr.set(new DynamicList<Type>());                                  \
        }                                                                      \
                                                                               \
        return ptr();                                                          \
    }                                                                          \
}
FOR_ALL_FIELD_TYPES(ACCESS_ACCUMULATING_POINT_TYPES)
#undef ACCESS_ACCUMULATING_POINT_TYPES


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cellPointLagrangianAccumulator::cellPointLagrangianAccumulator
(
    const polyMesh& mesh
)
:
    DemandDrivenMeshObject
    <
        polyMesh,
        MoveableMeshObject,
        cellPointLagrangianAccumulator
    >(mesh),
    pointCellWeights_(),
    pointAccumulatingPoint_(mesh.nPoints(), -1),
    accumulatingPointPoint_()
{
    calcPointCellWeights();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::cellPointLagrangianAccumulator::~cellPointLagrangianAccumulator()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::cellPointLagrangianAccumulator::movePoints()
{
    calcPointCellWeights();
    return true;
}


// ************************************************************************* //
