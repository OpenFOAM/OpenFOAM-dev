/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2018 OpenFOAM Foundation
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

#include "pointConstraints.H"
#include "pointFields.H"
#include "valuePointPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type, class CombineOp>
void Foam::pointConstraints::syncUntransformedData
(
    const polyMesh& mesh,
    List<Type>& pointData,
    const CombineOp& cop
)
{
    // Transfer onto coupled patch
    const globalMeshData& gmd = mesh.globalData();
    const indirectPrimitivePatch& cpp = gmd.coupledPatch();
    const labelList& meshPoints = cpp.meshPoints();

    const mapDistribute& slavesMap = gmd.globalCoPointSlavesMap();
    const labelListList& slaves = gmd.globalCoPointSlaves();

    List<Type> elems(slavesMap.constructSize());
    forAll(meshPoints, i)
    {
        elems[i] = pointData[meshPoints[i]];
    }

    // Pull slave data onto master. No need to update transformed slots.
    slavesMap.distribute(elems, false);

    // Combine master data with slave data
    forAll(slaves, i)
    {
        Type& elem = elems[i];

        const labelList& slavePoints = slaves[i];

        // Combine master with untransformed slave data
        forAll(slavePoints, j)
        {
            cop(elem, elems[slavePoints[j]]);
        }

        // Copy result back to slave slots
        forAll(slavePoints, j)
        {
            elems[slavePoints[j]] = elem;
        }
    }

    // Push slave-slot data back to slaves
    slavesMap.reverseDistribute(elems.size(), elems, false);

    // Extract back onto mesh
    forAll(meshPoints, i)
    {
        pointData[meshPoints[i]] = elems[i];
    }
}


template<class Type>
void Foam::pointConstraints::setPatchFields
(
    GeometricField<Type, pointPatchField, pointMesh>& pf
)
{
    typename GeometricField<Type, pointPatchField, pointMesh>::
        Boundary& pfbf = pf.boundaryFieldRef();

    forAll(pfbf, patchi)
    {
        pointPatchField<Type>& ppf = pfbf[patchi];

        if (isA<valuePointPatchField<Type>>(ppf))
        {
            refCast<valuePointPatchField<Type>>(ppf) =
                ppf.patchInternalField();
        }
    }
}


template<class Type>
void Foam::pointConstraints::constrainCorners
(
    GeometricField<Type, pointPatchField, pointMesh>& pf
) const
{
    forAll(patchPatchPointConstraintPoints_, pointi)
    {
        pf[patchPatchPointConstraintPoints_[pointi]] = transform
        (
            patchPatchPointConstraintTensors_[pointi],
            pf[patchPatchPointConstraintPoints_[pointi]]
        );
    }
}


template<class Type>
void Foam::pointConstraints::constrain
(
    GeometricField<Type, pointPatchField, pointMesh>& pf,
    const bool overrideFixedValue
) const
{
    // Override constrained pointPatchField types with the constraint value.
    // This relies on only constrained pointPatchField implementing the evaluate
    // function
    pf.correctBoundaryConditions();

    // Sync any dangling points
    syncUntransformedData
    (
        mesh()(),
        pf.primitiveFieldRef(),
        maxMagSqrEqOp<Type>()
    );

    // Apply multiple constraints on edge/corner points
    constrainCorners(pf);

    if (overrideFixedValue)
    {
        setPatchFields(pf);
    }
}


// ************************************************************************* //
