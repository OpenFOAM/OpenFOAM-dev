/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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

#include "zonalBoundaryIndices.H"
#include "surfaceFields.H"
#include "cyclicFvPatch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace blendingMethods
{
    defineTypeNameAndDebug(zonal::boundaryIndices, 0);
}
}


// * * * * * * * * * * * * * * Protected Constructors  * * * * * * * * * * * //

Foam::blendingMethods::zonal::boundaryIndices::boundaryIndices
(
    const word& name,
    const fvMesh& mesh,
    const labelList& zoneIndices
)
:
    DemandDrivenMeshObject<fvMesh, MoveableMeshObject, boundaryIndices>
    (
        name,
        mesh
    ),
    List<List<label>>(mesh.boundary().size())
{
    surfaceLabelField::Boundary zoneIndexOwnerBf
    (
        mesh.boundary(),
        surfaceLabelField::Internal::null(),
        calculatedFvsPatchField<label>::typeName
    );

    zoneIndexOwnerBf = -1;

    forAll(zoneIndices_, zonei)
    {
        const cellZone& zoneCells =
            mesh.cellZones()[zoneIndices_[zonei]];

        boolList cellInZone(mesh.nCells(), false);
        UIndirectList<bool>(cellInZone, zoneCells) = true;

        forAll(zoneIndexOwnerBf, patchi)
        {
            const fvPatch& fvp = mesh.boundary()[patchi];

            const cyclicFvPatch& cfvp =
                refCastNull<const cyclicFvPatch>(fvp);

            const labelList& patchFaceCells =
                notNull(cfvp) && cfvp.neighbour()
              ? cfvp.nbrPatch().faceCells()
              : fvp.faceCells();

            forAll(zoneIndexOwnerBf[patchi], patchFacei)
            {
                if (cellInZone[patchFaceCells[patchFacei]])
                {
                    zoneIndexOwnerBf[patchi][patchFacei] = zonei;
                }
            }
        }
    }

    PtrList<labelField> zoneIndexNeighbourBf =
        zoneIndexOwnerBf.coupledNeighbourField();

    forAll(mesh.boundary(), patchi)
    {
        if (zoneIndexNeighbourBf.set(patchi))
        {
            operator[](patchi).transfer(zoneIndexNeighbourBf[patchi]);
        }
        else
        {
            operator[](patchi).transfer(zoneIndexOwnerBf[patchi]);
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

const Foam::blendingMethods::zonal::boundaryIndices&
Foam::blendingMethods::zonal::boundaryIndices::New
(
    const fvMesh& mesh,
    const labelList& zoneIndices
)
{
    return DemandDrivenMeshObject
    <
        fvMesh,
        MoveableMeshObject,
        boundaryIndices
    >::New
    (
        Foam::name(Hash<labelList>()(zoneIndices)),
        mesh,
        zoneIndices
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::blendingMethods::zonal::boundaryIndices::~boundaryIndices()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::blendingMethods::zonal::boundaryIndices::movePoints()
{
    return true;
}


// ************************************************************************* //
