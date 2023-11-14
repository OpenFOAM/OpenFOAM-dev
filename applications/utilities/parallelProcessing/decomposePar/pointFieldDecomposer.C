/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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

#include "pointFieldDecomposer.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::labelList Foam::pointFieldDecomposer::patchFieldDecomposer::addressing
(
    const pointPatch& completePatch,
    const pointPatch& procPatch,
    const labelList& pointProcAddressing
)
{
    const labelList& completePatchPoints = completePatch.meshPoints();
    const labelList& procPatchPoints = procPatch.meshPoints();

    // Create a map from complete mesh point index to complete patch point index
    labelList map(completePatch.boundaryMesh().mesh().size(), -1);
    forAll(completePatchPoints, pointi)
    {
        map[completePatchPoints[pointi]] = pointi;
    }

    // Determine the complete patch point for every proc patch point, going via
    // the complete mesh point index and using the above map
    labelList result(procPatch.size(), -1);
    forAll(procPatchPoints, pointi)
    {
        result[pointi] = map[pointProcAddressing[procPatchPoints[pointi]]];
    }

    // Check that all the patch point addresses are set
    if (result.size() && min(result) < 0)
    {
        FatalErrorInFunction
            << "Incomplete patch point addressing"
            << abort(FatalError);
    }

    return result;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pointFieldDecomposer::patchFieldDecomposer::patchFieldDecomposer
(
    const pointPatch& completePatch,
    const pointPatch& procPatch,
    const labelList& pointProcAddressing
)
:
    labelList(addressing(completePatch, procPatch, pointProcAddressing)),
    forwardFieldMapper(static_cast<const labelList&>(*this))
{}


Foam::pointFieldDecomposer::pointFieldDecomposer
(
    const pointMesh& completeMesh,
    const PtrList<fvMesh>& procMeshes,
    const labelListList& pointProcAddressing
)
:
    completeMesh_(completeMesh),
    procMeshes_(procMeshes),
    pointProcAddressing_(pointProcAddressing),
    patchFieldDecomposers_(procMeshes_.size())
{
    forAll(procMeshes_, proci)
    {
        const pointMesh& procMesh = pointMesh::New(procMeshes_[proci]);

        patchFieldDecomposers_.set
        (
            proci,
            new PtrList<patchFieldDecomposer>(procMesh.boundary().size())
        );

        forAll(procMesh.boundary(), procPatchi)
        {
            if (procPatchi < completeMesh_.boundary().size())
            {
                patchFieldDecomposers_[proci].set
                (
                    procPatchi,
                    new patchFieldDecomposer
                    (
                        completeMesh_.boundary()[procPatchi],
                        procMesh.boundary()[procPatchi],
                        pointProcAddressing_[proci]
                    )
                );
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::pointFieldDecomposer::~pointFieldDecomposer()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::pointFieldDecomposer::decomposes(const IOobjectList& objects)
{
    bool result = false;

    #define DO_POINT_FIELDS_TYPE(Type, nullArg)                                \
        result = result                                                        \
         || !objects.lookupClass(PointField<Type>::typeName).empty();
    FOR_ALL_FIELD_TYPES(DO_POINT_FIELDS_TYPE)
    #undef DO_POINT_FIELDS_TYPE

    return result;
}


// ************************************************************************* //
