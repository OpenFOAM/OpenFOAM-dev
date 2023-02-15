/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2023 OpenFOAM Foundation
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

#include "fvMeshToFvMesh.H"
#include "directFvPatchFieldMapper.H"
#include "identityFvPatchFieldMapper.H"
#include "patchToPatchFvPatchFieldMapper.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::fvMeshToFvMesh::evaluateConstraintTypes(VolField<Type>& fld)
{
    typename VolField<Type>::Boundary& fldBf = fld.boundaryFieldRef();

    if
    (
        Pstream::defaultCommsType == Pstream::commsTypes::blocking
     || Pstream::defaultCommsType == Pstream::commsTypes::nonBlocking
    )
    {
        label nReq = Pstream::nRequests();

        forAll(fldBf, patchi)
        {
            fvPatchField<Type>& tgtField = fldBf[patchi];

            if
            (
                tgtField.type() == tgtField.patch().patch().type()
             && polyPatch::constraintType(tgtField.patch().patch().type())
            )
            {
                tgtField.initEvaluate(Pstream::defaultCommsType);
            }
        }

        // Block for any outstanding requests
        if
        (
            Pstream::parRun()
         && Pstream::defaultCommsType == Pstream::commsTypes::nonBlocking
        )
        {
            Pstream::waitRequests(nReq);
        }

        forAll(fldBf, patchi)
        {
            fvPatchField<Type>& tgtField = fldBf[patchi];

            if
            (
                tgtField.type() == tgtField.patch().patch().type()
             && polyPatch::constraintType(tgtField.patch().patch().type())
            )
            {
                tgtField.evaluate(Pstream::defaultCommsType);
            }
        }
    }
    else if (Pstream::defaultCommsType == Pstream::commsTypes::scheduled)
    {
        const lduSchedule& patchSchedule =
            fld.mesh().globalData().patchSchedule();

        forAll(patchSchedule, patchEvali)
        {
            label patchi = patchSchedule[patchEvali].patch;
            fvPatchField<Type>& tgtField = fldBf[patchi];

            if
            (
                tgtField.type() == tgtField.patch().patch().type()
             && polyPatch::constraintType(tgtField.patch().patch().type())
            )
            {
                if (patchSchedule[patchEvali].init)
                {
                    tgtField.initEvaluate(Pstream::commsTypes::scheduled);
                }
                else
                {
                    tgtField.evaluate(Pstream::commsTypes::scheduled);
                }
            }
        }
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::VolField<Type>> Foam::fvMeshToFvMesh::srcToTgt
(
    const VolField<Type>& srcFld
) const
{
    const fvMesh& tgtMesh = static_cast<const fvMesh&>(meshToMesh::tgtMesh());

    // Construct target patch fields as copies of source patch fields, but do
    // not map values yet
    PtrList<fvPatchField<Type>> tgtPatchFields(tgtMesh.boundary().size());
    forAll(patchIDs(), i)
    {
        const label srcPatchi = patchIDs()[i].first();
        const label tgtPatchi = patchIDs()[i].second();

        if (!tgtPatchFields.set(tgtPatchi))
        {
            tgtPatchFields.set
            (
                tgtPatchi,
                fvPatchField<Type>::New
                (
                    srcFld.boundaryField()[srcPatchi],
                    tgtMesh.boundary()[tgtPatchi],
                    DimensionedField<Type, volMesh>::null(),
                    directFvPatchFieldMapper
                    (
                        labelList(tgtMesh.boundary()[tgtPatchi].size(), -1)
                    )
                )
            );
        }
    }

    // Create calculated patch fields for any unset target patches. Use
    // fvPatchField<Type>::New to construct these, rather than using the
    // calculated constructor directly, so that constraints are maintained.
    labelList tgtPatchFieldIsUnMapped(tgtPatchFields.size(), false);
    forAll(tgtPatchFields, tgtPatchi)
    {
        if (!tgtPatchFields.set(tgtPatchi))
        {
            tgtPatchFields.set
            (
                tgtPatchi,
                fvPatchField<Type>::New
                (
                    calculatedFvPatchField<Type>::typeName,
                    tgtMesh.boundary()[tgtPatchi],
                    DimensionedField<Type, volMesh>::null()
                )
            );

            tgtPatchFieldIsUnMapped[tgtPatchi] =
                polyPatch::constraintType
                (
                    tgtMesh.boundary()[tgtPatchi].patch().type()
                );
        }
    }

    // Construct the result field
    tmp<VolField<Type>> ttgtFld =
        VolField<Type>::New
        (
            typedName("interpolate(" + srcFld.name() + ")"),
            srcToTgt<Type>(srcFld.internalField())(),
            tgtPatchFields
        );
    typename VolField<Type>::Boundary& tgtBfld =
        ttgtFld.ref().boundaryFieldRef();

    // Mapped patches
    forAll(patchIDs(), i)
    {
        const label srcPatchi = patchIDs()[i].first();
        const label tgtPatchi = patchIDs()[i].second();

        tgtBfld[tgtPatchi].map
        (
            srcFld.boundaryField()[srcPatchi],
            patchToPatchNormalisedFvPatchFieldMapper
            (
                patchInterpolation(i),
                tgtPatchStabilisation(i)
            )
        );
    }

    // Un-mapped patches. Set values to that of the internal cell field.
    forAll(tgtBfld, patchi)
    {
        if (tgtPatchFieldIsUnMapped[patchi])
        {
            fvPatchField<Type>& tgtPfld = tgtBfld[patchi];
            tgtPfld == tgtPfld.patchInternalField();
        }
    }

    // Evaluate constraints
    evaluateConstraintTypes(ttgtFld.ref());

    return ttgtFld;
}


template<class Type>
Foam::tmp<Foam::VolField<Type>> Foam::fvMeshToFvMesh::srcToTgt
(
    const VolField<Type>& srcFld,
    const VolField<Type>& leftOverTgtFld,
    const UList<wordRe>& tgtCuttingPatchNames
) const
{
    // Construct the result field
    tmp<VolField<Type>> ttgtFld =
        VolField<Type>::New
        (
            typedName("interpolate(" + srcFld.name() + ")"),
            srcToTgt<Type>(srcFld.v(), leftOverTgtFld.v())(),
            leftOverTgtFld.boundaryField()
        );
    typename VolField<Type>::Boundary& tgtBfld =
        ttgtFld.ref().boundaryFieldRef();

    // Mapped patches
    forAll(patchIDs(), i)
    {
        const label srcPatchi = patchIDs()[i].first();
        const label tgtPatchi = patchIDs()[i].second();

        tgtBfld[tgtPatchi].map
        (
            leftOverTgtFld.boundaryField()[tgtPatchi],
            identityFvPatchFieldMapper()
        );
        tgtBfld[tgtPatchi].map
        (
            srcFld.boundaryField()[srcPatchi],
            patchToPatchLeftOverFvPatchFieldMapper(patchInterpolation(i))
        );
    }

    // Cutting patches. Set values to that of the internal cell field.
    const labelHashSet tgtCuttingPatchIDs =
        leftOverTgtFld.mesh().boundaryMesh().patchSet(tgtCuttingPatchNames);
    forAllConstIter(labelHashSet, tgtCuttingPatchIDs, iter)
    {
        tgtBfld[iter.key()] == tgtBfld[iter.key()].patchInternalField();
    }

    // Evaluate constraints
    evaluateConstraintTypes(ttgtFld.ref());

    return ttgtFld;
}


template<class Type>
Foam::tmp<Foam::VolInternalField<Type>> Foam::fvMeshToFvMesh::srcToTgt
(
    const VolInternalField<Type>& srcFld
) const
{
    tmp<VolInternalField<Type>> ttgtFld =
        VolInternalField<Type>::New
        (
            typedName("interpolate(" + srcFld.name() + ")"),
            static_cast<const fvMesh&>(meshToMesh::tgtMesh()),
            srcFld.dimensions(),
            cellsInterpolation().srcToTgt(srcFld)
        );

    tgtCellsStabilisation().stabilise(ttgtFld.ref());

    return ttgtFld;
}


template<class Type>
Foam::tmp<Foam::VolInternalField<Type>> Foam::fvMeshToFvMesh::srcToTgt
(
    const VolInternalField<Type>& srcFld,
    const VolInternalField<Type>& leftOverTgtFld
) const
{
    return
        VolInternalField<Type>::New
        (
            typedName("interpolate(" + srcFld.name() + ")"),
            static_cast<const fvMesh&>(meshToMesh::tgtMesh()),
            leftOverTgtFld.dimensions(),
            cellsInterpolation().srcToTgt(srcFld, leftOverTgtFld)
        );
}


// ************************************************************************* //
