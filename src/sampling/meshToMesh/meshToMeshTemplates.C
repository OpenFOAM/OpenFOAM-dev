/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2022 OpenFOAM Foundation
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

#include "meshToMesh.H"
#include "fvMesh.H"
#include "volFields.H"
#include "directFvPatchFieldMapper.H"
#include "calculatedFvPatchField.H"
#include "patchToPatchFvPatchFieldMapper.H"

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Type>
void Foam::meshToMesh::add
(
    UList<Type>& fld,
    const label offset
) const
{
    forAll(fld, i)
    {
        fld[i] += offset;
    }
}


template<class Type>
void Foam::meshToMesh::mapSrcToTgt
(
    const UList<Type>& srcField,
    List<Type>& result
) const
{
    if (result.size() != tgtToSrcCellAddr_.size())
    {
        FatalErrorInFunction
            << "Supplied field size is not equal to target mesh size" << nl
            << "    source mesh    = " << srcToTgtCellAddr_.size() << nl
            << "    target mesh    = " << tgtToSrcCellAddr_.size() << nl
            << "    supplied field = " << result.size()
            << abort(FatalError);
    }

    if (singleMeshProc_ == -1)
    {
        const distributionMap& map = srcMapPtr_();

        List<Type> work(srcField);
        map.distribute(work);

        forAll(result, celli)
        {
            const labelList& srcAddress = tgtToSrcCellAddr_[celli];
            const scalarList& srcWeight = tgtToSrcCellWght_[celli];

            if (srcAddress.size())
            {
                result[celli] *= (1.0 - sum(srcWeight));
                forAll(srcAddress, i)
                {
                    result[celli] += srcWeight[i]*work[srcAddress[i]];
                }
            }
        }
    }
    else
    {
        forAll(result, celli)
        {
            const labelList& srcAddress = tgtToSrcCellAddr_[celli];
            const scalarList& srcWeight = tgtToSrcCellWght_[celli];

            if (srcAddress.size())
            {
                result[celli] *= (1.0 - sum(srcWeight));
                forAll(srcAddress, i)
                {
                    result[celli] += srcWeight[i]*srcField[srcAddress[i]];
                }
            }
        }
    }
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::meshToMesh::mapSrcToTgt
(
    const Field<Type>& srcField
) const
{
    tmp<Field<Type>> tresult
    (
        new Field<Type>
        (
            tgtToSrcCellAddr_.size(),
            Zero
        )
    );

    mapSrcToTgt(srcField, tresult.ref());

    return tresult;
}


template<class Type>
void Foam::meshToMesh::mapTgtToSrc
(
    const UList<Type>& tgtField,
    List<Type>& result
) const
{
    if (result.size() != srcToTgtCellAddr_.size())
    {
        FatalErrorInFunction
            << "Supplied field size is not equal to source mesh size" << nl
            << "    source mesh    = " << srcToTgtCellAddr_.size() << nl
            << "    target mesh    = " << tgtToSrcCellAddr_.size() << nl
            << "    supplied field = " << result.size()
            << abort(FatalError);
    }

    if (singleMeshProc_ == -1)
    {
        const distributionMap& map = tgtMapPtr_();

        List<Type> work(tgtField);
        map.distribute(work);

        forAll(result, celli)
        {
            const labelList& tgtAddress = srcToTgtCellAddr_[celli];
            const scalarList& tgtWeight = srcToTgtCellWght_[celli];

            if (tgtAddress.size())
            {
                result[celli] *= (1.0 - sum(tgtWeight));
                forAll(tgtAddress, i)
                {
                    result[celli] += tgtWeight[i]*work[tgtAddress[i]];
                }
            }
        }
    }
    else
    {
        forAll(result, celli)
        {
            const labelList& tgtAddress = srcToTgtCellAddr_[celli];
            const scalarList& tgtWeight = srcToTgtCellWght_[celli];

            if (tgtAddress.size())
            {
                result[celli] *= (1.0 - sum(tgtWeight));
                forAll(tgtAddress, i)
                {
                    result[celli] += tgtWeight[i]*tgtField[tgtAddress[i]];
                }
            }
        }
    }
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::meshToMesh::mapTgtToSrc
(
    const Field<Type>& tgtField
) const
{
    tmp<Field<Type>> tresult
    (
        new Field<Type>
        (
            srcToTgtCellAddr_.size(),
            Zero
        )
    );

    mapTgtToSrc(tgtField, tresult.ref());

    return tresult;
}


template<class Type>
void Foam::meshToMesh::mapSrcToTgt
(
    const GeometricField<Type, fvPatchField, volMesh>& field,
    GeometricField<Type, fvPatchField, volMesh>& result
) const
{
    mapSrcToTgt(field, result.primitiveFieldRef());

    typename GeometricField<Type, fvPatchField, volMesh>::
        Boundary& resultBf = result.boundaryFieldRef();

    forAll(patchToPatches(), i)
    {
        label srcPatchi = srcPatchID_[i];
        label tgtPatchi = tgtPatchID_[i];

        const fvPatchField<Type>& srcField = field.boundaryField()[srcPatchi];
        fvPatchField<Type>& tgtField = resultBf[tgtPatchi];

        // Clone and map (since rmap does not do general mapping)
        tmp<fvPatchField<Type>> tnewTgt
        (
            fvPatchField<Type>::New
            (
                srcField,
                tgtField.patch(),
                result(),
                patchToPatchFvPatchFieldMapper(patchToPatches()[i], true)
            )
        );

        // Transfer all mapped quantities (value and e.g. gradient) onto
        // tgtField. Value will get overwritten below.
        tgtField.rmap(tnewTgt(), identity(tgtField.size()));
    }

    forAll(cuttingPatches_, i)
    {
        label patchi = cuttingPatches_[i];
        fvPatchField<Type>& pf = resultBf[patchi];
        pf == pf.patchInternalField();
    }
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh>>
Foam::meshToMesh::mapSrcToTgt
(
    const GeometricField<Type, fvPatchField, volMesh>& field
) const
{
    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;

    const fvMesh& tgtMesh = static_cast<const fvMesh&>(tgtRegion_);

    const fvBoundaryMesh& tgtBm = tgtMesh.boundary();
    const typename fieldType::Boundary& srcBfld =
        field.boundaryField();

    PtrList<fvPatchField<Type>> tgtPatchFields(tgtBm.size());

    // construct tgt boundary patch types as copy of 'field' boundary types
    // note: this will provide place holders for fields with additional
    // entries, but these values will need to be reset
    forAll(tgtPatchID_, i)
    {
        label srcPatchi = srcPatchID_[i];
        label tgtPatchi = tgtPatchID_[i];

        if (!tgtPatchFields.set(tgtPatchi))
        {
            tgtPatchFields.set
            (
                tgtPatchi,
                fvPatchField<Type>::New
                (
                    srcBfld[srcPatchi],
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

    // Any unset tgtPatchFields become calculated
    forAll(tgtPatchFields, tgtPatchi)
    {
        if (!tgtPatchFields.set(tgtPatchi))
        {
            // Note: use factory New method instead of direct generation of
            //       calculated so we keep constraints
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
        }
    }

    tmp<fieldType> tresult
    (
        new fieldType
        (
            IOobject
            (
                typedName("interpolate(" + field.name() + ")"),
                tgtMesh.time().timeName(),
                tgtMesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            tgtMesh,
            field.dimensions(),
            Field<Type>(tgtMesh.nCells(), Zero),
            tgtPatchFields
        )
    );

    mapSrcToTgt(field, tresult.ref());

    return tresult;
}


// ************************************************************************* //
