/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2014 OpenFOAM Foundation
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

#include "fvMesh.H"
#include "volFields.H"
#include "directFvPatchFieldMapper.H"
#include "calculatedFvPatchField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    //- Helper class for list
    template<class Type>
    class ListPlusEqOp
    {
        public:
        void operator()(List<Type>& x, const List<Type> y) const
        {
            if (y.size())
            {
                if (x.size())
                {
                    label sz = x.size();
                    x.setSize(sz + y.size());
                    forAll(y, i)
                    {
                        x[sz++] = y[i];
                    }
                }
                else
                {
                    x = y;
                }
            }
        }
    };
}


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


template<class Type, class CombineOp>
void Foam::meshToMesh::mapSrcToTgt
(
    const UList<Type>& srcField,
    const CombineOp& cop,
    List<Type>& result
) const
{
    if (result.size() != tgtToSrcCellAddr_.size())
    {
        FatalErrorIn
        (
            "void Foam::meshToMesh::mapSrcToTgt"
            "("
                "const UList<Type>&, "
                "const CombineOp&, "
                "List<Type>&"
            ") const"
        )   << "Supplied field size is not equal to target mesh size" << nl
            << "    source mesh    = " << srcToTgtCellAddr_.size() << nl
            << "    target mesh    = " << tgtToSrcCellAddr_.size() << nl
            << "    supplied field = " << result.size()
            << abort(FatalError);
    }

    multiplyWeightedOp<Type, CombineOp> cbop(cop);

    if (singleMeshProc_ == -1)
    {
        const mapDistribute& map = srcMapPtr_();

        List<Type> work(srcField);
        map.distribute(work);

        forAll(result, cellI)
        {
            const labelList& srcAddress = tgtToSrcCellAddr_[cellI];
            const scalarList& srcWeight = tgtToSrcCellWght_[cellI];

            if (srcAddress.size())
            {
//                result[cellI] = pTraits<Type>::zero;
                result[cellI] *= (1.0 - sum(srcWeight));
                forAll(srcAddress, i)
                {
                    label srcI = srcAddress[i];
                    scalar w = srcWeight[i];
                    cbop(result[cellI], cellI, work[srcI], w);
                }
            }
        }
    }
    else
    {
        forAll(result, cellI)
        {
            const labelList& srcAddress = tgtToSrcCellAddr_[cellI];
            const scalarList& srcWeight = tgtToSrcCellWght_[cellI];

            if (srcAddress.size())
            {
//                result[cellI] = pTraits<Type>::zero;
                result[cellI] *= (1.0 - sum(srcWeight));
                forAll(srcAddress, i)
                {
                    label srcI = srcAddress[i];
                    scalar w = srcWeight[i];
                    cbop(result[cellI], cellI, srcField[srcI], w);
                }
            }
        }
    }
}


template<class Type, class CombineOp>
Foam::tmp<Foam::Field<Type> > Foam::meshToMesh::mapSrcToTgt
(
    const Field<Type>& srcField,
    const CombineOp& cop
) const
{
    tmp<Field<Type> > tresult
    (
        new Field<Type>
        (
            tgtToSrcCellAddr_.size(),
            pTraits<Type>::zero
        )
    );

    mapSrcToTgt(srcField, cop, tresult());

    return tresult;
}


template<class Type, class CombineOp>
Foam::tmp<Foam::Field<Type> > Foam::meshToMesh::mapSrcToTgt
(
    const tmp<Field<Type> >& tsrcField,
    const CombineOp& cop
) const
{
    return mapSrcToTgt(tsrcField(), cop);
}


template<class Type>
Foam::tmp<Foam::Field<Type> > Foam::meshToMesh::mapSrcToTgt
(
    const Field<Type>& srcField
) const
{
    return mapSrcToTgt(srcField, plusEqOp<Type>());
}


template<class Type>
Foam::tmp<Foam::Field<Type> > Foam::meshToMesh::mapSrcToTgt
(
    const tmp<Field<Type> >& tsrcField
) const
{
    return mapSrcToTgt(tsrcField());
}


template<class Type, class CombineOp>
void Foam::meshToMesh::mapTgtToSrc
(
    const UList<Type>& tgtField,
    const CombineOp& cop,
    List<Type>& result
) const
{
    if (result.size() != srcToTgtCellAddr_.size())
    {
        FatalErrorIn
        (
            "void Foam::meshToMesh::mapTgtToSrc"
            "("
                "const UList<Type>&, "
                "const CombineOp&, "
                "List<Type>&"
            ") const"
        )   << "Supplied field size is not equal to source mesh size" << nl
            << "    source mesh    = " << srcToTgtCellAddr_.size() << nl
            << "    target mesh    = " << tgtToSrcCellAddr_.size() << nl
            << "    supplied field = " << result.size()
            << abort(FatalError);
    }

    multiplyWeightedOp<Type, CombineOp> cbop(cop);

    if (singleMeshProc_ == -1)
    {
        const mapDistribute& map = tgtMapPtr_();

        List<Type> work(tgtField);
        map.distribute(work);

        forAll(result, cellI)
        {
            const labelList& tgtAddress = srcToTgtCellAddr_[cellI];
            const scalarList& tgtWeight = srcToTgtCellWght_[cellI];

            if (tgtAddress.size())
            {
                result[cellI] *= (1.0 - sum(tgtWeight));
                forAll(tgtAddress, i)
                {
                    label tgtI = tgtAddress[i];
                    scalar w = tgtWeight[i];
                    cbop(result[cellI], cellI, work[tgtI], w);
                }
            }
        }
    }
    else
    {
        forAll(result, cellI)
        {
            const labelList& tgtAddress = srcToTgtCellAddr_[cellI];
            const scalarList& tgtWeight = srcToTgtCellWght_[cellI];

            if (tgtAddress.size())
            {
                result[cellI] *= (1.0 - sum(tgtWeight));
                forAll(tgtAddress, i)
                {
                    label tgtI = tgtAddress[i];
                    scalar w = tgtWeight[i];
                    cbop(result[cellI], cellI, tgtField[tgtI], w);
                }
            }
        }
    }
}


template<class Type, class CombineOp>
Foam::tmp<Foam::Field<Type> > Foam::meshToMesh::mapTgtToSrc
(
    const Field<Type>& tgtField,
    const CombineOp& cop
) const
{
    tmp<Field<Type> > tresult
    (
        new Field<Type>
        (
            srcToTgtCellAddr_.size(),
            pTraits<Type>::zero
        )
    );

    mapTgtToSrc(tgtField, cop, tresult());

    return tresult;
}


template<class Type, class CombineOp>
Foam::tmp<Foam::Field<Type> > Foam::meshToMesh::mapTgtToSrc
(
    const tmp<Field<Type> >& ttgtField,
    const CombineOp& cop
) const
{
    return mapTgtToSrc(ttgtField(), cop);
}


template<class Type>
Foam::tmp<Foam::Field<Type> > Foam::meshToMesh::mapTgtToSrc
(
    const Field<Type>& tgtField
) const
{
    return mapTgtToSrc(tgtField, plusEqOp<Type>());
}


template<class Type>
Foam::tmp<Foam::Field<Type> > Foam::meshToMesh::mapTgtToSrc
(
    const tmp<Field<Type> >& ttgtField
) const
{
    return mapTgtToSrc(ttgtField(), plusEqOp<Type>());
}


template<class Type, class CombineOp>
void Foam::meshToMesh::mapSrcToTgt
(
    const GeometricField<Type, fvPatchField, volMesh>& field,
    const CombineOp& cop,
    GeometricField<Type, fvPatchField, volMesh>& result
) const
{
    mapSrcToTgt(field, cop, result.internalField());

    const PtrList<AMIPatchToPatchInterpolation>& AMIList = patchAMIs();

    forAll(AMIList, i)
    {
        label srcPatchI = srcPatchID_[i];
        label tgtPatchI = tgtPatchID_[i];

        const Field<Type>& srcField = field.boundaryField()[srcPatchI];
        Field<Type>& tgtField = result.boundaryField()[tgtPatchI];

        tgtField = pTraits<Type>::zero;

        AMIList[i].interpolateToTarget
        (
            srcField,
            multiplyWeightedOp<Type, CombineOp>(cop),
            tgtField,
            UList<Type>::null()
        );
    }

    forAll(cuttingPatches_, i)
    {
        label patchI = cuttingPatches_[i];
        fvPatchField<Type>& pf = result.boundaryField()[patchI];
        pf == pf.patchInternalField();
    }
}


template<class Type, class CombineOp>
Foam::tmp<Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh> >
Foam::meshToMesh::mapSrcToTgt
(
    const GeometricField<Type, fvPatchField, volMesh>& field,
    const CombineOp& cop
) const
{
    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;

    const fvMesh& tgtMesh = static_cast<const fvMesh&>(tgtRegion_);

    const fvBoundaryMesh& tgtBm = tgtMesh.boundary();
    const typename fieldType::GeometricBoundaryField& srcBfld =
        field.boundaryField();

    PtrList<fvPatchField<Type> > tgtPatchFields(tgtBm.size());

    // constuct tgt boundary patch types as copy of 'field' boundary types
    // note: this will provide place holders for fields with additional
    // entries, but these values will need to be reset
    forAll(tgtPatchID_, i)
    {
        label srcPatchI = srcPatchID_[i];
        label tgtPatchI = tgtPatchID_[i];

        if (!tgtPatchFields.set(tgtPatchI))
        {
            tgtPatchFields.set
            (
                tgtPatchI,
                fvPatchField<Type>::New
                (
                    srcBfld[srcPatchI],
                    tgtMesh.boundary()[tgtPatchI],
                    DimensionedField<Type, volMesh>::null(),
                    directFvPatchFieldMapper
                    (
                        labelList(tgtMesh.boundary()[tgtPatchI].size(), -1)
                    )
                )
            );
        }
    }

    // Any unset tgtPatchFields become calculated
    forAll(tgtPatchFields, tgtPatchI)
    {
        if (!tgtPatchFields.set(tgtPatchI))
        {
            // Note: use factory New method instead of direct generation of
            //       calculated so we keep constraints
            tgtPatchFields.set
            (
                tgtPatchI,
                fvPatchField<Type>::New
                (
                    calculatedFvPatchField<Type>::typeName,
                    tgtMesh.boundary()[tgtPatchI],
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
                type() + ":interpolate(" + field.name() + ")",
                tgtMesh.time().timeName(),
                tgtMesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            tgtMesh,
            field.dimensions(),
            Field<Type>(tgtMesh.nCells(), pTraits<Type>::zero),
            tgtPatchFields
        )
    );

    mapSrcToTgt(field, cop, tresult());

    return tresult;
}


template<class Type, class CombineOp>
Foam::tmp<Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh> >
Foam::meshToMesh::mapSrcToTgt
(
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tfield,
    const CombineOp& cop
) const
{
    return mapSrcToTgt(tfield(), cop);
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh> >
Foam::meshToMesh::mapSrcToTgt
(
    const GeometricField<Type, fvPatchField, volMesh>& field
) const
{
    return mapSrcToTgt(field, plusEqOp<Type>());
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh> >
Foam::meshToMesh::mapSrcToTgt
(
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tfield
) const
{
    return mapSrcToTgt(tfield(), plusEqOp<Type>());
}


template<class Type, class CombineOp>
void Foam::meshToMesh::mapTgtToSrc
(
    const GeometricField<Type, fvPatchField, volMesh>& field,
    const CombineOp& cop,
    GeometricField<Type, fvPatchField, volMesh>& result
) const
{
    mapTgtToSrc(field, cop, result.internalField());

    const PtrList<AMIPatchToPatchInterpolation>& AMIList = patchAMIs();

    forAll(AMIList, i)
    {
        label srcPatchI = srcPatchID_[i];
        label tgtPatchI = tgtPatchID_[i];

        Field<Type>& srcField = result.boundaryField()[srcPatchI];
        const Field<Type>& tgtField = field.boundaryField()[tgtPatchI];

        srcField = pTraits<Type>::zero;

        AMIList[i].interpolateToSource
        (
            tgtField,
            multiplyWeightedOp<Type, CombineOp>(cop),
            srcField
        );
    }

    forAll(cuttingPatches_, i)
    {
        label patchI = cuttingPatches_[i];
        fvPatchField<Type>& pf = result.boundaryField()[patchI];
        pf == pf.patchInternalField();
    }
}


template<class Type, class CombineOp>
Foam::tmp<Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh> >
Foam::meshToMesh::mapTgtToSrc
(
    const GeometricField<Type, fvPatchField, volMesh>& field,
    const CombineOp& cop
) const
{
    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;

    const fvMesh& srcMesh = static_cast<const fvMesh&>(srcRegion_);

    const fvBoundaryMesh& srcBm = srcMesh.boundary();
    const typename fieldType::GeometricBoundaryField& tgtBfld =
        field.boundaryField();

    PtrList<fvPatchField<Type> > srcPatchFields(srcBm.size());

    // constuct src boundary patch types as copy of 'field' boundary types
    // note: this will provide place holders for fields with additional
    // entries, but these values will need to be reset
    forAll(srcPatchID_, i)
    {
        label srcPatchI = srcPatchID_[i];
        label tgtPatchI = tgtPatchID_[i];

        if (!srcPatchFields.set(tgtPatchI))
        {
            srcPatchFields.set
            (
                srcPatchI,
                fvPatchField<Type>::New
                (
                    tgtBfld[srcPatchI],
                    srcMesh.boundary()[tgtPatchI],
                    DimensionedField<Type, volMesh>::null(),
                    directFvPatchFieldMapper
                    (
                        labelList(srcMesh.boundary()[srcPatchI].size(), -1)
                    )
                )
            );
        }
    }

    // Any unset srcPatchFields become calculated
    forAll(srcPatchFields, srcPatchI)
    {
        if (!srcPatchFields.set(srcPatchI))
        {
            // Note: use factory New method instead of direct generation of
            //       calculated so we keep constraints
            srcPatchFields.set
            (
                srcPatchI,
                fvPatchField<Type>::New
                (
                    calculatedFvPatchField<Type>::typeName,
                    srcMesh.boundary()[srcPatchI],
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
                type() + ":interpolate(" + field.name() + ")",
                srcMesh.time().timeName(),
                srcMesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            srcMesh,
            field.dimensions(),
            Field<Type>(srcMesh.nCells(), pTraits<Type>::zero),
            srcPatchFields
        )
    );

    mapTgtToSrc(field, cop, tresult());

    return tresult;
}


template<class Type, class CombineOp>
Foam::tmp<Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh> >
Foam::meshToMesh::mapTgtToSrc
(
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tfield,
    const CombineOp& cop
) const
{
    return mapTgtToSrc(tfield(), cop);
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh> >
Foam::meshToMesh::mapTgtToSrc
(
    const GeometricField<Type, fvPatchField, volMesh>& field
) const
{
    return mapTgtToSrc(field, plusEqOp<Type>());
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh> >
Foam::meshToMesh::mapTgtToSrc
(
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tfield
) const
{
    return mapTgtToSrc(tfield(), plusEqOp<Type>());
}


// ************************************************************************* //
