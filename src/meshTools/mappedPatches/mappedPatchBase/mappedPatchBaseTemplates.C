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

#include "mappedPatchBase.H"
#include "stringOps.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class PatchField>
void Foam::mappedPatchBase::validateMapForField
(
    const PatchField& field,
    const dictionary& context,
    const label froms
)
{
    const polyPatch& pp = field.patch().patch();

    if (!isA<mappedPatchBase>(pp))
    {
        OStringStream str;
        str << "Field " << field.internalField().name() << " of type "
            << field.type() << " cannot apply to patch " << pp.name()
            << " because the patch is not of " << typeName << " type";
        FatalIOErrorInFunction(context)
            << stringOps::breakIntoIndentedLines(str.str()).c_str()
            << exit(FatalIOError);
    }

    refCast<const mappedPatchBase>(pp).validateForField(field, context, froms);
}


template<class PatchField>
void Foam::mappedPatchBase::validateForField
(
    const PatchField& field,
    const dictionary& context,
    const label froms
) const
{
    const bool isNotRegion = !sameRegion() && (froms & from::sameRegion);
    const bool isRegion = sameRegion() && (froms & from::differentRegion);
    const bool isPatch = samePatch() && (froms & from::differentPatch);

    OStringStream str;

    if (isNotRegion || isRegion || isPatch)
    {
        str << "Field " << field.internalField().name() << " of type "
            << field.type() << " cannot apply to patch " << patch_.name()
            << " because values are mapped from ";
    }

    if (isNotRegion)
    {
        str << "a different region";
    }
    else if (isRegion)
    {
        str << "within the same region";
    }
    else if (isPatch)
    {
        str << "the same patch";
    }

    if (isNotRegion || isRegion || isPatch)
    {
        FatalIOErrorInFunction(context)
            << stringOps::breakIntoIndentedLines(str.str()).c_str()
            << exit(FatalIOError);
    }
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::mappedPatchBase::distribute(const Field<Type>& fld) const
{
    if (sameUntransformedPatch())
    {
        return fld;
    }

    if (!patchToPatchIsUsed_)
    {
        if (mapPtr_.empty())
        {
            calcMapping();
        }

        tmp<Field<Type>> tResult(new Field<Type>(fld, nbrPatchFaceIndices_));
        mapPtr_->distribute(tResult.ref());
        return transform_.transform().transform(tResult);
    }
    else
    {
        if
        (
            !patchToPatchIsValid_
         && !(symmetric() && nbrMappedPatch().patchToPatchIsValid_)
        )
        {
            calcMapping();
        }

        if (!patchToPatchIsValid_ && !symmetric())
        {
            calcMapping();
        }

        return
            transform_.transform().transform
            (
                patchToPatchIsValid_
              ? patchToPatchPtr_->tgtToSrc(fld)
              : nbrMappedPatch().patchToPatchPtr_->srcToTgt(fld)
            );
    }
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::mappedPatchBase::distribute(const tmp<Field<Type>>& fld) const
{
    tmp<Field<Type>> tResult = distribute(fld());
    fld.clear();
    return tResult;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::mappedPatchBase::reverseDistribute(const Field<Type>& fld) const
{
    if (sameUntransformedPatch())
    {
        return fld;
    }

    if (!patchToPatchIsUsed_)
    {
        if (mapPtr_.empty())
        {
            calcMapping();
        }

        Field<Type> nbrFld(fld);
        mapPtr_->reverseDistribute(nbrPatchFaceIndices_.size(), nbrFld);
        tmp<Field<Type>> tResult(new Field<Type>(nbrPolyPatch().size()));
        tResult.ref().rmap(nbrFld, nbrPatchFaceIndices_);
        return transform_.transform().invTransform(tResult);
    }
    else
    {
        if
        (
            !patchToPatchIsValid_
         && !(symmetric() && nbrMappedPatch().patchToPatchIsValid_)
        )
        {
            calcMapping();
        }

        return
            transform_.transform().invTransform
            (
                patchToPatchIsValid_
              ? patchToPatchPtr_->srcToTgt(fld)
              : nbrMappedPatch().patchToPatchPtr_->tgtToSrc(fld)
            );
    }
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::mappedPatchBase::reverseDistribute(const tmp<Field<Type>>& fld) const
{
    tmp<Field<Type>> tResult = reverseDistribute(fld());
    fld.clear();
    return tResult;
}


// ************************************************************************* //
