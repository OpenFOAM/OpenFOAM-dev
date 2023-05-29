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

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

template<class Type>
Foam::autoPtr<Foam::pointPatchField<Type>> Foam::pointPatchField<Type>::New
(
    const word& patchFieldType,
    const word& actualPatchType,
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF
)
{
    if (debug)
    {
        InfoInFunction
            << "patchFieldType = " << patchFieldType
            << ", actualPatchType = " << actualPatchType
            << ", patchType = " << p.type()
            << endl;
    }

    typename pointPatchConstructorTable::iterator cstrIter =
        pointPatchConstructorTablePtr_->find(patchFieldType);

    if (cstrIter == pointPatchConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown null-constructable patchField type " << patchFieldType
            << " for patch " << p.name() << " of type " << p.type()
            << " for field " << iF.name() << nl << nl
            << "Valid null-constructable patchField types are :" << endl
            << pointPatchConstructorTablePtr_->sortedToc()
            << abort(FatalError);
    }

    if
    (
        actualPatchType == word::null
     || actualPatchType != p.type()
    )
    {
        typename pointPatchConstructorTable::iterator patchTypeCstrIter =
            pointPatchConstructorTablePtr_->find(p.type());

        if (patchTypeCstrIter != pointPatchConstructorTablePtr_->end())
        {
            return patchTypeCstrIter()(p, iF);
        }
        else
        {
            return cstrIter()(p, iF);
        }
    }
    else
    {
        return cstrIter()(p, iF);
    }
}


template<class Type>
Foam::autoPtr<Foam::pointPatchField<Type>> Foam::pointPatchField<Type>::New
(
    const word& patchFieldType,
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF
)
{
    return New(patchFieldType, word::null, p, iF);
}


template<class Type>
Foam::autoPtr<Foam::pointPatchField<Type>> Foam::pointPatchField<Type>::New
(
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF,
    const dictionary& dict
)
{
    word patchFieldType(dict.lookup("type"));

    if (debug)
    {
        InfoInFunction
            << "patchFieldType = " << patchFieldType
            << ", actualPatchType = "
            << dict.lookupOrDefault<word>("patchType", word::null)
            << ", patchType = " << p.type()
            << endl;
    }

    typename dictionaryConstructorTable::iterator cstrIter
        = dictionaryConstructorTablePtr_->find(patchFieldType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        if (!disallowGenericPointPatchField)
        {
            cstrIter = dictionaryConstructorTablePtr_->find("generic");
        }

        if (cstrIter == dictionaryConstructorTablePtr_->end())
        {
            FatalIOErrorInFunction
            (
                dict
            )   << "Unknown patchField type " << patchFieldType
                << " for patch " << p.name() << " of type " << p.type()
                << " for field " << iF.name() << nl << nl
                << "Valid patchField types are :" << endl
                << dictionaryConstructorTablePtr_->sortedToc()
                << exit(FatalIOError);
        }
    }

    if
    (
       !dict.found("patchType")
     || word(dict.lookup("patchType")) != p.type()
    )
    {
        typename dictionaryConstructorTable::iterator patchTypeCstrIter =
            dictionaryConstructorTablePtr_->find(p.type());

        if
        (
            patchTypeCstrIter != dictionaryConstructorTablePtr_->end()
         && patchTypeCstrIter() != cstrIter()
        )
        {
            FatalIOErrorInFunction
            (
                dict
            )   << "inconsistent patch and patchField types for \n"
                << "    patch type " << p.type()
                << " and patchField type " << patchFieldType
                << " for field " << iF.name()
                << exit(FatalIOError);
        }
    }

    return cstrIter()(p, iF, dict);
}


template<class Type>
Foam::autoPtr<Foam::pointPatchField<Type>> Foam::pointPatchField<Type>::New
(
    const pointPatchField<Type>& ptf,
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF,
    const pointPatchFieldMapper& pfMapper
)
{
    if (debug)
    {
        InfoInFunction << "Constructing pointPatchField<Type>" << endl;
    }

    typename patchMapperConstructorTable::iterator cstrIter =
        patchMapperConstructorTablePtr_->find(ptf.type());

    if (cstrIter == patchMapperConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown patchField type " << ptf.type()
            << " for patch " << p.name() << " of type " << p.type()
            << " for field " << iF.name() << nl << nl
            << "Valid patchField types are :" << endl
            << patchMapperConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return cstrIter()(ptf, p, iF, pfMapper);
}


// ************************************************************************* //
