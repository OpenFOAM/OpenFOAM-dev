/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023-2025 OpenFOAM Foundation
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

#include "fvFieldSource.H"
#include "fvSource.H"
#include "volMesh.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::fvFieldSource<Type>::fvFieldSource
(
    const DimensionedField<Type, volMesh>& iF
)
:
    libs_(fileNameList::null()),
    internalField_(iF)
{}


template<class Type>
Foam::fvFieldSource<Type>::fvFieldSource
(
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    libs_
    (
        dict.lookupOrDefault
        (
            "libs",
            fileNameList::null(),
            dictionary::writeOptionalEntries > 1
        )
    ),
    internalField_(iF)
{}


template<class Type>
Foam::fvFieldSource<Type>::fvFieldSource
(
    const fvFieldSource<Type>& field,
    const DimensionedField<Type, volMesh>& iF
)
:
    libs_(field.libs_),
    internalField_(iF)
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

template<class Type>
Foam::autoPtr<Foam::fvFieldSource<Type>>
Foam::fvFieldSource<Type>::New
(
    const word& fieldSourceType,
    const DimensionedField<Type, volMesh>& iF
)
{
    typename nullConstructorTable::iterator cstrIter =
        nullConstructorTablePtr_->find(fieldSourceType);

    if (cstrIter == nullConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown null-constructable fieldSource type " << fieldSourceType
            << nl << nl << "Valid fieldSource types are :" << endl
            << nullConstructorTablePtr_->sortedToc()
            << exit(FatalIOError);
    }

    return cstrIter()(iF);
}


template<class Type>
Foam::autoPtr<Foam::fvFieldSource<Type>>
Foam::fvFieldSource<Type>::New
(
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
{
    const word fieldSourceType(dict.lookup("type"));

    libs.open(dict, "libs", dictionaryConstructorTablePtr_);

    typename dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(fieldSourceType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        if (!disallowGenericFvFieldSource)
        {
            cstrIter = dictionaryConstructorTablePtr_->find("generic");
        }

        if (cstrIter == dictionaryConstructorTablePtr_->end())
        {
            FatalIOErrorInFunction(dict)
                << "Unknown fieldSource type " << fieldSourceType
                << " for model " << dict.dictName() << nl << nl
                << "Valid fieldSource types are :" << endl
                << dictionaryConstructorTablePtr_->sortedToc()
                << exit(FatalIOError);
        }
    }

    return cstrIter()(iF, dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::fvFieldSource<Type>::~fvFieldSource()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
const Foam::objectRegistry& Foam::fvFieldSource<Type>::db() const
{
    return internalField_.mesh();
}


template<class Type>
const Foam::DimensionedField<Type, Foam::volMesh>&
Foam::fvFieldSource<Type>::internalField() const
{
    return internalField_;
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::fvFieldSource<Type>::sourceValue
(
    const fvSource& model,
    const scalarField& source,
    const labelUList& cells
) const
{
    FatalErrorInFunction
        << "Condition of type " << type() << " cannot be used for source "
        << model.name() << " of field " << internalField().name()
        << " in file " << internalField().objectPath() << " as it is not "
        << "defined for cell sets. It can only be used for sources that "
        << "span the entire mesh." << exit(FatalError);

    return tmp<Field<Type>>(nullptr);
}


template<class Type>
Foam::tmp<Foam::scalarField> Foam::fvFieldSource<Type>::internalCoeff
(
    const fvSource& model,
    const scalarField& source,
    const labelUList& cells
) const
{
    FatalErrorInFunction
        << "Condition of type " << type() << " cannot be used for source "
        << model.name() << " of field " << internalField().name()
        << " in file " << internalField().objectPath() << " as it is not "
        << "defined for cell sets. It can only be used for sources that "
        << "span the entire mesh." << exit(FatalError);

    return tmp<scalarField>(nullptr);
}


template<class Type>
Foam::tmp<Foam::DimensionedField<Type, Foam::volMesh>>
Foam::fvFieldSource<Type>::sourceCoeff
(
    const fvSource& model,
    const DimensionedField<scalar, volMesh>& source
) const
{
    return
        (1 - internalCoeff(model, source))
       *sourceValue(model, source);
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::fvFieldSource<Type>::sourceCoeff
(
    const fvSource& model,
    const scalarField& source,
    const labelUList& cells
) const
{
    return
        (1 - internalCoeff(model, source, cells))
       *sourceValue(model, source, cells);
}


template<class Type>
Foam::tmp<Foam::DimensionedField<Type, Foam::volMesh>>
Foam::fvFieldSource<Type>::value
(
    const fvSource& model,
    const DimensionedField<scalar, volMesh>& source
) const
{
    return
        sourceCoeff(model, source)
      + internalCoeff(model, source)*internalField();
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::fvFieldSource<Type>::value
(
    const fvSource& model,
    const scalarField& source,
    const labelUList& cells
) const
{
    return
        sourceCoeff(model, source, cells)
      + internalCoeff(model, source, cells)*Field<Type>(internalField(), cells);
}


template<class Type>
template<class OtherType>
const Foam::fvFieldSource<OtherType>& Foam::fvFieldSource<Type>::fieldSource
(
    const word& name,
    const fvSource& model
) const
{
    const VolField<OtherType>& vf =
        db().template lookupObject<VolField<OtherType>>(name);

    return vf.sources()[model.name()];
}


template<class Type>
template<class OtherType>
Foam::tmp<Foam::DimensionedField<OtherType, Foam::volMesh>>
Foam::fvFieldSource<Type>::value
(
    const word& name,
    const fvSource& model,
    const DimensionedField<scalar, volMesh>& source
) const
{
    return fieldSource<OtherType>(name, model).value(model, source);
}


template<class Type>
template<class OtherType>
Foam::tmp<Foam::Field<OtherType>> Foam::fvFieldSource<Type>::value
(
    const word& name,
    const fvSource& model,
    const scalarField& source,
    const labelUList& cells
) const
{
    return fieldSource<OtherType>(name, model).value(model, source, cells);
}


template<class Type>
void Foam::fvFieldSource<Type>::write(Ostream& os) const
{
    writeEntry(os, "type", type());

    if (libs_.size())
    {
        writeEntry(os, "libs", libs_);
    }
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class Type>
Foam::Ostream& Foam::operator<<(Ostream& os, const fvFieldSource<Type>& ptf)
{
    ptf.write(os);

    os.check("Ostream& operator<<(Ostream&, const fvFieldSource<Type>&");

    return os;
}


// ************************************************************************* //
