/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025-2026 OpenFOAM Foundation
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

#include "LagrangianFieldSource.H"
#include "LagrangianMesh.H"
#include "LagrangianSubMesh.H"
#include "LagrangianSource.H"
#include "LagrangianInjection.H"
#include "dictionary.H"
#include "dlLibraryTable.H"
#include "regIOobject.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::LagrangianFieldSource<Type>::LagrangianFieldSource
(
    const regIOobject& iIo
)
:
    LagrangianFieldSourceBase(iIo),
    libs_(fileNameList::null()),
    internalField_
    (
        refCastNull<const LagrangianInternalDynamicField<Type>>(iIo)
    ),
    internalNonDynamicField_
    (
        refCastNull<const LagrangianInternalField<Type>>(iIo)
    )
{}


template<class Type>
Foam::LagrangianFieldSource<Type>::LagrangianFieldSource
(
    const regIOobject& iIo,
    const dictionary& dict
)
:
    LagrangianFieldSourceBase(iIo),
    libs_(dict.lookupOrDefault("libs", fileNameList::null())),
    internalField_
    (
        refCastNull<const LagrangianInternalDynamicField<Type>>(iIo)
    ),
    internalNonDynamicField_
    (
        refCastNull<const LagrangianInternalField<Type>>(iIo)
    )
{}


template<class Type>
Foam::LagrangianFieldSource<Type>::LagrangianFieldSource
(
    const LagrangianFieldSource<Type>& field,
    const regIOobject& iIo
)
:
    LagrangianFieldSourceBase(iIo),
    libs_(field.libs_),
    internalField_
    (
        refCastNull<const LagrangianInternalDynamicField<Type>>(iIo)
    ),
    internalNonDynamicField_
    (
        refCastNull<const LagrangianInternalField<Type>>(iIo)
    )
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

template<class Type>
Foam::autoPtr<Foam::LagrangianFieldSource<Type>>
Foam::LagrangianFieldSource<Type>::New
(
    const word& fieldSourceType,
    const regIOobject& iIo
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
            << exit(FatalError);
    }

    return cstrIter()(iIo);
}


template<class Type>
Foam::autoPtr<Foam::LagrangianFieldSource<Type>>
Foam::LagrangianFieldSource<Type>::New
(
    const regIOobject& iIo,
    const dictionary& dict
)
{
    const word fieldSourceType(dict.lookup("type"));

    libs.open(dict, "libs", dictionaryConstructorTablePtr_);

    typename dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(fieldSourceType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        if (!disallowGenericLagrangianFieldSource)
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

    return cstrIter()(iIo, dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::LagrangianFieldSource<Type>::~LagrangianFieldSource()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Type>
const Foam::dimensionSet&
Foam::LagrangianFieldSource<Type>::internalDimensions() const
{
    if (notNull(internalField_))
    {
        return internalField_.dimensions();
    }

    if (notNull(internalNonDynamicField_))
    {
        return internalNonDynamicField_.dimensions();
    }

    FatalErrorInFunction
        << "Dimensions of internal object " << internalIo_.name()
        << " could not be determined"
        << exit(FatalError);

    return NullObjectRef<dimensionSet>();
}


template<class Type>
const Foam::LagrangianInternalDynamicField<Type>&
Foam::LagrangianFieldSource<Type>::internalField() const
{
    if (notNull(internalField_))
    {
        return internalField_;
    }

    FatalErrorInFunction
        << "Internal field " << internalIo_.name() << " is not of type "
        << LagrangianInternalDynamicField<Type>::typeName
        << exit(FatalError);

    return NullObjectRef<LagrangianInternalDynamicField<Type>>();
}


template<class Type>
Foam::tmp<Foam::LagrangianSubField<Type>>
Foam::LagrangianFieldSource<Type>::sourceValue
(
    const LagrangianSource& source,
    const LagrangianSubMesh& subMesh
) const
{
    FatalErrorInFunction
        << "The " << type() << " condition for source " << source.name()
        << " of field " << this->internalName()
        << " does not define a value for a continuous source "
        << exit(FatalError);

    return tmp<LagrangianSubField<Type>>(nullptr);
}


template<class Type>
Foam::tmp<Foam::LagrangianSubScalarField>
Foam::LagrangianFieldSource<Type>::internalCoeff
(
    const LagrangianSource& source,
    const LagrangianSubMesh& subMesh
) const
{
    FatalErrorInFunction
        << "The " << type() << " condition for source " << source.name()
        << " of field " << this->internalName()
        << " does not define a value for a continuous source "
        << exit(FatalError);

    return tmp<LagrangianSubScalarField>(nullptr);
}


template<class Type>
Foam::tmp<Foam::LagrangianSubField<Type>>
Foam::LagrangianFieldSource<Type>::sourceCoeff
(
    const LagrangianSource& source,
    const LagrangianSubMesh& subMesh
) const
{
    return (1 - internalCoeff(source, subMesh))*sourceValue(source, subMesh);
}


template<class Type>
Foam::tmp<Foam::LagrangianSubField<Type>>
Foam::LagrangianFieldSource<Type>::value
(
    const LagrangianSource& source,
    const LagrangianSubMesh& subMesh
) const
{
    return
        sourceCoeff(source, subMesh)
      + internalCoeff(source, subMesh)*subMesh.sub(internalField());
}


template<class Type>
Foam::tmp<Foam::LagrangianSubField<Type>>
Foam::LagrangianFieldSource<Type>::value
(
    const LagrangianInjection& injection,
    const LagrangianSubMesh& subMesh
) const
{
    FatalErrorInFunction
        << "The " << type() << " condition for injection " << injection.name()
        << " of field " << this->internalName()
        << " does not define a value for an instantaneous injection"
        << exit(FatalError);

    return tmp<LagrangianSubField<Type>>(nullptr);
}


template<class Type>
Foam::tmp<Foam::LagrangianSubField<Type>>
Foam::LagrangianFieldSource<Type>::value
(
    const LagrangianSubMesh& subMesh
) const
{
    FatalErrorInFunction
        << "The " << type() << " condition for field "
        << this->internalName() << " does not define a value without a model"
        << exit(FatalError);

    return tmp<LagrangianSubField<Type>>(nullptr);
}


template<class Type>
void Foam::LagrangianFieldSource<Type>::write(Ostream& os) const
{
    writeEntry(os, "type", type());
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class Type>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const LagrangianFieldSource<Type>& ptf
)
{
    ptf.write(os);

    os.check
    (
        "Ostream& operator<<(Ostream&, const LagrangianFieldSource<Type>&)"
    );

    return os;
}


// ************************************************************************* //
