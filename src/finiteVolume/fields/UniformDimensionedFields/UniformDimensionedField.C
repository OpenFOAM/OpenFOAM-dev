/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2026 OpenFOAM Foundation
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

#include "UniformDimensionedField.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class Type>
Foam::UniformDimensionedField<Type>::UniformDimensionedField
(
    const IOobject& io,
    const dimensionSet& dims,
    const bool
)
:
    regIOobject(io),
    dimensioned<Type>(regIOobject::name(), dims, Zero),
    OldTimeField<UniformDimensionedField>(this->time().timeIndex())
{}


template<class Type>
Foam::UniformDimensionedField<Type>::UniformDimensionedField
(
    const IOobject& io,
    const dimensioned<Type>& defaultDt,
    const bool
)
:
    regIOobject(io),
    dimensioned<Type>(defaultDt),
    OldTimeField<UniformDimensionedField>(this->time().timeIndex())
{}


template<class Type>
void Foam::UniformDimensionedField<Type>::read
(
    const IOobject::readOption& ro,
    const dimensionSet& dims
)
{
    if
    (
        ro == IOobject::MUST_READ
     || ro == IOobject::MUST_READ_IF_MODIFIED
     || (ro == IOobject::READ_IF_PRESENT && headerOk())
    )
    {
        const dictionary dict(readStream(type()));

        this->dimensions().read(dict.lookup("dimensions"));
        this->value() = dict.lookup<Type>("value", this->dimensions());

        if (this->dimensions() != dims)
        {
            FatalIOErrorInFunction(dict)
                << "Inconsistent dimensions specified for " << this->name()
                << " " << this->dimensions()
                << ", required dimensions are " << dims
                << exit(FatalIOError);
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::UniformDimensionedField<Type>::UniformDimensionedField
(
    const IOobject& io
)
:
    UniformDimensionedField(io, dimensions::invalid, false)
{
    const dictionary dict(readStream(type()));
    this->dimensions().read(dict.lookup("dimensions"));
    this->value() = dict.lookup<Type>("value", this->dimensions());
}


template<class Type>
Foam::UniformDimensionedField<Type>::UniformDimensionedField
(
    const IOobject& io,
    const dimensionSet& dims
)
:
    UniformDimensionedField(io, dims, false)
{
    this->read(IOobject::MUST_READ, dims);
}


template<class Type>
Foam::UniformDimensionedField<Type>::UniformDimensionedField
(
    const IOobject& io,
    const dimensioned<Type>& defaultDt
)
:
    UniformDimensionedField(io, defaultDt, false)
{
    this->read(io.readOpt(), defaultDt.dimensions());
}


template<class Type>
Foam::UniformDimensionedField<Type>::UniformDimensionedField
(
    const UniformDimensionedField<Type>& udt
)
:
    regIOobject(udt),
    dimensioned<Type>(udt),
    OldTimeField<UniformDimensionedField>(udt)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::UniformDimensionedField<Type>::~UniformDimensionedField()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Type& Foam::UniformDimensionedField<Type>::value()
{
    this->storeOldTimes();
    return dimensioned<Type>::value();
}


template<class Type>
const Type& Foam::UniformDimensionedField<Type>::value() const
{
    return dimensioned<Type>::value();
}


template<class Type>
void Foam::UniformDimensionedField<Type>::reset
(
    const UniformDimensionedField<Type>& rhs
)
{
    dimensioned<Type>::operator=(rhs);
}


template<class Type>
bool Foam::UniformDimensionedField<Type>::writeData(Ostream& os) const
{
    writeKeyword(os, "dimensions");
    this->dimensions().write(os) << token::END_STATEMENT << nl;
    writeEntry(os, "value", this->value());
    os << nl;

    return (os.good());
}


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

template<class Type>
void Foam::UniformDimensionedField<Type>::operator==
(
    const UniformDimensionedField<Type>& rhs
)
{
    dimensioned<Type>::operator=(rhs);
}


template<class Type>
void Foam::UniformDimensionedField<Type>::operator=
(
    const UniformDimensionedField<Type>& rhs
)
{
    dimensioned<Type>::operator=(rhs);
}


template<class Type>
void Foam::UniformDimensionedField<Type>::operator=
(
    const dimensioned<Type>& rhs
)
{
    dimensioned<Type>::operator=(rhs);
}


template<class Type>
const Type& Foam::UniformDimensionedField<Type>::operator[](const label) const
{
    return this->value();
}


// ************************************************************************* //
