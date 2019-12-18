/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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

#include "DimensionedField.H"
#include "IOstreams.H"


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type, class GeoMesh>
void Foam::DimensionedField<Type, GeoMesh>::readField
(
    const dictionary& fieldDict,
    const word& fieldDictEntry
)
{
    dimensions_.reset(dimensionSet(fieldDict.lookup("dimensions")));

    Field<Type> f(fieldDictEntry, fieldDict, GeoMesh::size(mesh_));
    this->transfer(f);
}


template<class Type, class GeoMesh>
void Foam::DimensionedField<Type, GeoMesh>::readIfPresent
(
    const word& fieldDictEntry
)
{
    if
    (
        (this->readOpt() == IOobject::READ_IF_PRESENT && this->headerOk())
     || this->readOpt() == IOobject::MUST_READ
     || this->readOpt() == IOobject::MUST_READ_IF_MODIFIED
    )
    {
        readField(dictionary(readStream(typeName)), fieldDictEntry);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type, class GeoMesh>
Foam::DimensionedField<Type, GeoMesh>::DimensionedField
(
    const IOobject& io,
    const Mesh& mesh,
    const word& fieldDictEntry
)
:
    regIOobject(io),
    Field<Type>(0),
    mesh_(mesh),
    dimensions_(dimless)
{
    readField(dictionary(readStream(typeName)), fieldDictEntry);
}


template<class Type, class GeoMesh>
Foam::DimensionedField<Type, GeoMesh>::DimensionedField
(
    const IOobject& io,
    const Mesh& mesh,
    const dictionary& fieldDict,
    const word& fieldDictEntry
)
:
    regIOobject(io),
    Field<Type>(0),
    mesh_(mesh),
    dimensions_(dimless)
{
    readField(fieldDict, fieldDictEntry);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type, class GeoMesh>
bool Foam::DimensionedField<Type, GeoMesh>::writeData
(
    Ostream& os,
    const word& fieldDictEntry
) const
{
    writeEntry(os, "dimensions", dimensions());
    os << nl;

    writeEntry(os, fieldDictEntry, static_cast<const Field<Type>&>(*this));

    // Check state of Ostream
    os.check
    (
        "bool DimensionedField<Type, GeoMesh>::writeData"
        "(Ostream& os, const word& fieldDictEntry) const"
    );

    return (os.good());
}


template<class Type, class GeoMesh>
bool Foam::DimensionedField<Type, GeoMesh>::writeData(Ostream& os) const
{
    return writeData(os, "value");
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class Type, class GeoMesh>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const DimensionedField<Type, GeoMesh>& df
)
{
    df.writeData(os);

    return os;
}


template<class Type, class GeoMesh>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const tmp<DimensionedField<Type, GeoMesh>>& tdf
)
{
    tdf().writeData(os);
    tdf.clear();

    return os;
}


// ************************************************************************* //
