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

#define DEFINE_FIELD_MAPPER_MAP_OPERATOR(Type, Modifier)                       \
                                                                               \
    virtual void operator()                                                    \
    (                                                                          \
        Field<Type>& f,                                                        \
        const Field<Type>& mapF                                                \
    ) const Modifier;                                                          \
                                                                               \
    virtual tmp<Field<Type>> operator()                                        \
    (                                                                          \
        const Field<Type>& mapF                                                \
    ) const Modifier;


#define DEFINE_FIELD_MAPPER_MAP_OR_ASSIGN_OPERATOR(Type, Modifier)             \
                                                                               \
    virtual void operator()                                                    \
    (                                                                          \
        Field<Type>& f,                                                        \
        const Field<Type>& mapF,                                               \
        const Type& unmappedVal                                                \
    ) const Modifier;                                                          \
                                                                               \
    virtual tmp<Field<Type>> operator()                                        \
    (                                                                          \
        const Field<Type>& mapF,                                               \
        const Type& unmappedVal                                                \
    ) const Modifier;                                                          \
                                                                               \
    virtual void operator()                                                    \
    (                                                                          \
        Field<Type>& f,                                                        \
        const Field<Type>& mapF,                                               \
        const FieldFunctor<Type>& unmappedFunc                                 \
    ) const Modifier;                                                          \
                                                                               \
    virtual tmp<Field<Type>> operator()                                        \
    (                                                                          \
        const Field<Type>& mapF,                                               \
        const FieldFunctor<Type>& unmappedFunc                                 \
    ) const Modifier;


#define IMPLEMENT_FIELD_MAPPER_MAP_OPERATOR(Type, mapperType)                  \
                                                                               \
    void Foam::mapperType::operator()                                          \
    (                                                                          \
        Field<Type>& f,                                                        \
        const Field<Type>& mapF                                                \
    ) const                                                                    \
    {                                                                          \
        map(f, mapF);                                                          \
    }                                                                          \
                                                                               \
    Foam::tmp<Foam::Field<Foam::Type>> Foam::mapperType::operator()            \
    (                                                                          \
        const Field<Type>& mapF                                                \
    ) const                                                                    \
    {                                                                          \
        return map(mapF);                                                      \
    }


#define IMPLEMENT_FIELD_MAPPER_MAP_OR_ASSIGN_OPERATOR(Type, mapperType)        \
                                                                               \
    void Foam::mapperType::operator()                                          \
    (                                                                          \
        Field<Type>& f,                                                        \
        const Field<Type>& mapF,                                               \
        const Type& unmappedVal                                                \
    ) const                                                                    \
    {                                                                          \
        mapOrAssign(f, mapF, unmappedVal);                                     \
    }                                                                          \
                                                                               \
    Foam::tmp<Foam::Field<Foam::Type>> Foam::mapperType::operator()            \
    (                                                                          \
        const Field<Type>& mapF,                                               \
        const Type& unmappedVal                                                \
    ) const                                                                    \
    {                                                                          \
        return mapOrAssign(mapF, unmappedVal);                                 \
    }                                                                          \
                                                                               \
    void Foam::mapperType::operator()                                          \
    (                                                                          \
        Field<Type>& f,                                                        \
        const Field<Type>& mapF,                                               \
        const FieldFunctor<Type>& unmappedFunc                                 \
    ) const                                                                    \
    {                                                                          \
        mapOrAssign(f, mapF, unmappedFunc);                                    \
    }                                                                          \
                                                                               \
    Foam::tmp<Foam::Field<Foam::Type>> Foam::mapperType::operator()            \
    (                                                                          \
        const Field<Type>& mapF,                                               \
        const FieldFunctor<Type>& unmappedFunc                                 \
    ) const                                                                    \
    {                                                                          \
        return mapOrAssign(mapF, unmappedFunc);                                \
    }


// ************************************************************************* //
