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

Class
    Foam::fvFieldSource

Description
    Base class for finite-volume field sources

SourceFiles
    fvFieldSource.C

\*---------------------------------------------------------------------------*/

#ifndef fvFieldSource_H
#define fvFieldSource_H

#include "DimensionedField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class volMesh;
class dictionary;
class fvSource;

template<class Type>
class fvFieldSource;

// Forward declaration of friend functions and operators
template<class Type>
Ostream& operator<<(Ostream&, const fvFieldSource<Type>&);

/*---------------------------------------------------------------------------*\
                         Class fvFieldSource Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class fvFieldSource
{
    // Private Data

        //- Optional list of libraries required for this field source value
        fileNameList libs_;

        //- Reference to internal field
        const DimensionedField<Type, volMesh>& internalField_;


public:

    //- Runtime type information
    TypeName("fvFieldSource");

    //- Debug switch to disallow the use of genericFvFieldSource
    static int disallowGenericFvFieldSource;


    // Declare run-time constructor selection tables

        //- Select given internal field
        declareRunTimeSelectionTable
        (
            autoPtr,
            fvFieldSource,
            null,
            (const DimensionedField<Type, volMesh>& iF),
            (iF)
        );

        //- Select given internal field and dictionary
        declareRunTimeSelectionTable
        (
            autoPtr,
            fvFieldSource,
            dictionary,
            (
                const DimensionedField<Type, volMesh>& iF,
                const dictionary& dict
            ),
            (iF, dict)
        );


    // Constructors

        //- Construct from internal field
        fvFieldSource(const DimensionedField<Type, volMesh>&);

        //- Construct from internal field and dictionary
        fvFieldSource
        (
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Disallow copy without setting internal field reference
        fvFieldSource(const fvFieldSource<Type>&) = delete;

        //- Disallow clone without setting internal field reference
        autoPtr<fvFieldSource<Type>> clone() const
        {
            NotImplemented;
            return autoPtr<fvFieldSource<Type>>(nullptr);
        }

        //- Copy constructor setting internal field reference
        fvFieldSource
        (
            const fvFieldSource<Type>&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual autoPtr<fvFieldSource<Type>> clone
        (
            const DimensionedField<Type, volMesh>&
        ) const = 0;


    // Selectors

        //- Return a pointer to a new field source
        static autoPtr<fvFieldSource<Type>> New
        (
            const word& fieldSourceType,
            const DimensionedField<Type, volMesh>&
        );

        //- Return a pointer to a new field source created from a dictionary
        static autoPtr<fvFieldSource<Type>> New
        (
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );


    //- Destructor
    virtual ~fvFieldSource();


    // Member Functions

        //- Return the local object registry
        const objectRegistry& db() const;

        //- Return the internal field reference
        const DimensionedField<Type, volMesh>& internalField() const;

        //- Return the source value
        virtual tmp<Field<Type>> sourceValue(const fvSource&) const = 0;

        //- Return the internal coefficient
        virtual tmp<scalarField> internalCoeff(const fvSource&) const = 0;

        //- Return the source coefficient
        virtual tmp<Field<Type>> sourceCoeff(const fvSource&) const;

        //- Return the value
        tmp<Field<Type>> value(const fvSource&) const;

        //- Lookup and return another field source
        template<class OtherType>
        const fvFieldSource<OtherType>& fieldSource
        (
            const word& name,
            const fvSource&
        ) const;

        //- Lookup and return the value of another field source
        template<class OtherType>
        tmp<Field<OtherType>> value
        (
            const word& name,
            const fvSource&
        ) const;

        //- Write
        virtual void write(Ostream&) const;


    //- Ostream operator
    friend Ostream& operator<< <Type>(Ostream&, const fvFieldSource<Type>&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fvFieldSource.C"
#endif


#define makeFvFieldSource(fvTypeFieldSource)                                   \
    defineNamedTemplateTypeNameAndDebug(fvTypeFieldSource, 0);                 \
    template<>                                                                 \
    int fvTypeFieldSource::disallowGenericFvFieldSource                        \
    (                                                                          \
        debug::debugSwitch("disallowGenericFvFieldSource", 0)                  \
    );                                                                         \
    defineTemplateRunTimeSelectionTable(fvTypeFieldSource, null);              \
    defineTemplateRunTimeSelectionTable(fvTypeFieldSource, dictionary)


#define addToFieldSourceRunTimeSelection(TypeFieldSource, typeTypeFieldSource) \
    addToRunTimeSelectionTable                                                 \
    (                                                                          \
        TypeFieldSource,                                                       \
        typeTypeFieldSource,                                                   \
        dictionary                                                             \
    )


#define addNullConstructableToFieldSourceRunTimeSelection(                     \
    TypeFieldSource, typeTypeFieldSource)                                      \
    addToRunTimeSelectionTable                                                 \
    (                                                                          \
        TypeFieldSource,                                                       \
        typeTypeFieldSource,                                                   \
        null                                                                   \
    );                                                                         \
    addToFieldSourceRunTimeSelection(TypeFieldSource, typeTypeFieldSource)


#define makeTypeFieldSource(TypeFieldSource, typeTypeFieldSource)              \
    defineTypeNameAndDebug(typeTypeFieldSource, 0);                            \
    addToFieldSourceRunTimeSelection(TypeFieldSource, typeTypeFieldSource)


#define makeNullConstructableTypeFieldSource(                                  \
    TypeFieldSource, typeTypeFieldSource)                                      \
    defineTypeNameAndDebug(typeTypeFieldSource, 0);                            \
    addNullConstructableToFieldSourceRunTimeSelection                          \
    (                                                                          \
        TypeFieldSource,                                                       \
        typeTypeFieldSource                                                    \
    )


#define makeTemplateTypeFieldSource(fieldType, type)                           \
    defineNamedTemplateTypeNameAndDebug                                        \
    (                                                                          \
        CAT4(type, Fv, CAPITALIZE(fieldType), FieldSource),                    \
        0                                                                      \
    );                                                                         \
    addToFieldSourceRunTimeSelection                                           \
    (                                                                          \
        CAT3(fv, CAPITALIZE(fieldType), FieldSource),                          \
        CAT4(type, Fv, CAPITALIZE(fieldType), FieldSource)                     \
    );


#define makeNullConstructableTemplateTypeFieldSource(fieldType, type)          \
    defineNamedTemplateTypeNameAndDebug                                        \
    (                                                                          \
        CAT4(type, Fv, CAPITALIZE(fieldType), FieldSource),                    \
        0                                                                      \
    );                                                                         \
    addNullConstructableToFieldSourceRunTimeSelection                          \
    (                                                                          \
        CAT3(fv, CAPITALIZE(fieldType), FieldSource),                          \
        CAT4(type, Fv, CAPITALIZE(fieldType), FieldSource)                     \
    );


#define makeFieldSources(type)                                                 \
    FOR_ALL_FIELD_TYPES(makeTemplateTypeFieldSource, type)


#define makeNullConstructableFieldSources(type)                                \
    FOR_ALL_FIELD_TYPES(makeNullConstructableTemplateTypeFieldSource, type)


#define makeFieldSourceTypeName(fieldType, type)                               \
    defineNamedTemplateTypeNameAndDebug                                        \
    (                                                                          \
        CAT4(type, Fv, CAPITALIZE(fieldType), FieldSource),                    \
        0                                                                      \
    );


#define makeFieldSourceTypeNames(type)                                         \
    FOR_ALL_FIELD_TYPES(makeFieldSourceTypeName, type)


#define makeTypeFieldSourceTypedef(fieldType, type)                            \
    typedef type##FvFieldSource<fieldType>                                     \
        CAT4(type, Fv, CAPITALIZE(fieldType), FieldSource);


#define makeTypeFieldSourceTypedefs(type)                                      \
    FOR_ALL_FIELD_TYPES(makeTypeFieldSourceTypedef, type)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
