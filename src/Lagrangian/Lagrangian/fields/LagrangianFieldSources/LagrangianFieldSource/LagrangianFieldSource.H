/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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
    Foam::LagrangianFieldSource

Description
    Base class for Lagrangian source conditions

SourceFiles
    LagrangianFieldSource.C

\*---------------------------------------------------------------------------*/

#ifndef LagrangianFieldSource_H
#define LagrangianFieldSource_H

#include "dlLibraryTable.H"
#include "LagrangianFieldsFwd.H"
#include "LagrangianSubFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class regIOobject;
class dimensionSet;
class objectRegistry;
class dictionary;
class LagrangianMesh;
class LagrangianSubMesh;
class LagrangianModel;
class LagrangianSource;
class LagrangianInjection;

template<class Type>
class LagrangianFieldSource;

// Forward declaration of friend functions and operators
template<class Type>
Ostream& operator<<(Ostream&, const LagrangianFieldSource<Type>&);

/*---------------------------------------------------------------------------*\
                         Class LagrangianFieldSource Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class LagrangianFieldSource
{
    // Private Data

        //- Optional list of libraries required for this field source value
        fileNameList libs_;

        //- Reference to the internal IO object
        const regIOobject& internalIo_;

        //- Reference to the internal field
        const LagrangianInternalDynamicField<Type>& internalField_;

        //- Reference to the non-dynamic internal field
        const LagrangianInternalField<Type>& internalNonDynamicField_;


public:

    //- Runtime type information
    TypeName("LagrangianFieldSource");


    //- Debug switch to disallow the use of genericLagrangianFieldSource
    static int disallowGenericLagrangianFieldSource;


    // Declare run-time constructor selection tables

        //- Select given internal field
        declareRunTimeSelectionTable
        (
            autoPtr,
            LagrangianFieldSource,
            null,
            (
                const regIOobject& iIo
            ),
            (iIo)
        );

        //- Select given internal field and dictionary
        declareRunTimeSelectionTable
        (
            autoPtr,
            LagrangianFieldSource,
            dictionary,
            (
                const regIOobject& iIo,
                const dictionary& dict
            ),
            (iIo, dict)
        );


    // Constructors

        //- Construct from internal field
        LagrangianFieldSource(const regIOobject&);

        //- Construct from internal field and dictionary
        LagrangianFieldSource(const regIOobject&, const dictionary&);

        //- Disallow copy without setting the internal field reference
        LagrangianFieldSource(const LagrangianFieldSource<Type>&) = delete;

        //- Disallow clone without setting the internal field reference
        autoPtr<LagrangianFieldSource<Type>> clone() const
        {
            NotImplemented;
            return autoPtr<LagrangianFieldSource<Type>>(nullptr);
        }

        //- Construct as copy setting the internal field reference
        LagrangianFieldSource
        (
            const LagrangianFieldSource<Type>&,
            const regIOobject&
        );

        //- Construct and return a clone setting the internal field reference
        virtual autoPtr<LagrangianFieldSource<Type>> clone
        (
            const regIOobject&
        ) const = 0;


    // Selectors

        //- Return a pointer to a new field source
        static autoPtr<LagrangianFieldSource<Type>> New
        (
            const word& fieldSourceType,
            const regIOobject&
        );

        //- Return a pointer to a new field source created from a dictionary
        static autoPtr<LagrangianFieldSource<Type>> New
        (
            const regIOobject&,
            const dictionary&
        );


    //- Destructor
    virtual ~LagrangianFieldSource();


    // Member Functions

        // Access

            //- Return local objectRegistry
            const objectRegistry& db() const;

            //- Return internal dimensions reference
            const dimensionSet& internalDimensions() const;

            //- Return internal field reference
            const LagrangianInternalDynamicField<Type>& internalField() const;


        // Evaluation

            //- Return the source value
            virtual tmp<LagrangianSubField<Type>> sourceValue
            (
                const LagrangianSource&,
                const LagrangianSubMesh&
            ) const;

            //- Return the internal coefficient
            virtual tmp<LagrangianSubScalarField> internalCoeff
            (
                const LagrangianSource&,
                const LagrangianSubMesh&
            ) const;

            //- Return the source coefficient
            tmp<LagrangianSubField<Type>> sourceCoeff
            (
                const LagrangianSource&,
                const LagrangianSubMesh&
            ) const;

            //- Return the value for a continuous source
            tmp<LagrangianSubField<Type>> value
            (
                const LagrangianSource&,
                const LagrangianSubMesh&
            ) const;

            //- Return the value for an instantaneous injection
            virtual tmp<LagrangianSubField<Type>> value
            (
                const LagrangianInjection&,
                const LagrangianSubMesh&
            ) const;

            //- Return the value for either a continuous source or an
            //  instantaneous injection
            tmp<LagrangianSubField<Type>> value
            (
                const LagrangianModel&,
                const LagrangianSubMesh&
            ) const;

            //- Lookup and return another field source
            template<class OtherType>
            const LagrangianFieldSource<OtherType>& fieldSource
            (
                const word& name,
                const LagrangianModel&
            ) const;

            //- Lookup and return another field source and then cast to a
            //  specified type. Handle errors.
            template<class OtherType, class OtherFieldSourceType>
            const OtherFieldSourceType& fieldSourceCast
            (
                const word& name,
                const LagrangianModel&
            ) const;

            //- Cast a model to the given type. Handle errors.
            template<class OtherModelType>
            const OtherModelType& modelCast
            (
                const LagrangianModel&
            ) const;

            //- Lookup and return the value of another field source
            template<class OtherType>
            tmp<LagrangianSubField<OtherType>> value
            (
                const word& name,
                const LagrangianSource&,
                const LagrangianSubMesh&
            ) const;

            //- Lookup and return the value of another field source
            template<class OtherType>
            tmp<LagrangianSubField<OtherType>> value
            (
                const word& name,
                const LagrangianInjection&,
                const LagrangianSubMesh&
            ) const;


        // I-O

            //- Write
            virtual void write(Ostream&) const;


    // Member Operators

        //- Ostream operator
        friend Ostream& operator<< <Type>
        (
            Ostream&,
            const LagrangianFieldSource<Type>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "LagrangianFieldSource.C"
#endif


#define addToLagrangianFieldSourceRunTimeSelection(                            \
    TypeFieldSource, typeTypeFieldSource)                                      \
    addToRunTimeSelectionTable                                                 \
    (                                                                          \
        TypeFieldSource,                                                       \
        typeTypeFieldSource,                                                   \
        dictionary                                                             \
    )


#define addNullConstructableToLagrangianFieldSourceRunTimeSelection(           \
    TypeFieldSource, typeTypeFieldSource)                                      \
    addToRunTimeSelectionTable                                                 \
    (                                                                          \
        TypeFieldSource,                                                       \
        typeTypeFieldSource,                                                   \
        null                                                                   \
    );                                                                         \
    addToLagrangianFieldSourceRunTimeSelection                                 \
    (                                                                          \
        TypeFieldSource,                                                       \
        typeTypeFieldSource                                                    \
    )


#define makeLagrangianTypeFieldSource(                                         \
    TypeFieldSource, typeTypeFieldSource)                                      \
    defineTypeNameAndDebug(typeTypeFieldSource, 0);                            \
    addToLagrangianFieldSourceRunTimeSelection                                 \
    (                                                                          \
        TypeFieldSource,                                                       \
        typeTypeFieldSource                                                    \
    )


#define makeNullConstructableLagrangianTypeFieldSource(                        \
    TypeFieldSource, typeTypeFieldSource)                                      \
    defineTypeNameAndDebug(typeTypeFieldSource, 0);                            \
    addNullConstructableToLagrangianFieldSourceRunTimeSelection                \
    (                                                                          \
        TypeFieldSource,                                                       \
        typeTypeFieldSource                                                    \
    )


#define makeTemplateLagrangianTypeFieldSource(fieldType, type)                 \
    defineNamedTemplateTypeNameAndDebug                                        \
    (                                                                          \
        CAT4(type, Lagrangian, CAPITALIZE(fieldType), FieldSource),            \
        0                                                                      \
    );                                                                         \
    addToLagrangianFieldSourceRunTimeSelection                                 \
    (                                                                          \
        CAT3(Lagrangian, CAPITALIZE(fieldType), FieldSource),                  \
        CAT4(type, Lagrangian, CAPITALIZE(fieldType), FieldSource)             \
    );


#define makeNullConstructableTemplateLagrangianTypeFieldSource(fieldType, type)\
    defineNamedTemplateTypeNameAndDebug                                        \
    (                                                                          \
        CAT4(type, Lagrangian, CAPITALIZE(fieldType), FieldSource),            \
        0                                                                      \
    );                                                                         \
    addNullConstructableToLagrangianFieldSourceRunTimeSelection                \
    (                                                                          \
        CAT3(Lagrangian, CAPITALIZE(fieldType), FieldSource),                  \
        CAT4(type, Lagrangian, CAPITALIZE(fieldType), FieldSource)             \
    );


#define makeLagrangianFieldSources(type)                                       \
    FOR_ALL_FIELD_TYPES                                                        \
    (                                                                          \
        makeTemplateLagrangianTypeFieldSource,                                 \
        type                                                                   \
    )


#define makeNullConstructableLagrangianFieldSources(type)                      \
    FOR_ALL_FIELD_TYPES                                                        \
    (                                                                          \
        makeNullConstructableTemplateLagrangianTypeFieldSource,                \
        type                                                                   \
    )


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
