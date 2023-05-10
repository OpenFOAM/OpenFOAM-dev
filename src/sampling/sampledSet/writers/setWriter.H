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

Class
    Foam::setWriter

Description
    Base class for writing coordinate sets with data

    Example:
    \verbatim
        setWriter::New("vtk")->write
        (
            "myDirectory",
            "mySet",
            coordSet(true, "position", points),
            "p", p,
            "U", U
        );
    \endverbatim

SourceFiles
    setWriter.C

\*---------------------------------------------------------------------------*/

#ifndef setWriter_H
#define setWriter_H

#include "runTimeSelectionTables.H"
#include "autoPtr.H"
#include "fieldTypes.H"
#include "coordSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class setWriter Declaration
\*---------------------------------------------------------------------------*/

class setWriter
{
public:

    // Public Static Functions

        //- Helper for variadic write
        template<class Types, class Type>
        static inline void appendTypeValueSet
        (
            UPtrList<const Field<Types>>& valueSets,
            const Field<Type>& valueSet
        )
        {
            valueSets.resize(valueSets.size() + 1);
            valueSets.set(valueSets.size() - 1, nullptr);
        }

        //- Helper for variadic write
        template<class Type>
        static inline void appendTypeValueSet
        (
            UPtrList<const Field<Type>>& valueSets,
            const Field<Type>& valueSet
        )
        {
            valueSets.resize(valueSets.size() + 1);
            valueSets.set(valueSets.size() - 1, &valueSet);
        }

        //- Helper for variadic write
        static inline void unpackTypeValueSets
        (
            wordList& valueSetNames
            #define TypeValueSetsNonConstArg(Type, nullArg) \
                , UPtrList<const Field<Type>>& Type##ValueSets
            FOR_ALL_FIELD_TYPES(TypeValueSetsNonConstArg)
            #undef TypeValueSetsNonConstArg
        )
        {}

        //- Helper for variadic write
        template<class Type, class ... Args>
        static inline void unpackTypeValueSets
        (
            wordList& valueSetNames
            #define TypeValueSetsNonConstArg(Type, nullArg) \
                , UPtrList<const Field<Type>>& Type##ValueSets
            FOR_ALL_FIELD_TYPES(TypeValueSetsNonConstArg),
            #undef TypeValueSetsNonConstArg
            const word& valueSetName,
            const Field<Type>& valueSet,
            Args& ... args
        )
        {
            valueSetNames.append(valueSetName);
            #define AppendTypeValueSet(Type, nullArg) \
                appendTypeValueSet(Type##ValueSets, valueSet);
            FOR_ALL_FIELD_TYPES(AppendTypeValueSet);
            #undef AppendTypeValueSet

            unpackTypeValueSets
            (
                valueSetNames
                #define TypeValueSetsParameter(Type, nullArg) \
                    , Type##ValueSets
                FOR_ALL_FIELD_TYPES(TypeValueSetsParameter),
                #undef TypeValueSetsParameter
                args ...
            );
        }

        //- Helper for variadic write
        template<class Type, class ... Args>
        static inline void unpackTypeValueSets
        (
            wordList& valueSetNames
            #define TypeValueSetsNonConstArg(Type, nullArg) \
                , UPtrList<const Field<Type>>& Type##ValueSets
            FOR_ALL_FIELD_TYPES(TypeValueSetsNonConstArg),
            #undef TypeValueSetsNonConstArg
            const wordList& valueSetNamesPart,
            const UPtrList<const Field<Type>>& valueSetsPart,
            Args& ... args
        )
        {
            forAll(valueSetNamesPart, i)
            {
                valueSetNames.append(valueSetNamesPart[i]);
                #define AppendTypeValueSet(Type, nullArg) \
                    appendTypeValueSet(Type##ValueSets, valueSetsPart[i]);
                FOR_ALL_FIELD_TYPES(AppendTypeValueSet);
                #undef AppendTypeValueSet
            }

            unpackTypeValueSets
            (
                valueSetNames
                #define TypeValueSetsParameter(Type, nullArg) \
                    , Type##ValueSets
                FOR_ALL_FIELD_TYPES(TypeValueSetsParameter),
                #undef TypeValueSetsParameter
                args ...
            );
        }

        //- Helper for variadic write
        template<class Type, class ... Args>
        static inline void unpackTypeValueSets
        (
            wordList& valueSetNames
            #define TypeValueSetsNonConstArg(Type, nullArg) \
                , UPtrList<const Field<Type>>& Type##ValueSets
            FOR_ALL_FIELD_TYPES(TypeValueSetsNonConstArg),
            #undef TypeValueSetsNonConstArg
            const wordList& valueSetNamesPart,
            const PtrList<Field<Type>>& valueSetsPart,
            Args& ... args
        )
        {
            unpackTypeValueSets
            (
                valueSetNames
                #define TypeValueSetsParameter(Type, nullArg) \
                    , Type##ValueSets
                FOR_ALL_FIELD_TYPES(TypeValueSetsParameter),
                #undef TypeValueSetsParameter
                valueSetNamesPart,
                reinterpret_cast<const UPtrList<const Field<Type>>&>
                (valueSetsPart),
                args ...
            );
        }


protected:

    // Protected Data

        //- Write format
        IOstream::streamFormat writeFormat_;

        //- Write compression
        IOstream::compressionType writeCompression_;

        //- Possible delimiters that need quoting. Constructed on demand from
        //  the write.*Separator functions.
        mutable autoPtr<List<string>> delimiters_;


    // Protected Member Functions

        //- Write a value separator
        virtual void writeValueSeparator(Ostream& os) const;

        //- Write a coordinate separator
        virtual void writeCoordSeparator(Ostream& os) const;

        //- Write a segment separator
        virtual void writeSegmentSeparator(Ostream& os) const;

        //- Possible delimiters that need quoting
        const List<string>& delimiters() const;

        //- Width of columns in tabulated output
        static inline unsigned long columnWidth(Ostream& os)
        {
            return os.precision() + 7;
        }

        //- Write a word
        inline Ostream& writeWord
        (
            const word& w,
            Ostream& os,
            const bool align = false,
            const unsigned long alignPad = 0
        ) const;

        //- Write a value
        template<class Type>
        Ostream& writeValue
        (
            const Type& value,
            Ostream& os,
            const bool align = false
        ) const;

        //- Write multi-column table header
        void writeTableHeader
        (
            const coordSet& set,
            const wordList& valueSetNames,
            #define TypeValueSetsConstArg(Type, nullArg) \
                const UPtrList<const Field<Type>>& Type##ValueSets ,
            FOR_ALL_FIELD_TYPES(TypeValueSetsConstArg)
            #undef TypeValueSetsConstArg
            Ostream& os,
            const bool align = false,
            const unsigned long alignPad = 0
        ) const;

        //- Write multi-column table of data
        void writeTable
        (
            const coordSet& set,
            #define TypeValueSetsConstArg(Type, nullArg) \
                const UPtrList<const Field<Type>>& Type##ValueSets ,
            FOR_ALL_FIELD_TYPES(TypeValueSetsConstArg)
            #undef TypeValueSetsConstArg
            Ostream& os,
            const bool align = false
        ) const;


public:

    //- Runtime type information
    TypeName("setWriter");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            setWriter,
            word,
            (
                const IOstream::streamFormat writeFormat,
                const IOstream::compressionType writeCompression
            ),
            (writeFormat, writeCompression)
        );

        declareRunTimeSelectionTable
        (
            autoPtr,
            setWriter,
            dict,
            (
                const dictionary& dict
            ),
            (dict)
        );


    // Selectors

        //- Select given write options
        static autoPtr<setWriter> New
        (
            const word& writeType,
            const IOstream::streamFormat writeFormat =
                IOstream::ASCII,
            const IOstream::compressionType writeCompression =
                IOstream::UNCOMPRESSED
        );

        //- Select given a dictionary
        static autoPtr<setWriter> New
        (
            const word& writeType,
            const dictionary& dict
        );


    // Constructors

        //- Construct given write options
        setWriter
        (
            const IOstream::streamFormat writeFormat,
            const IOstream::compressionType writeCompression
        );

        //- Construct from dictionary
        setWriter(const dictionary& dict);

        //- Construct copy
        setWriter(const setWriter& writer);

        //- Construct and return a clone
        virtual autoPtr<setWriter> clone() const = 0;


    //- Destructor
    virtual ~setWriter() = 0;


    // Member Functions

        //- Write a coordSet and associated data
        virtual void write
        (
            const fileName& outputDir,
            const fileName& setName,
            const coordSet& set,
            const wordList& valueSetNames
            #define TypeValueSetsConstArg(Type, nullArg) \
                , const UPtrList<const Field<Type>>& Type##ValueSets
            FOR_ALL_FIELD_TYPES(TypeValueSetsConstArg)
            #undef TypeValueSetsConstArg
        ) const = 0;

        //- Write a coordSet and associated data
        virtual void write
        (
            const fileName& outputDir,
            const fileName& setName,
            const coordSet& set,
            const wordList& valueSetNames
            #define TypeValueSetsConstArg(Type, nullArg) \
                , const PtrList<Field<Type>>& Type##ValueSets
            FOR_ALL_FIELD_TYPES(TypeValueSetsConstArg)
            #undef TypeValueSetsConstArg
        ) const
        {
            write
            (
                outputDir,
                setName,
                set,
                valueSetNames
                #define TypeValueSetsParameter(Type, nullArg)              \
                    , reinterpret_cast<const UPtrList<const Field<Type>>&> \
                    (Type##ValueSets)
                FOR_ALL_FIELD_TYPES(TypeValueSetsParameter)
                #undef TypeValueSetsParameter
            );
        }

        //- Write fields for a single surface to file. For use in code where
        //  the fields that are to be written are known. Takes any number of
        //  name, values arguments at the end. E.g.:
        //
        //  write
        //  (
        //      // Output options
        //      "myDirectory", "mySet",
        //
        //      // Geometry
        //      coordSet(true, "position", points),
        //
        //      // Fields
        //      "p", Field<scalar>(points.size(), ...),
        //      "U", Field<vector>(points.size(), ...)
        //  );
        //
        template<class ... Args>
        void inline write
        (
            const fileName& outputDir,
            const fileName& setName,
            const coordSet& set,
            const Args& ... args
        ) const
        {
            wordList valueSetNames;
            #define DeclareTypeValueSets(Type, nullArg) \
                UPtrList<const Field<Type>> Type##ValueSets;
            FOR_ALL_FIELD_TYPES(DeclareTypeValueSets);
            #undef DeclareTypeValueSets

            unpackTypeValueSets
            (
                valueSetNames
                #define TypeValueSetsParameter(Type, nullArg) \
                    , Type##ValueSets
                FOR_ALL_FIELD_TYPES(TypeValueSetsParameter),
                #undef TypeValueSetsParameter
                args ...
            );

            write
            (
                outputDir,
                setName,
                set,
                valueSetNames
                #define TypeValueSetsParameter(Type, nullArg) \
                    , Type##ValueSets
                FOR_ALL_FIELD_TYPES(TypeValueSetsParameter)
                #undef TypeValueSetsParameter
            );
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "setWriterTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
