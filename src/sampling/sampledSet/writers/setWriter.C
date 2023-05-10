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

#include "setWriter.H"
#include "coordSet.H"
#include "OFstream.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(setWriter, 0);
    defineRunTimeSelectionTable(setWriter, word);
    defineRunTimeSelectionTable(setWriter, dict);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::setWriter::writeValueSeparator(Ostream& os) const
{
    os << token::SPACE;
}


void Foam::setWriter::writeCoordSeparator(Ostream& os) const
{
    os << nl;
}


void Foam::setWriter::writeSegmentSeparator(Ostream& os) const
{
    os << nl << nl;
}


const Foam::List<Foam::string>& Foam::setWriter::delimiters() const
{
    if (!delimiters_.valid())
    {
        delimiters_.set(new List<string>(3));
        OStringStream oss;
        writeValueSeparator(oss);
        delimiters_()[0] = oss.str();
        oss.rewind();
        writeCoordSeparator(oss);
        delimiters_()[1] = oss.str();
        oss.rewind();
        writeSegmentSeparator(oss);
        delimiters_()[2] = oss.str();
        oss.rewind();
    }

    return delimiters_;
}


inline Foam::Ostream& Foam::setWriter::writeWord
(
    const word& w,
    Ostream& os,
    const bool align,
    const unsigned long alignPad
) const
{
    string s = w;
    forAll(delimiters(), i)
    {
        if (w.find(delimiters()[i]) != string::npos)
        {
            s = "\"" + w + "\"";
            break;
        }
    }

    if (!align)
    {
        os  << s.c_str();
    }
    else if (s.size() < columnWidth(os) - alignPad)
    {
        os  << string(columnWidth(os) - alignPad - s.size(), ' ').c_str()
            << s.c_str();
    }
    else
    {
        os  << s(columnWidth(os) - alignPad - 3).c_str() << "...";
    }

    return os;
}


void Foam::setWriter::writeTableHeader
(
    const coordSet& set,
    const wordList& valueSetNames,
    #define TypeValueSetsConstArg(Type, nullArg) \
        const UPtrList<const Field<Type>>& Type##ValueSets ,
    FOR_ALL_FIELD_TYPES(TypeValueSetsConstArg)
    #undef TypeValueSetsConstArg
    Ostream& os,
    const bool align,
    const unsigned long alignPad
) const
{
    bool first = true;

    // Write coordinate names
    if (set.hasScalarAxis())
    {
        if (!first) writeValueSeparator(os);

        writeWord(set.scalarName(), os, align, first*alignPad);
        first = false;
    }
    if (set.hasPointAxis())
    {
        for (direction cmpt = 0; cmpt < pTraits<point>::nComponents; ++ cmpt)
        {
            if (!first) writeValueSeparator(os);

            const bool separator =
                !set.pointName().empty()
             && strlen(pTraits<point>::componentNames[cmpt]) > 0;

            const word& w =
                set.pointName()
              + (separator ? "_" : "")
              + pTraits<point>::componentNames[cmpt];

            writeWord(w, os, align, first*alignPad);
            first = false;
        }
    }

    // Write value names
    forAll(scalarValueSets, fieldi)
    {
        #define WriteTypeValueSetNames(Type, nullArg)                       \
            if (Type##ValueSets.set(fieldi))                                \
            {                                                               \
                const label nCmpt = pTraits<Type>::nComponents;             \
                                                                            \
                for (direction cmpt = 0; cmpt < nCmpt; ++ cmpt)             \
                {                                                           \
                    if (!first) writeValueSeparator(os);                    \
                                                                            \
                    const bool separator =                                  \
                        !valueSetNames[fieldi].empty()                      \
                     && strlen(pTraits<Type>::componentNames[cmpt]) > 0;    \
                                                                            \
                    const word w =                                          \
                        valueSetNames[fieldi]                               \
                      + (separator ? "_" : "")                              \
                      + pTraits<Type>::componentNames[cmpt];                \
                                                                            \
                    writeWord(w, os, align, first*alignPad);                \
                    first = false;                                          \
                }                                                           \
            }
        FOR_ALL_FIELD_TYPES(WriteTypeValueSetNames);
        #undef WriteTypeValueSetNames
    }
}


void Foam::setWriter::writeTable
(
    const coordSet& set,
    #define TypeValueSetsConstArg(Type, nullArg) \
        const UPtrList<const Field<Type>>& Type##ValueSets ,
    FOR_ALL_FIELD_TYPES(TypeValueSetsConstArg)
    #undef TypeValueSetsConstArg
    Ostream& os,
    const bool align
) const
{
    forAll(set, pointi)
    {
        if (pointi != 0)
        {
            writeCoordSeparator(os);

            if (set.segments()[pointi] != set.segments()[pointi - 1])
            {
                writeSegmentSeparator(os);
            }
        }

        bool first = true;

        if (set.hasScalarAxis())
        {
            if (!first) writeValueSeparator(os);
            writeValue(set.scalarCoord(pointi), os, align);
            first = false;
        }
        if (set.hasPointAxis())
        {
            if (!first) writeValueSeparator(os);
            writeValue(set.pointCoord(pointi), os, align);
            first = false;
        }

        forAll(scalarValueSets, fieldi)
        {
            #define WriteTypeValueSets(Type, nullArg)                       \
                if (Type##ValueSets.set(fieldi))                            \
                {                                                           \
                    if (!first) writeValueSeparator(os);                    \
                    writeValue(Type##ValueSets[fieldi][pointi], os, align); \
                    first = false;                                          \
                }
            FOR_ALL_FIELD_TYPES(WriteTypeValueSets);
            #undef WriteTypeValueSets
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::setWriter::setWriter
(
    const IOstream::streamFormat writeFormat,
    const IOstream::compressionType writeCompression
)
:
    writeFormat_(writeFormat),
    writeCompression_(writeCompression)
{}


Foam::setWriter::setWriter(const dictionary& dict)
:
    writeFormat_
    (
        dict.found("writeFormat")
      ? IOstream::formatEnum(dict.lookup("writeFormat"))
      : IOstream::ASCII
    ),
    writeCompression_
    (
        dict.found("writeCompression")
      ? IOstream::compressionEnum(dict.lookup("writeCompression"))
      : IOstream::UNCOMPRESSED
    )
{}


Foam::setWriter::setWriter(const setWriter& writer)
:
    writeFormat_(writer.writeFormat_),
    writeCompression_(writer.writeCompression_)
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::setWriter> Foam::setWriter::New
(
    const word& writeType,
    const IOstream::streamFormat writeFormat,
    const IOstream::compressionType writeCompression
)
{
    typename wordConstructorTable::iterator cstrIter =
        wordConstructorTablePtr_->find(writeType);

    if (cstrIter == wordConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown write type "
            << writeType << nl << nl
            << "Valid write types : " << endl
            << wordConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<setWriter>(cstrIter()(writeFormat, writeCompression));
}


Foam::autoPtr<Foam::setWriter> Foam::setWriter::New
(
    const word& writeType,
    const dictionary& dict
)
{
    // find constructors with dictionary options
    dictConstructorTable::iterator cstrIter =
        dictConstructorTablePtr_->find(writeType);

    if (cstrIter == dictConstructorTablePtr_->end())
    {
        const IOstream::streamFormat writeFormat =
            dict.found("writeFormat")
          ? IOstream::formatEnum(dict.lookup("writeFormat"))
          : IOstream::ASCII;

        const IOstream::compressionType writeCompression =
            dict.found("writeCompression")
          ? IOstream::compressionEnum(dict.lookup("writeCompression"))
          : IOstream::UNCOMPRESSED;

        // Revert to versions without options
        return
            Foam::setWriter::New
            (
                writeType,
                writeFormat,
                writeCompression
            );
    }

    return autoPtr<setWriter>(cstrIter()(dict));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::setWriter::~setWriter()
{}


// ************************************************************************* //
