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

#include "csvSetWriter.H"
#include "coordSet.H"
#include "fileName.H"
#include "OFstream.H"
#include "OSspecific.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(csvSetWriter, 0);
    addToRunTimeSelectionTable(setWriter, csvSetWriter, word);
    addToRunTimeSelectionTable(setWriter, csvSetWriter, dict);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::csvSetWriter::writeValueSeparator(Ostream& os) const
{
    os << separator_;
}


void Foam::csvSetWriter::writeSegmentSeparator(Ostream& os) const
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::csvSetWriter::csvSetWriter
(
    const IOstream::streamFormat writeFormat,
    const IOstream::compressionType writeCompression
)
:
    setWriter(writeFormat, writeCompression),
    separator_(',')
{}


Foam::csvSetWriter::csvSetWriter(const dictionary& dict)
:
    setWriter(dict),
    separator_(dict.lookupOrDefault<string>("separator", string(","))[0])
{}


Foam::csvSetWriter::csvSetWriter(const csvSetWriter& writer)
:
    setWriter(writer),
    separator_(writer.separator_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::csvSetWriter::~csvSetWriter()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::csvSetWriter::write
(
    const fileName& outputDir,
    const fileName& setName,
    const coordSet& set,
    const wordList& valueSetNames
    #define TypeValueSetsConstArg(Type, nullArg) \
        , const UPtrList<const Field<Type>>& Type##ValueSets
    FOR_ALL_FIELD_TYPES(TypeValueSetsConstArg)
    #undef TypeValueSetsConstArg
) const
{
    if (!isDir(outputDir))
    {
        mkDir(outputDir);
    }

    OFstream os
    (
        outputDir/setName + ".csv",
        IOstream::ASCII,
        IOstream::currentVersion,
        writeCompression_
    );

    writeTableHeader
    (
        set,
        valueSetNames,
        #define TypeValueSetsParameter(Type, nullArg) Type##ValueSets,
        FOR_ALL_FIELD_TYPES(TypeValueSetsParameter)
        #undef TypeValueSetsParameter
        os
    );

    os << nl;

    writeTable
    (
        set,
        #define TypeValueSetsParameter(Type, nullArg) Type##ValueSets,
        FOR_ALL_FIELD_TYPES(TypeValueSetsParameter)
        #undef TypeValueSetsParameter
        os
    );

    os << nl;
}


// ************************************************************************* //
