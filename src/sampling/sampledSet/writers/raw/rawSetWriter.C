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

#include "rawSetWriter.H"
#include "coordSet.H"
#include "OFstream.H"
#include "OSspecific.H"
#include "SubList.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(rawSetWriter, 0);
    addToRunTimeSelectionTable(setWriter, rawSetWriter, word);
    addToRunTimeSelectionTable(setWriter, rawSetWriter, dict);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::rawSetWriter::writeSegmentSeparator(Ostream& os) const
{
    if (separateSegments_)
    {
        os << nl << nl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::rawSetWriter::rawSetWriter
(
    const IOstream::streamFormat writeFormat,
    const IOstream::compressionType writeCompression
)
:
    setWriter(writeFormat, writeCompression),
    separateSegments_(true)
{}


Foam::rawSetWriter::rawSetWriter(const dictionary& dict)
:
    setWriter(dict),
    separateSegments_(true)
{}


Foam::rawSetWriter::rawSetWriter(const rawSetWriter& writer)
:
    setWriter(writer),
    separateSegments_(true)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::rawSetWriter::~rawSetWriter()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::rawSetWriter::write
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
        outputDir/setName + ".xy",
        IOstream::ASCII,
        IOstream::currentVersion,
        writeCompression_
    );

    separateSegments_ = set.segments() != identityMap(set.size());

    os << '#';

    OStringStream oss;
    writeValueSeparator(oss);

    os << oss.str().c_str();

    writeTableHeader
    (
        set,
        valueSetNames,
        #define TypeValueSetsParameter(Type, nullArg) Type##ValueSets,
        FOR_ALL_FIELD_TYPES(TypeValueSetsParameter)
        #undef TypeValueSetsParameter
        os,
        true,
        1 + oss.str().size()
    );

    os << nl;

    writeTable
    (
        set,
        #define TypeValueSetsParameter(Type, nullArg) Type##ValueSets,
        FOR_ALL_FIELD_TYPES(TypeValueSetsParameter)
        #undef TypeValueSetsParameter
        os,
        true
    );
}


// ************************************************************************* //
