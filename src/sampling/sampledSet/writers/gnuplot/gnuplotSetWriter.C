/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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

#include "gnuplotSetWriter.H"
#include "coordSet.H"
#include "fileName.H"
#include "OFstream.H"
#include "OSspecific.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(gnuplotSetWriter, 0);
    addToRunTimeSelectionTable(setWriter, gnuplotSetWriter, word);
    addToRunTimeSelectionTable(setWriter, gnuplotSetWriter, dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::gnuplotSetWriter::~gnuplotSetWriter()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::gnuplotSetWriter::write
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
    if (!set.hasScalarAxis())
    {
        FatalErrorInFunction
            << "Cannot write " << setName << " in " << typeName
            << " format as it does not have a scalar axis"
            << exit(FatalError);
    }

    if (!isDir(outputDir))
    {
        mkDir(outputDir);
    }

    OFstream os
    (
        outputDir/setName + ".gplt",
        IOstream::ASCII,
        IOstream::currentVersion,
        writeCompression_
    );

    os << "$data << EOD" << nl;

    writeTable
    (
        set.axis() == coordSet::axisTypeNames_[coordSet::axisType::DEFAULT]
      ? coordSet
        (
            set.segments(),
            word::null,
            pointField::null(),
            set.scalarName(),
            set.scalarCoords()
        )
      : set,
        #define TypeValueSetsParameter(Type, nullArg) Type##ValueSets,
        FOR_ALL_FIELD_TYPES(TypeValueSetsParameter)
        #undef TypeValueSetsParameter
        os
    );

    os << nl << "EOD" << nl << nl;

    os  << "set term postscript color" << nl
        << "set output \"" << setName.c_str() << ".ps\"" << nl
        << "set xlabel \"" << set.scalarName() << "\"" << nl;

    if (valueSetNames.size() == 2)
    {
        os  << "set ytics nomirror" << nl
            << "set ylabel \"" << valueSetNames[0] << "\"" << nl
            << "set y2tics nomirror" << nl
            << "set y2label \"" << valueSetNames[1] << "\"" << nl;
    }

    os << nl << "plot ";

    label columni = 0;

    forAll(valueSetNames, fieldi)
    {
        #define PlotTypeValueSets(Type, nullArg)                            \
            if (Type##ValueSets.set(fieldi))                                \
            {                                                               \
                for                                                         \
                (                                                           \
                    direction cmpt = 0;                                     \
                    cmpt < pTraits<Type>::nComponents;                      \
                    cmpt++                                                  \
                )                                                           \
                {                                                           \
                    const bool separator =                                  \
                        !valueSetNames[fieldi].empty()                      \
                     && strlen(pTraits<Type>::componentNames[cmpt]) > 0;    \
                                                                            \
                    const word w =                                          \
                        valueSetNames[fieldi]                               \
                      + (separator ? "_" : "")                              \
                      + pTraits<Type>::componentNames[cmpt];                \
                                                                            \
                    if (columni != 0) os << ", ";                           \
                                                                            \
                    os  << "$data us 1:" << 2 + columni;                    \
                                                                            \
                    if (valueSetNames.size() == 2)                          \
                    {                                                       \
                        os  << " axis x1y" << 1 + fieldi;                   \
                    }                                                       \
                                                                            \
                    os  << " title \"" << w << "\" with lines";             \
                                                                            \
                    ++ columni;                                             \
                }                                                           \
            }
        FOR_ALL_FIELD_TYPES(PlotTypeValueSets)
        #undef PlotTypeValueSets
    }
}


// ************************************************************************* //
