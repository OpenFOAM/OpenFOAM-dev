/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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

#include "CsvTableReader.H"
#include "DynamicList.H"
#include "Field.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

namespace Foam
{
namespace TableReaders
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<>
label Csv<label>::readValue(const List<string>& split) const
{
    if (component(columns_.second(), 0) >= split.size())
    {
        FatalErrorInFunction
            << "No column " << component(columns_.second(), 0) << " in "
            << split << endl
            << exit(FatalError);
    }

    return readLabel(IStringStream(split[component(columns_.second(), 0)])());
}


template<>
scalar Csv<scalar>::readValue(const List<string>& split) const
{
    if (component(columns_.second(), 0) >= split.size())
    {
        FatalErrorInFunction
            << "No column " << component(columns_.second(), 0) << " in "
            << split << endl
            << exit(FatalError);
    }

    return readScalar(IStringStream(split[component(columns_.second(), 0)])());
}


template<class Type>
Type Csv<Type>::readValue(const List<string>& split) const
{
    Type result;

    for (label i = 0; i < pTraits<Type>::nComponents; i++)
    {
        if (component(columns_.second(), i) >= split.size())
        {
            FatalErrorInFunction
                << "No column " << component(columns_.second(), i) << " in "
                << split << endl
                << exit(FatalError);
        }

        result[i] = readScalar
        (
            IStringStream(split[component(columns_.second(), i)])()
        );
    }

    return result;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace TableReaders
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void Foam::TableReaders::Csv<Type>::read
(
    ISstream& is,
    List<Tuple2<scalar, Type>>& data
) const
{
    DynamicList<Tuple2<scalar, Type>> values;

    // Skip header
    for (label i = 0; i < nHeaderLine_; i++)
    {
        string line;
        is.getLine(line);
    }

    const label nEntries = max(columns_.first(), cmptMax(columns_.second()));

    // Read data
    while (is.good())
    {
        string line;
        is.getLine(line);

        label n = 0;
        std::size_t pos = 0;
        DynamicList<string> split;

        if (mergeSeparators_)
        {
            std::size_t nPos = 0;

            while ((pos != std::string::npos) && (n <= nEntries))
            {
                bool found = false;
                while (!found)
                {
                    nPos = line.find(separator_, pos);

                    if ((nPos != std::string::npos) && (nPos - pos == 0))
                    {
                        pos = nPos + 1;
                    }
                    else
                    {
                        found = true;
                    }
                }

                nPos = line.find(separator_, pos);

                if (nPos == std::string::npos)
                {
                    split.append(line.substr(pos));
                    pos = nPos;
                    n++;
                }
                else
                {
                    split.append(line.substr(pos, nPos - pos));
                    pos = nPos + 1;
                    n++;
                }
            }
        }
        else
        {
            while ((pos != std::string::npos) && (n <= nEntries))
            {
                std::size_t nPos = line.find(separator_, pos);

                if (nPos == std::string::npos)
                {
                    split.append(line.substr(pos));
                    pos = nPos;
                    n++;
                }
                else
                {
                    split.append(line.substr(pos, nPos - pos));
                    pos = nPos + 1;
                    n++;
                }
            }
        }

        if (split.size() <= 1)
        {
            break;
        }

        scalar x = readScalar(IStringStream(split[columns_.first()])());
        Type value = readValue(split);

        values.append(Tuple2<scalar, Type>(x, value));
    }

    data.transfer(values);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::TableReaders::Csv<Type>::Csv
(
    const word& name,
    const Function1s::unitConversions& units,
    const dictionary& dict
)
:
    TableFileReader<Type>(units, dict),
    nHeaderLine_(dict.lookup<label>("nHeaderLine")),
    columns_
    (
        !dict.found("refColumn") || !dict.found("componentColumns")
      ? columnIndices(dict.lookup("columns"))
      : columnIndices
        (
            dict.lookup<label>("refColumn"),
            CsvLabelType<Type>()
            (
                dict.lookup<typename CsvLabelType<Type>::oldType>
                (
                    "componentColumns"
                )
            )
        )
    ),
    separator_(dict.lookupOrDefault<string>("separator", string(","))[0]),
    mergeSeparators_(readBool(dict.lookup("mergeSeparators")))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::TableReaders::Csv<Type>::~Csv()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::TableReaders::Csv<Type>::write
(
    Ostream& os,
    const Function1s::unitConversions& units,
    const List<Tuple2<scalar, Type>>& table,
    const word&
) const
{
    TableFileReader<Type>::write(os, units, table);

    writeEntry(os, "nHeaderLine", nHeaderLine_);
    writeEntry(os, "columns", columns_);
    writeEntry(os, "separator", string(separator_));
    writeEntry(os, "mergeSeparators", mergeSeparators_);
}


// ************************************************************************* //
