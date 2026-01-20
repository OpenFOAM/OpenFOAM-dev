/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2026 OpenFOAM Foundation
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

template<class Type>
Type CsvReadValue
(
    const typename CsvLabelType<Type>::type& columns,
    const List<string>& split
)
{
    Type result;

    for (label i = 0; i < pTraits<Type>::nComponents; i++)
    {
        if (component(columns, i) >= split.size())
        {
            FatalErrorInFunction
                << "No column " << component(columns, i) << " in "
                << split << endl
                << exit(FatalError);
        }

        result[i] = readScalar
        (
            IStringStream(split[component(columns, i)])()
        );
    }

    return result;
}


template<>
inline label CsvReadValue<label>
(
    const typename CsvLabelType<scalar>::type& columns,
    const List<string>& split
)
{
    if (component(columns, 0) >= split.size())
    {
        FatalErrorInFunction
            << "No column " << component(columns, 0) << " in "
            << split << endl
            << exit(FatalError);
    }

    return readLabel(IStringStream(split[component(columns, 0)])());
}

template<>
inline scalar CsvReadValue<scalar>
(
    const typename CsvLabelType<scalar>::type& columns,
    const List<string>& split
)
{
    if (component(columns, 0) >= split.size())
    {
        FatalErrorInFunction
            << "No column " << component(columns, 0) << " in "
            << split << endl
            << exit(FatalError);
    }

    return readScalar(IStringStream(split[component(columns, 0)])());
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace TableReaders
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Coordinate, class Value>
void Foam::TableReaders::Csv<Coordinate, Value>::read
(
    ISstream& is,
    List<Tuple2<Coordinate, Value>>& data
) const
{
    DynamicList<Tuple2<Coordinate, Value>> values;

    // Skip header
    for (label i = 0; i < nHeaderLine_; i++)
    {
        string line;
        is.getLine(line);
    }

    const label nEntries =
        max(cmptMax(columns_.first()), cmptMax(columns_.second()));

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

        values.append
        (
            Tuple2<Coordinate , Value>
            (
                CsvReadValue<Coordinate>(columns_.first(), split),
                CsvReadValue<Value>(columns_.second(), split)
            )
        );
    }

    data.transfer(values);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Coordinate, class Value>
Foam::TableReaders::Csv<Coordinate, Value>::Csv
(
    const word& name,
    const Function1s::unitConversions& units,
    const dictionary& dict
)
:
    TableFileReader<Coordinate, Value>(units, dict),
    nHeaderLine_(dict.lookup<label>("nHeaderLine")),
    columns_
    (
        dict.found("column")
      ? columnIndices(dict.lookup("columns"))
      : dict.found("refColumn") || dict.found("componentColumns")
      ? columnIndices
        (
            CsvLabelType<Coordinate>()
            (
                dict.lookup<typename CsvLabelType<Coordinate>::oldType>
                (
                    "refColumn"
                )
            ),
            CsvLabelType<Value>()
            (
                dict.lookup<typename CsvLabelType<Value>::oldType>
                (
                    "componentColumns"
                )
            )
        )
      : columnIndices(dict.lookup("columns"))
    ),
    separator_(dict.lookupOrDefault<string>("separator", string(","))[0]),
    mergeSeparators_(readBool(dict.lookup("mergeSeparators")))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Coordinate, class Value>
Foam::TableReaders::Csv<Coordinate, Value>::~Csv()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Coordinate, class Value>
void Foam::TableReaders::Csv<Coordinate, Value>::write
(
    Ostream& os,
    const Function1s::unitConversions& units,
    const List<Tuple2<Coordinate, Value>>& table,
    const word&
) const
{
    TableFileReader<Coordinate, Value>::write(os, units, table);

    writeEntry(os, "nHeaderLine", nHeaderLine_);
    writeEntry(os, "columns", columns_);
    writeEntry(os, "separator", string(separator_));
    writeEntry(os, "mergeSeparators", mergeSeparators_);
}


// ************************************************************************* //
