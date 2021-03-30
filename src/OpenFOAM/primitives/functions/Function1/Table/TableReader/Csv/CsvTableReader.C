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
    if (componentColumns_[0] >= split.size())
    {
        FatalErrorInFunction
            << "No column " << componentColumns_[0] << " in "
            << split << endl
            << exit(FatalError);
    }

    return readLabel(IStringStream(split[componentColumns_[0]])());
}


template<>
scalar Csv<scalar>::readValue(const List<string>& split) const
{
    if (componentColumns_[0] >= split.size())
    {
        FatalErrorInFunction
            << "No column " << componentColumns_[0] << " in "
            << split << endl
            << exit(FatalError);
    }

    return readScalar(IStringStream(split[componentColumns_[0]])());
}


template<class Type>
Type Csv<Type>::readValue(const List<string>& split) const
{
    Type result;

    for(label i = 0;i < pTraits<Type>::nComponents; i++)
    {
        if (componentColumns_[i] >= split.size())
        {
            FatalErrorInFunction
                << "No column " << componentColumns_[i] << " in "
                << split << endl
                << exit(FatalError);
        }

        result[i] = readScalar
        (
            IStringStream(split[componentColumns_[i]])()
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

    const label nEntries = max(refColumn_, max(componentColumns_));

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

        scalar x = readScalar(IStringStream(split[refColumn_])());
        Type value = readValue(split);

        values.append(Tuple2<scalar,Type>(x, value));
    }

    data.transfer(values);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::TableReaders::Csv<Type>::Csv
(
    const word& name,
    const dictionary& dict,
    List<Tuple2<scalar, Type>>& table
)
:
    TableFileReader<Type>(dict),
    nHeaderLine_(dict.lookup<label>("nHeaderLine")),
    refColumn_(dict.lookup<label>("refColumn")),
    componentColumns_(dict.lookup("componentColumns")),
    separator_(dict.lookupOrDefault<string>("separator", string(","))[0]),
    mergeSeparators_(readBool(dict.lookup("mergeSeparators")))
{
    if (componentColumns_.size() != pTraits<Type>::nComponents)
    {
        FatalErrorInFunction
            << componentColumns_ << " does not have the expected length "
            << pTraits<Type>::nComponents << endl
            << exit(FatalError);
    }

    TableFileReader<Type>::read(dict, table);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::TableReaders::Csv<Type>::~Csv()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::TableReaders::Csv<Type>::write
(
    Ostream& os,
    const List<Tuple2<scalar, Type>>& table
) const
{
    TableFileReader<Type>::write(os, table);

    writeEntry(os, "nHeaderLine", nHeaderLine_);
    writeEntry(os, "refColumn", refColumn_);
    writeEntry(os, "componentColumns", componentColumns_);
    writeEntry(os, "separator", string(separator_));
    writeEntry(os, "mergeSeparators", mergeSeparators_);
}


// ************************************************************************* //
