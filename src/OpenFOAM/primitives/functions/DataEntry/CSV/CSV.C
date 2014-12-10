/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "CSV.H"
#include "DynamicList.H"
#include "IFstream.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

namespace Foam
{
    // doesn't recognize specialization otherwise
    template<>
    scalar CSV<scalar>::readValue(const List<string>& splitted)
    {
        if (componentColumns_[0] >= splitted.size())
        {
            FatalErrorIn("CSV<scalar>::readValue(const List<string>&)")
                << "No column " << componentColumns_[0] << " in "
                << splitted << endl
                << exit(FatalError);
        }

        return readScalar(IStringStream(splitted[componentColumns_[0]])());
    }


    template<class Type>
    Type CSV<Type>::readValue(const List<string>& splitted)
    {
        Type result;

        for (label i = 0; i < pTraits<Type>::nComponents; i++)
        {
            if (componentColumns_[i] >= splitted.size())
            {
                FatalErrorIn("CSV<Type>::readValue(const List<string>&)")
                    << "No column " << componentColumns_[i] << " in "
                    << splitted << endl
                    << exit(FatalError);
            }

            result[i] =
                readScalar(IStringStream(splitted[componentColumns_[i]])());
        }

        return result;
    }
}


template<class Type>
void Foam::CSV<Type>::read()
{
    fileName expandedFile(fName_);
    IFstream is(expandedFile.expand());

    if (!is.good())
    {
        FatalIOErrorIn("CSV<Type>::read()", is)
            << "Cannot open CSV file for reading."
            << exit(FatalIOError);
    }

    DynamicList<Tuple2<scalar, Type> > values;

    // skip header
    for (label i = 0; i < nHeaderLine_; i++)
    {
        string line;
        is.getLine(line);
    }

    label nEntries = max(componentColumns_);

    // read data
    while (is.good())
    {
        string line;
        is.getLine(line);


        label n = 0;
        std::size_t pos = 0;
        DynamicList<string> splitted;

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
                    splitted.append(line.substr(pos));
                    pos = nPos;
                    n++;
                }
                else
                {
                    splitted.append(line.substr(pos, nPos - pos));
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
                    splitted.append(line.substr(pos));
                    pos = nPos;
                    n++;
                }
                else
                {
                    splitted.append(line.substr(pos, nPos - pos));
                    pos = nPos + 1;
                    n++;
                }
            }
        }


        if (splitted.size() <= 1)
        {
            break;
        }

        scalar x = readScalar(IStringStream(splitted[refColumn_])());
        Type value = readValue(splitted);

        values.append(Tuple2<scalar,Type>(x, value));
    }

    this->table_.transfer(values);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::CSV<Type>::CSV
(
    const word& entryName,
    const dictionary& dict,
    const word& ext
)
:
    DataEntry<Type>(entryName),
    TableBase<Type>(entryName, dict.subDict(type() + ext)),
    coeffs_(dict.subDict(type() + ext)),
    nHeaderLine_(readLabel(coeffs_.lookup("nHeaderLine"))),
    refColumn_(readLabel(coeffs_.lookup("refColumn"))),
    componentColumns_(coeffs_.lookup("componentColumns")),
    separator_(coeffs_.lookupOrDefault<string>("separator", string(","))[0]),
    mergeSeparators_(readBool(coeffs_.lookup("mergeSeparators"))),
    fName_(coeffs_.lookup("fileName"))
{
    if (componentColumns_.size() != pTraits<Type>::nComponents)
    {
        FatalErrorIn("Foam::CSV<Type>::CSV(const word&, Istream&)")
            << componentColumns_ << " does not have the expected length of "
            << pTraits<Type>::nComponents << endl
            << exit(FatalError);
    }

    read();

    TableBase<Type>::check();
}


template<class Type>
Foam::CSV<Type>::CSV(const CSV<Type>& tbl)
:
    DataEntry<Type>(tbl),
    TableBase<Type>(tbl),
    nHeaderLine_(tbl.nHeaderLine_),
    refColumn_(tbl.refColumn_),
    componentColumns_(tbl.componentColumns_),
    separator_(tbl.separator_),
    mergeSeparators_(tbl.mergeSeparators_),
    fName_(tbl.fName_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::CSV<Type>::~CSV()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
const Foam::fileName& Foam::CSV<Type>::fName() const
{
    return fName_;
}


// * * * * * * * * * * * * * *  IOStream operators * * * * * * * * * * * * * //

#include "CSVIO.C"


// ************************************************************************* //
