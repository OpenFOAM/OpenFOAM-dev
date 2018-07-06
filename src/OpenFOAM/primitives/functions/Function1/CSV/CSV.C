/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
//#include "IFstream.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<>
Foam::label Foam::Function1Types::CSV<Foam::label>::readValue
(
    const List<string>& split
)
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
Foam::scalar Foam::Function1Types::CSV<Foam::scalar>::readValue
(
    const List<string>& split
)
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
Type Foam::Function1Types::CSV<Type>::readValue(const List<string>& split)
{
    Type result;

    for (label i = 0; i < pTraits<Type>::nComponents; i++)
    {
        if (componentColumns_[i] >= split.size())
        {
            FatalErrorInFunction
            << "No column " << componentColumns_[i] << " in "
                << split << endl
                << exit(FatalError);
        }

        result[i] =
        readScalar(IStringStream(split[componentColumns_[i]])());
    }

    return result;
}


template<class Type>
void Foam::Function1Types::CSV<Type>::read()
{
    fileName expandedFile(fName_);
    // IFstream is(expandedFile.expand());
    autoPtr<ISstream> isPtr(fileHandler().NewIFstream(expandedFile.expand()));
    ISstream& is = isPtr();

    if (!is.good())
    {
        FatalIOErrorInFunction(is)
            << "Cannot open CSV file for reading."
            << exit(FatalIOError);
    }

    DynamicList<Tuple2<scalar, Type>> values;

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

    this->table_.transfer(values);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::Function1Types::CSV<Type>::CSV
(
    const word& entryName,
    const dictionary& dict
)
:
    TableBase<Type>(entryName, dict),
    nHeaderLine_(readLabel(dict.lookup("nHeaderLine"))),
    refColumn_(readLabel(dict.lookup("refColumn"))),
    componentColumns_(dict.lookup("componentColumns")),
    separator_(dict.lookupOrDefault<string>("separator", string(","))[0]),
    mergeSeparators_(readBool(dict.lookup("mergeSeparators"))),
    fName_(dict.lookup("file"))
{
    if (componentColumns_.size() != pTraits<Type>::nComponents)
    {
        FatalErrorInFunction
            << componentColumns_ << " does not have the expected length of "
            << pTraits<Type>::nComponents << endl
            << exit(FatalError);
    }

    read();

    TableBase<Type>::check();
}


template<class Type>
Foam::Function1Types::CSV<Type>::CSV(const CSV<Type>& tbl)
:
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
Foam::Function1Types::CSV<Type>::~CSV()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
const Foam::fileName& Foam::Function1Types::CSV<Type>::fName() const
{
    return fName_;
}


template<class Type>
void Foam::Function1Types::CSV<Type>::writeData(Ostream& os) const
{
    Function1<Type>::writeData(os);
    os  << token::END_STATEMENT << nl;
    os  << indent << word(this->name() + "Coeffs") << nl;
    os  << indent << token::BEGIN_BLOCK << incrIndent << nl;

    // Note: for TableBase write the dictionary entries it needs but not
    // the values themselves
    TableBase<Type>::writeEntries(os);

    os.writeKeyword("nHeaderLine") << nHeaderLine_ << token::END_STATEMENT
        << nl;

    os.writeKeyword("refColumn") << refColumn_ << token::END_STATEMENT << nl;

    componentColumns_.writeEntry("componentColumns", os);

    os.writeKeyword("separator") << string(separator_)
        << token::END_STATEMENT << nl;
    os.writeKeyword("mergeSeparators") << mergeSeparators_
        << token::END_STATEMENT << nl;
    os.writeKeyword("file") << fName_ << token::END_STATEMENT << nl;
    os  << decrIndent << indent << token::END_BLOCK << endl;
}


// ************************************************************************* //
