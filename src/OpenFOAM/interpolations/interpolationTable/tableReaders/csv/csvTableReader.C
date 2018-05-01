/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
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

#include "csvTableReader.H"
#include "fileOperation.H"
#include "DynamicList.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::csvTableReader<Type>::csvTableReader(const dictionary& dict)
:
    tableReader<Type>(dict),
    headerLine_(readBool(dict.lookup("hasHeaderLine"))),
    timeColumn_(readLabel(dict.lookup("timeColumn"))),
    componentColumns_(dict.lookup("valueColumns")),
    separator_(dict.lookupOrDefault<string>("separator", string(","))[0])
{
    if (componentColumns_.size() != pTraits<Type>::nComponents)
    {
        FatalErrorInFunction
            << componentColumns_ << " does not have the expected length "
            << pTraits<Type>::nComponents << endl
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::csvTableReader<Type>::~csvTableReader()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

namespace Foam
{
    // doesn't recognize specialization otherwise
    template<>
    scalar csvTableReader<scalar>::readValue(const List<string>& split)
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
    Type csvTableReader<Type>::readValue(const List<string>& split)
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
}


template<class Type>
void Foam::csvTableReader<Type>::operator()
(
    const fileName& fName,
    List<Tuple2<scalar, Type>>& data
)
{
    // IFstream in(fName);
    autoPtr<ISstream> inPtr(fileHandler().NewIFstream(fName));
    ISstream& in = inPtr();

    DynamicList<Tuple2<scalar, Type>> values;

    // Skip header
    if (headerLine_)
    {
        string line;
        in.getLine(line);
    }

    while (in.good())
    {
        string line;
        in.getLine(line);

        DynamicList<string> split;

        std::size_t pos = 0;
        while (pos != std::string::npos)
        {
            std::size_t nPos = line.find(separator_, pos);

            if (nPos == std::string::npos)
            {
                split.append(line.substr(pos));
                pos=nPos;
            }
            else
            {
                split.append(line.substr(pos, nPos-pos));
                pos=nPos+1;
            }
        }

        if (split.size() <= 1)
        {
            break;
        }

        scalar time = readScalar(IStringStream(split[timeColumn_])());
        Type value = readValue(split);

        values.append(Tuple2<scalar,Type>(time, value));
    }

    data.transfer(values);
}


template<class Type>
void Foam::csvTableReader<Type>::operator()
(
    const fileName& fName,
    List<Tuple2<scalar, List<Tuple2<scalar, Type>>>>& data
)
{
    NotImplemented;
}


template<class Type>
void Foam::csvTableReader<Type>::write(Ostream& os) const
{
    tableReader<Type>::write(os);

    os.writeKeyword("hasHeaderLine")
        << headerLine_ << token::END_STATEMENT << nl;
    os.writeKeyword("timeColumn")
        << timeColumn_ << token::END_STATEMENT << nl;

    // Force writing labelList in ascii
    os.writeKeyword("valueColumns");
    if (os.format() == IOstream::BINARY)
    {
        os.format(IOstream::ASCII);
        os  << componentColumns_;
        os.format(IOstream::BINARY);
    }
    os  << token::END_STATEMENT << nl;

    os.writeKeyword("separator")
        << string(separator_) << token::END_STATEMENT << nl;
}


// ************************************************************************* //
