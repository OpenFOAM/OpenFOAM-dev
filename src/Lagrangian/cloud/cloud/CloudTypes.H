/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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

Struct
    Foam::CloudTypes

Description
    Struct to aid casting of the cloud to a derived type

\*---------------------------------------------------------------------------*/

#ifndef CloudTypes_H
#define CloudTypes_H

#include "cloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                            Class CloudTypes Declaration
\*---------------------------------------------------------------------------*/

template<class ... Clouds>
struct CloudTypes;


template<class Cloud, class ... Clouds>
struct CloudTypes<Cloud, Clouds ...>
{
    static bool isA(const class cloud& c)
    {
        return Foam::isA<Cloud>(c) || CloudTypes<Clouds ...>::isA(c);
    }

    template<unsigned Size>
    static void types(FixedList<word, Size>& result)
    {
        result[Size - sizeof...(Clouds) - 1] = Cloud::typeName;
        CloudTypes<Clouds ...>::types(result);
    }

    static FixedList<word, sizeof...(Clouds) + 1> types()
    {
        FixedList<word, sizeof...(Clouds) + 1> result;
        types(result);
        return result;
    }

    static string typesString(const word& conjunction)
    {
        string result = '\'' + Cloud::typeName + '\'';
        if (sizeof...(Clouds))
        {
            FixedList<word, sizeof...(Clouds) + 1> types =
                CloudTypes<Cloud, Clouds ...>::types();
            for (label i = 1; i < types.size() - 1; ++ i)
            {
                result += ", '" + types[i] + '\'';
            }
            result += ' ' + conjunction + " '" + types.last() + '\'';
        }
        return result;
    }
};


template<>
struct CloudTypes<>
{
    static bool isA(const class cloud& c)
    {
        return false;
    }

    template<unsigned Size>
    static void types(FixedList<word, Size>& result)
    {}
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
