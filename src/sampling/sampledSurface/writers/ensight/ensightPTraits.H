/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2018 OpenFOAM Foundation
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

Class
    Foam::ensightPTraits

Description
    Conversion of OpenFOAM pTraits into the Ensight equivalent

\*---------------------------------------------------------------------------*/

#ifndef ensightPTraits_H
#define ensightPTraits_H

#include "pTraits.H"
#include "fieldTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class ensightPTraits Declaration
\*---------------------------------------------------------------------------*/

template<class PrimitiveType>
class ensightPTraits
{
public:

    // Static data members

        static const char* const typeName;

};


template<>
class ensightPTraits<scalar>
{
public:

    static const char* const typeName;
};

template<>
class ensightPTraits<vector>
{
public:

    static const char* const typeName;
};

template<>
class ensightPTraits<sphericalTensor>
{
public:

    static const char* const typeName;
};

template<>
class ensightPTraits<symmTensor>
{
public:

    static const char* const typeName;
};

template<>
class ensightPTraits<tensor>
{
public:

    static const char* const typeName;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
