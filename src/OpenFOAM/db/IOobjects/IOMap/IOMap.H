/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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
    Foam::IOMap

Description
    A Map of objects of type \<Type\> with automated input and output.
    Is a global object; i.e. can be read from undecomposed case.

\*---------------------------------------------------------------------------*/

#ifndef IOMap_H
#define IOMap_H

#include "IOList.H"
#include "Map.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                            Class IOMap Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class IOMap
:
    public IOListBase<Map, IOMap, Type>
{

public:

    // Public Static Data

        //- Does this object support re-reading?
        static const bool rereading = true;


    //- Runtime type information
    TypeName("Map");


    // Constructors

        //- Inherit constructors
        using IOListBase<Map, IOMap, Type>::IOListBase;


    // Member Operators

        //- Inherit assignment operators
        using IOListBase<Map, IOMap, Type>::operator=;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
