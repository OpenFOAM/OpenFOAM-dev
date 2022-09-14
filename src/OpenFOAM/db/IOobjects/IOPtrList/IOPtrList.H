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
    Foam::IOPtrList

Description
    A PtrList of objects of type \<Type\> with automated input and output.

SourceFiles
    IOPtrList.C

\*---------------------------------------------------------------------------*/

#ifndef IOPtrList_H
#define IOPtrList_H

#include "IOList.H"
#include "PtrList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class IOPtrList Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class IOPtrList
:
    public IOListBase<PtrList, IOPtrList, Type>
{

public:

    // Public Static Data

        //- Does this object support re-reading?
        static const bool rereading = true;


    //- Runtime type information
    TypeName("PtrList");


    // Constructors

        //- Inherit constructors
        using IOListBase<PtrList, IOPtrList, Type>::IOListBase;

        //- Construct from IOobject using given Istream constructor class
        template<class INew>
        IOPtrList(const IOobject&, const INew&);


    // Member Operators

        //- Inherit assignment operators
        using IOListBase<PtrList, IOPtrList, Type>::operator=;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "IOPtrList.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
