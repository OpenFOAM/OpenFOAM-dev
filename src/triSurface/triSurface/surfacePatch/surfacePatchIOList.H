/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
    Foam::surfacePatchIOList

Description
    IOobject for a surfacePatchList

SourceFiles
    surfacePatchIOList.C

\*---------------------------------------------------------------------------*/

#ifndef surfacePatchIOList_H
#define surfacePatchIOList_H

#include "surfacePatchList.H"
#include "regIOobject.H"
#include "faceList.H"
#include "className.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

class surfacePatchIOList;

Ostream& operator<<(Ostream&, const surfacePatchIOList&);


/*---------------------------------------------------------------------------*\
                     Class surfacePatchIOList Declaration
\*---------------------------------------------------------------------------*/

class surfacePatchIOList
:
    public surfacePatchList,
    public regIOobject
{
public:

    //- Runtime type information
    TypeName("surfacePatchIOList");


    // Constructors

        //- Construct from IOobject
        explicit surfacePatchIOList(const IOobject& io);

        //- Construct from IOobject
        surfacePatchIOList(const IOobject& io, const surfacePatchList&);

        //- Disallow default bitwise copy construction
        surfacePatchIOList(const surfacePatchIOList&) = delete;


    //- Destructor
    ~surfacePatchIOList();


    // Member Functions

        //- writeData member function required by regIOobject
        bool writeData(Ostream&) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const surfacePatchIOList&) = delete;


    // IOstream Operators

        friend Ostream& operator<<(Ostream&, const surfacePatchIOList&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
