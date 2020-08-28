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
    Foam::surfZoneIOList

Description
    IOobject for a surfZoneList

SourceFiles
    surfZoneIOList.C

\*---------------------------------------------------------------------------*/

#ifndef surfZoneIOList_H
#define surfZoneIOList_H

#include "surfZoneList.H"
#include "regIOobject.H"
#include "faceList.H"
#include "className.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

class surfZoneIOList;

Ostream& operator<<(Ostream&, const surfZoneIOList&);


/*---------------------------------------------------------------------------*\
                       Class surfZoneIOList Declaration
\*---------------------------------------------------------------------------*/

class surfZoneIOList
:
    public surfZoneList,
    public regIOobject
{
public:

    //- Runtime type information
    TypeName("surfZoneList");


    // Constructors

        //- Construct from IOobject
        explicit surfZoneIOList(const IOobject&);

        //- Construct from IOobject and surfZoneList
        surfZoneIOList(const IOobject&, const surfZoneList&);

        //- Construct from IOobject and surfZoneList
        surfZoneIOList(const IOobject&, surfZoneList&&);

        //- Disallow default bitwise copy construction
        surfZoneIOList(const surfZoneIOList&) = delete;


    //- Destructor
    ~surfZoneIOList();


    // Member Functions

        //- writeData member function required by regIOobject
        bool writeData(Ostream&) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const surfZoneIOList&) = delete;


    // IOstream Operators

        friend Ostream& operator<<(Ostream&, const surfZoneIOList&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
