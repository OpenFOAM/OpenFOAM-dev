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
    Foam::surfZoneIdentifier

Description
    An identifier for a surface zone on a meshed surface.

    Similar in concept to a faceZone on the surface, but can also have a
    "geometricType" as well.

SourceFiles
    surfZoneIdentifier.C

\*---------------------------------------------------------------------------*/

#ifndef surfZoneIdentifier_H
#define surfZoneIdentifier_H

#include "word.H"
#include "label.H"
#include "typeInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class dictionary;

// Forward declaration of friend functions and operators

class surfZoneIdentifier;
Ostream& operator<<(Ostream&, const surfZoneIdentifier&);

/*---------------------------------------------------------------------------*\
                     Class surfZoneIdentifier Declaration
\*---------------------------------------------------------------------------*/

class surfZoneIdentifier
{
    // Private Data

        //- Name of zone
        word name_;

        //- Index of zone in surface mesh
        label index_;

        //- Type name of zone
        mutable word geometricType_;


public:

    // Constructors

        //- Construct null
        surfZoneIdentifier();

        //- Construct from components
        surfZoneIdentifier
        (
            const word& name,
            const label index,
            const word& geometricType = word::null

        );

        //- Construct from dictionary
        surfZoneIdentifier
        (
            const word& name,
            const dictionary&,
            const label index
        );

        //- Construct from another zone identifier, resetting the index
        surfZoneIdentifier
        (
            const surfZoneIdentifier&,
            const label index
        );


    //- Destructor
    virtual ~surfZoneIdentifier();


    // Member Functions

        //- Return name
        const word& name() const
        {
            return name_;
        }

        //- Return name for modification
        word& name()
        {
            return name_;
        }

        //- Return the geometric type of the zone
        const word& geometricType() const
        {
            return geometricType_;
        }

        //- Return the geometric type of the zone for modification
        word& geometricType()
        {
            return geometricType_;
        }

        //- Return the index of this zone in the surface mesh
        label index() const
        {
            return index_;
        }

        //- Write surfZoneIdentifier as a dictionary
        void write(Ostream&) const;

        //- Write surfZoneIdentifier as a dictionary
//        void writeDict(Ostream&) const;


    // Member Operators

//        bool operator!=(const surfZoneIdentifier&) const;
//
//        //- compare.
//        bool operator==(const surfZoneIdentifier&) const;

    // Ostream Operator

        friend Ostream& operator<<(Ostream&, const surfZoneIdentifier&);
//        friend Istream& operator>>(Istream&, surfZoneIdentifier&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
