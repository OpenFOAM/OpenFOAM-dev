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
    Foam::objectMap

Description
    An objectMap is a pair of labels defining the mapping of an object from
    another object, e.g. a cell mapped from a point.

SourceFiles
    objectMapI.H

\*---------------------------------------------------------------------------*/

#ifndef objectMap_H
#define objectMap_H

#include "labelList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

class objectMap;
inline bool operator==(const objectMap& a, const objectMap& b);
inline bool operator!=(const objectMap& a, const objectMap& b);
inline Ostream& operator<<(Ostream&, const objectMap&);
inline Istream& operator>>(Istream&, objectMap&);


/*---------------------------------------------------------------------------*\
                          Class objectMap Declaration
\*---------------------------------------------------------------------------*/

class objectMap
{
    // Private Data

        //- Object index
        label index_;

        //- Master object index
        labelList masterObjects_;


public:

    // Constructors

        //- Null constructor for lists
        inline objectMap();

        //- Construct from components
        inline objectMap(const label index, const labelList& master);

        //- Construct from Istream
        inline objectMap(Istream&);


    // Member Functions

        //- Return object index
        inline label& index();
        inline label index() const;

        //- Return master object index
        inline labelList& masterObjects();
        inline const labelList& masterObjects() const;


    // Friend Operators

        friend bool operator==(const objectMap& a, const objectMap& b);
        friend bool operator!=(const objectMap& a, const objectMap& b);


    // IOstream Operators

        friend Ostream& operator<<(Ostream&, const objectMap&);

        friend Istream& operator>>(Istream&, objectMap&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "objectMapI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
