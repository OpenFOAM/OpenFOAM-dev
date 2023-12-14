/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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
    Foam::DynamicID

Description
    A class that holds the data needed to identify patches and zones.

    The object is identified by name and the index is updated if the mesh has
    changed.

\*---------------------------------------------------------------------------*/

#ifndef DynamicID_H
#define DynamicID_H

#include "word.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators
template<class> class DynamicID;

template<class ObjectType>
Ostream& operator<<(Ostream&, const DynamicID<ObjectType>&);

/*---------------------------------------------------------------------------*\
                          Class DynamicID Declaration
\*---------------------------------------------------------------------------*/

template<class ObjectType>
class DynamicID
{
    // Private Data

        //- Object name
        word name_;

        //- Object index
        label index_;


public:

    // Constructors

        //- Construct from name
        DynamicID(const word& name, const ObjectType& obj)
        :
            name_(name),
            index_(obj.findIndex(name_))
        {}

        //- Construct from Istream
        DynamicID(Istream& is, const ObjectType& obj)
        :
            name_(is),
            index_(obj.findIndex(name_))
        {}


    // Member Functions

        // Access

            //- Return name
            const word& name() const
            {
                return name_;
            }

            //- Return index of first matching zone
            label index() const
            {
                return index_;
            }

            //- Has the zone been found
            bool active() const
            {
                return index_ != -1;
            }


        // Edit

            //- Update
            void update(const ObjectType& obj)
            {
                index_ = obj.findIndex(name_);
            }


    // IOstream Operators

        friend Ostream& operator<< <ObjectType>
        (
            Ostream&,
            const DynamicID<ObjectType>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class ObjectType>
Ostream& operator<<(Ostream& os, const DynamicID<ObjectType>& dynId)
{
    os  << token::BEGIN_LIST
        << dynId.name() << token::SPACE << dynId.index()
        << token::END_LIST;

    // Check state of Ostream
    os.check("Ostream& operator<<(Ostream&, const DynamicID<ObjectType>&)");

    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
