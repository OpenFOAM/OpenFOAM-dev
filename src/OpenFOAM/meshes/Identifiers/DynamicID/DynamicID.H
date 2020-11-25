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
    Foam::DynamicID

Description
    A class that holds the data needed to identify things (zones, patches)
    in a dynamic mesh.

    The thing is identified by name.
    Its indices are updated if the mesh has changed.

\*---------------------------------------------------------------------------*/

#ifndef DynamicID_H
#define DynamicID_H

#include "wordRe.H"
#include "labelList.H"

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

        //- Zone name
        wordRe key_;

        //- Zone indices
        labelList indices_;


public:

    // Constructors

        //- Construct from name
        DynamicID(const wordRe& key, const ObjectType& obj)
        :
            key_(key),
            indices_(obj.findIndices(key_))
        {}

        //- Construct from Istream
        DynamicID(Istream& is, const ObjectType& obj)
        :
            key_(is),
            indices_(obj.findIndices(key_))
        {}


    // Destructor - default


    // Member Functions

        // Access

            //- Return name
            const wordRe& name() const
            {
                return key_;
            }

            //- Return indices of matching zones
            const labelList& indices() const
            {
                return indices_;
            }

            //- Return index of first matching zone
            label index() const
            {
                return indices_.empty() ? -1 : indices_[0];
            }

            //- Has the zone been found
            bool active() const
            {
                return !indices_.empty();
            }


        // Edit

            //- Update
            void update(const ObjectType& obj)
            {
                indices_ = obj.findIndices(key_);
            }


    // IOstream Operators

        friend Ostream& operator<< <ObjectType>
        (Ostream&, const DynamicID<ObjectType>&);
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
