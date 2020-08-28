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
    Foam::CollisionModel

Description
    Templated collision model class.

SourceFiles
    CollisionModel.C
    CollisionModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef CollisionModel_H
#define CollisionModel_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "CloudSubModelBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class CollisionModel Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class CollisionModel
:
    public CloudSubModelBase<CloudType>
{
    // Protected data

        //- Convenience typedef for parcel type
        typedef typename CloudType::parcelType parcelType;


public:

    //- Runtime type information
    TypeName("collisionModel");

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        CollisionModel,
        dictionary,
        (
            const dictionary& dict,
            CloudType& owner
        ),
        (dict, owner)
    );


    // Constructors

        //- Construct null from owner
        CollisionModel(CloudType& owner);

        //- Construct from components
        CollisionModel
        (
            const dictionary& dict,
            CloudType& owner,
            const word& type
        );

        //- Construct copy
        CollisionModel(const CollisionModel<CloudType>& cm);

        //- Construct and return a clone
        virtual autoPtr<CollisionModel<CloudType>> clone() const = 0;


    //- Destructor
    virtual ~CollisionModel();


    //- Selector
    static autoPtr<CollisionModel<CloudType>> New
    (
        const dictionary& dict,
        CloudType& owner
    );


    // Member Functions

        //- Return the number of times to subcycle the current
        //  timestep to meet the criteria of the collision model
        virtual label nSubCycles() const = 0;

        // Collision function
        virtual void collide() = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeCollisionModel(CloudType)                                          \
                                                                               \
    typedef Foam::CloudType::collidingCloudType collidingCloudType;            \
    defineNamedTemplateTypeNameAndDebug                                        \
    (                                                                          \
        Foam::CollisionModel<collidingCloudType>,                              \
        0                                                                      \
    );                                                                         \
    namespace Foam                                                             \
    {                                                                          \
        defineTemplateRunTimeSelectionTable                                    \
        (                                                                      \
            CollisionModel<collidingCloudType>,                                \
            dictionary                                                         \
        );                                                                     \
    }


#define makeCollisionModelType(SS, CloudType)                                  \
                                                                               \
    typedef Foam::CloudType::collidingCloudType collidingCloudType;            \
    defineNamedTemplateTypeNameAndDebug(Foam::SS<collidingCloudType>, 0);      \
                                                                               \
    Foam::CollisionModel<collidingCloudType>::                                 \
        adddictionaryConstructorToTable<Foam::SS<collidingCloudType>>          \
            add##SS##CloudType##collidingCloudType##ConstructorToTable_;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "CollisionModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
