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
    Foam::PairCollision

Description

SourceFiles
    PairCollision.C

\*---------------------------------------------------------------------------*/

#ifndef PairCollision_H
#define PairCollision_H

#include "CollisionModel.H"
#include "InteractionLists.H"
#include "WallSiteData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

template<class CloudType>
class PairModel;

template<class CloudType>
class WallModel;


/*---------------------------------------------------------------------------*\
                        Class PairCollision Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class PairCollision
:
    public CollisionModel<CloudType>
{
    // Static data

        //- Tolerance to determine flat wall interactions
        static scalar cosPhiMinFlatWall;

        //- Distance to multiply the distance to a wall with within
        //  which a duplicate flat wall interaction is defined.  If
        //  cosPhiFlatWall = 1 - e, this = sqrt(2e - e^2) < sqrt(3e)
        //  for all e < 1, giving a safe tolerance (i.e. no double
        //  interactions next to each other.)
        static scalar flatWallDuplicateExclusion;


    // Private Data

        //- PairModel to calculate the interaction between two parcels
        autoPtr<PairModel<CloudType>> pairModel_;

        //- WallModel to calculate the interaction between the parcel and walls
        autoPtr<WallModel<CloudType>> wallModel_;

        //- Interactions lists determining which cells are in
        //  interaction range of each other
        InteractionLists<typename CloudType::parcelType> il_;


    // Private Member Functions

        //- Pre collision tasks
        void preInteraction();

        //- Interactions between parcels
        void parcelInteraction();

        //- Interactions between real (on-processor) particles
        void realRealInteraction();

        //- Interactions between real and referred (off processor) particles
        void realReferredInteraction();

        //- Interactions with walls
        void wallInteraction();

        bool duplicatePointInList
        (
            const DynamicList<point>& existingPoints,
            const point& pointToTest,
            scalar duplicateRangeSqr
        ) const;

        bool duplicatePointInList
        (
            const DynamicList<point>& existingPoints,
            const point& pointToTest,
            const scalarList& duplicateRangeSqr
        ) const;

        //- Post collision tasks
        void postInteraction();

        //- Calculate the pair force between parcels
        void evaluatePair
        (
            typename CloudType::parcelType& pA,
            typename CloudType::parcelType& pB
        ) const;

        //- Calculate the wall forces on a parcel
        void evaluateWall
        (
            typename CloudType::parcelType& p,
            const List<point>& flatSitePoints,
            const List<WallSiteData<vector>>& flatSiteData,
            const List<point>& sharpSitePoints,
            const List<WallSiteData<vector>>& sharpSiteData
        ) const;


public:

    //- Runtime type information
    TypeName("pairCollision");


    // Constructors

        //- Construct from components
        PairCollision(const dictionary& dict, CloudType& owner);

        //- Construct copy
        PairCollision(const PairCollision<CloudType>& cm);

        //- Construct and return a clone
        virtual autoPtr<CollisionModel<CloudType>> clone() const
        {
            return autoPtr<CollisionModel<CloudType>>
            (
                new PairCollision<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~PairCollision();


    // Member Functions

        //- Return the number of times to subcycle the current
        //  timestep to meet the criteria of the collision model.
        virtual label nSubCycles() const;

        // Collision function
        virtual void collide();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "PairCollision.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
