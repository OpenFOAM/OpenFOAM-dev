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
    Foam::NoStochasticCollision

Description
    Dummy collision model for 'none'

\*---------------------------------------------------------------------------*/

#ifndef NoStochasticCollision_H
#define NoStochasticCollision_H

#include "StochasticCollisionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                    Class NoStochasticCollision Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class NoStochasticCollision
:
    public StochasticCollisionModel<CloudType>
{
protected:

    // Protected Member Functions

        //- Update the model
        virtual void collide
        (
            typename CloudType::parcelType::trackingData& td,
            const scalar dt
        );


public:

    //- Runtime type information
    TypeName("none");


    // Constructors

        //- Construct from dictionary
        NoStochasticCollision(const dictionary& dict, CloudType& owner);

        //- Construct copy
        NoStochasticCollision(const NoStochasticCollision<CloudType>& cm);

        //- Construct and return a clone
        virtual autoPtr<StochasticCollisionModel<CloudType>> clone() const
        {
            return autoPtr<StochasticCollisionModel<CloudType>>
            (
                new NoStochasticCollision<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~NoStochasticCollision();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "NoStochasticCollision.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
