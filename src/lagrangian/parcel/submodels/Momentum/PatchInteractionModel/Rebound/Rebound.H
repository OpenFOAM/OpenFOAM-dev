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
    Foam::Rebound

Description
    Simple rebound patch interaction model

\*---------------------------------------------------------------------------*/

#ifndef Rebound_H
#define Rebound_H

#include "PatchInteractionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                           Class Rebound Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class Rebound
:
    public PatchInteractionModel<CloudType>
{
    // Private Data

       //- Factor applied to velocity on rebound
       //  Normal rebound = 1
       scalar UFactor_;


public:

    //- Runtime type information
    TypeName("rebound");


    // Constructors

        //- Construct from dictionary
        Rebound(const dictionary& dict, CloudType& cloud);

        //- Construct copy
        Rebound(const Rebound<CloudType>& pim);

        //- Construct and return a clone
        virtual autoPtr<PatchInteractionModel<CloudType>> clone() const
        {
            return autoPtr<PatchInteractionModel<CloudType>>
            (
                new Rebound<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~Rebound();


    // Member Functions
        //- Apply velocity correction
        //  Returns true if particle remains in same cell
        virtual bool correct
        (
            typename CloudType::parcelType& p,
            const polyPatch& pp,
            bool& keepParticle
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "Rebound.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
