/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018 OpenFOAM Foundation
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
    Foam::PatchCollisionDensity

Description
    Function object which generates fields of the number and rate of collisions
    per unit area on all patches. Can optionally take a minimum speed below
    which a collision is not counted.

    Example usage:
    \verbatim
    patchCollisionDensity1
    {
        type        patchCollisionDensity;
        minSpeed    1e-3;
    }
    \endverbatim

SourceFiles
    PatchCollisionDensity.C

\*---------------------------------------------------------------------------*/

#ifndef PatchCollisionDensity_H
#define PatchCollisionDensity_H

#include "CloudFunctionObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class PatchCollisionDensity Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class PatchCollisionDensity
:
    public CloudFunctionObject<CloudType>
{
    // Private data

        typedef typename CloudType::particleType parcelType;

        //- The threshold for a collision
        const scalar minSpeed_;

        //- The field of the number of collisions per unit area
        volScalarField::Boundary collisionDensity_;

        //- The field of the number of collisions per unit area at the last
        //  output
        volScalarField::Boundary collisionDensity0_;

        //- The time at the last output
        scalar time0_;


protected:

    // Protected Member Functions

        //- Write post-processing info
        void write();


public:

    //- Runtime type information
    TypeName("patchCollisionDensity");


    // Constructors

        //- Construct from dictionary
        PatchCollisionDensity
        (
            const dictionary& dict,
            CloudType& owner,
            const word& modelName
        );

        //- Construct copy
        PatchCollisionDensity(const PatchCollisionDensity<CloudType>& ppm);

        //- Construct and return a clone
        virtual autoPtr<CloudFunctionObject<CloudType>> clone() const
        {
            return autoPtr<CloudFunctionObject<CloudType>>
            (
                new PatchCollisionDensity<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~PatchCollisionDensity();


    // Member Functions

        // Evaluation

            //- Post-patch hook
            virtual void postPatch
            (
                const parcelType& p,
                const polyPatch& pp,
                bool& keepParticle
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "PatchCollisionDensity.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
