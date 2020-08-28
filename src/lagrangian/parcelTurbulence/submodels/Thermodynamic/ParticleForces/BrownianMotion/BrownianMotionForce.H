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
    Foam::BrownianMotionForce

Description
    Calculates particle Brownian motion force.

    Reference:
    \verbatim
        Li, A., & Ahmadi, G. (1992).
        Dispersion and deposition of spherical particles from point sources
        in a turbulent channel flow.
        Aerosol science and technology,
        16(4), 209-226.
    \endverbatim

SourceFiles
    BrownianMotionForceI.H
    BrownianMotionForce.C

\*---------------------------------------------------------------------------*/

#ifndef BrownianMotionForce_H
#define BrownianMotionForce_H

#include "ParticleForce.H"
#include "Random.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class BrownianMotionForce Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class BrownianMotionForce
:
    public ParticleForce<CloudType>
{
    // Private Data

        //- Reference to the cloud random number generator
        Random& rndGen_;

        //- Molecular free path length [m]
        const scalar lambda_;

        //- Turbulence flag
        bool turbulence_;

        //- Pointer to the turbulence kinetic energy field
        const volScalarField* kPtr_;

        //- Flag that indicates ownership of turbulence k field
        bool ownK_;


    // Private Member Functions

        //- Inverse erf for Gaussian distribution
        scalar erfInv(const scalar y) const;

        //- Return the k field from the turbulence model
        tmp<volScalarField> kModel() const;


public:

    //- Runtime type information
    TypeName("BrownianMotion");


    // Constructors

        //- Construct from mesh
        BrownianMotionForce
        (
            CloudType& owner,
            const fvMesh& mesh,
            const dictionary& dict
        );

        //- Construct copy
        BrownianMotionForce(const BrownianMotionForce& bmf);

        //- Construct and return a clone
        virtual autoPtr<ParticleForce<CloudType>> clone() const
        {
            return autoPtr<ParticleForce<CloudType>>
            (
                new BrownianMotionForce<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~BrownianMotionForce();


    // Member Functions

        // Access

            //- Return const access to the molecular free path length [m]
            inline scalar lambda() const;

            //- Return const access to the turbulence flag
            inline bool turbulence() const;


        // Evaluation

            //- Cache fields
            virtual void cacheFields(const bool store);

            //- Calculate the coupled force
            virtual forceSuSp calcCoupled
            (
                const typename CloudType::parcelType& p,
                const typename CloudType::parcelType::trackingData& td,
                const scalar dt,
                const scalar mass,
                const scalar Re,
                const scalar muc
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "BrownianMotionForceI.H"

#ifdef NoRepository
    #include "BrownianMotionForce.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
