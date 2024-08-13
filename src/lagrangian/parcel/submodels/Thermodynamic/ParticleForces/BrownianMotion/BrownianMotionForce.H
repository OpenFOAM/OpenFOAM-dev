/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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
#include "randomGenerator.H"

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

        //- Molecular free path length [m]
        const scalar lambda_;


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


        // Evaluation

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
