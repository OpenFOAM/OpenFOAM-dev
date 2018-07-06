/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    Foam::ParticleForceList

Description
    List of particle forces

SourceFiles
    ParticleForceListI.H
    ParticleForceList.C

\*---------------------------------------------------------------------------*/

#ifndef ParticleForceList_H
#define ParticleForceList_H

#include "ParticleForce.H"
#include "forceSuSp.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class ParticleForceList Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class ParticleForceList
:
    public PtrList<ParticleForce<CloudType>>
{
    // Private data

        //- Reference to the owner cloud
        CloudType& owner_;

        //- Reference to the mesh database
        const fvMesh& mesh_;

        //- Forces dictionary
        const dictionary dict_;

        //- Calculate coupled forces flag
        bool calcCoupled_;

        //- Calculate non-coupled forces flag
        bool calcNonCoupled_;


public:

    // Constructors

        //- Null constructor
        ParticleForceList(CloudType& owner, const fvMesh& mesh);

        //- Construct from mesh
        ParticleForceList
        (
            CloudType& owner,
            const fvMesh& mesh,
            const dictionary& dict,
            const bool readFields
        );

        //- Construct copy
        ParticleForceList(const ParticleForceList& pfl);


    //- Destructor
    virtual ~ParticleForceList();


    // Member Functions

        // Access

            //- Return const access to the cloud owner
            inline const CloudType& owner() const;

            //- Return references to the cloud owner
            inline CloudType& owner();

            //- Return the mesh database
            inline const fvMesh& mesh() const;

            //- Return the forces dictionary
            inline const dictionary& dict() const;

            //- Set the calcCoupled flag
            inline void setCalcCoupled(bool flag);

            //- Set the calcNonCoupled flag
            inline void setCalcNonCoupled(bool flag);


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

            //- Calculate the non-coupled force
            virtual forceSuSp calcNonCoupled
            (
                const typename CloudType::parcelType& p,
                const typename CloudType::parcelType::trackingData& td,
                const scalar dt,
                const scalar mass,
                const scalar Re,
                const scalar muc
            ) const;

            //- Return the effective mass
            virtual scalar massEff
            (
                const typename CloudType::parcelType& p,
                const typename CloudType::parcelType::trackingData& td,
                const scalar mass
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ParticleForceListI.H"

#ifdef NoRepository
    #include "ParticleForceList.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
