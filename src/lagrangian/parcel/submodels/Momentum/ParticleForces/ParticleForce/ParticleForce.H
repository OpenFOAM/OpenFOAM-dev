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
    Foam::ParticleForce

Description
    Abstract base class for particle forces

SourceFiles
    ParticleForceI.H
    ParticleForce.C
    ParticleForceNew.C

\*---------------------------------------------------------------------------*/

#ifndef ParticleForce_H
#define ParticleForce_H

#include "dictionary.H"
#include "forceSuSp.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class ParticleForce Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class ParticleForce
{
    // Private Data

        //- Reference to the owner cloud
        CloudType& owner_;

        //- Reference to the mesh database
        const fvMesh& mesh_;

        //- Force coefficients dictionary
        const dictionary coeffs_;


public:

    //- Runtime type information
    TypeName("particleForce");

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        ParticleForce,
        dictionary,
        (
            CloudType& owner,
            const fvMesh& mesh,
            const dictionary& dict
        ),
        (owner, mesh, dict)
    );


    //- Convenience typedef for return type
    typedef VectorSpace<Vector<vector>, vector, 2> returnType;


    // Constructors

        //- Construct from mesh
        ParticleForce
        (
            CloudType& owner,
            const fvMesh& mesh,
            const dictionary& dict,
            const word& forceType,
            const bool readCoeffs
        );

        //- Construct copy
        ParticleForce(const ParticleForce& pf);

        //- Construct and return a clone
        virtual autoPtr<ParticleForce<CloudType>> clone() const
        {
            return autoPtr<ParticleForce<CloudType>>
            (
                new ParticleForce<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~ParticleForce();


    //- Selector
    static autoPtr<ParticleForce<CloudType>> New
    (
        CloudType& owner,
        const fvMesh& mesh,
        const dictionary& dict,
        const word& forceType
    );


    // Member Functions

        // Access

            //- Return const access to the cloud owner
            inline const CloudType& owner() const;

            //- Return references to the cloud owner
            inline CloudType& owner();

            //- Return the mesh database
            inline const fvMesh& mesh() const;

            //- Return the force coefficients dictionary
            inline const dictionary& coeffs() const;


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

            //- Return the added mass
            virtual scalar massAdd
            (
                const typename CloudType::parcelType& p,
                const typename CloudType::parcelType::trackingData& td,
                const scalar mass
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ParticleForceI.H"

#ifdef NoRepository
    #include "ParticleForce.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeParticleForceModel(CloudType)                                      \
                                                                               \
    typedef Foam::CloudType::momentumCloudType momentumCloudType;            \
    defineNamedTemplateTypeNameAndDebug                                        \
        (Foam::ParticleForce<momentumCloudType>, 0);                          \
                                                                               \
    namespace Foam                                                             \
    {                                                                          \
        defineTemplateRunTimeSelectionTable                                    \
        (                                                                      \
            ParticleForce<momentumCloudType>,                                 \
            dictionary                                                         \
        );                                                                     \
    }


#define makeParticleForceModelType(SS, CloudType)                              \
                                                                               \
    typedef Foam::CloudType::momentumCloudType momentumCloudType;            \
    defineNamedTemplateTypeNameAndDebug(Foam::SS<momentumCloudType>, 0);      \
                                                                               \
    Foam::ParticleForce<momentumCloudType>::                                  \
        adddictionaryConstructorToTable<Foam::SS<momentumCloudType>>          \
            add##SS##CloudType##momentumCloudType##ConstructorToTable_;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //
