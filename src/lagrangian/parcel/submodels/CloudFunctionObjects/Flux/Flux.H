/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022 OpenFOAM Foundation
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
    Foam::Flux

Description
    Function objects which generate the number, volume or mass flux of
    particles in a cloud.

    Example usage:
    \verbatim
    massFlux1
    {
        type        numberFlux; // numberFlux, volumeFlux, or massFlux
    }
    \endverbatim

SourceFiles
    Flux.C

\*---------------------------------------------------------------------------*/

#ifndef Flux_H
#define Flux_H

#include "NamedEnum.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "cloud.H"
#include "CloudFunctionObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                            Class Flux Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType, class Derived>
class Flux
:
    public CloudFunctionObject<CloudType>
{
    // Private Typedefs

        //- The parcel type
        typedef typename CloudType::particleType parcelType;


    // Private Data

        //- Whether or not to write the flux field
        const bool write_;

        //- The flux field
        surfaceScalarField phi_;


protected:

    // Protected Member Functions

        //- Write post-processing info
        void write();

        //- Accumulate the particle into the flux
        void accumulate(const parcelType& p, const bool isPre);


public:

    // Constructors

        //- Construct from dictionary
        Flux
        (
            const dictionary& dict,
            CloudType& owner,
            const word& modelName
        );

        //- Construct copy
        Flux(const Flux<CloudType, Derived>& ppm);

        //- Construct and return a clone
        virtual autoPtr<CloudFunctionObject<CloudType>> clone() const
        {
            return autoPtr<CloudFunctionObject<CloudType>>
            (
                new Derived(static_cast<const Derived&>(*this))
            );
        }


    //- Destructor
    virtual ~Flux();


    // Member Functions

        // Evaluation

            //- Pre-evolve hook
            virtual void preEvolve();

            //- Pre-face hook
            virtual void preFace(const parcelType& p);

            //- Post-face hook
            virtual void postFace(const parcelType& p);
};


/*---------------------------------------------------------------------------*\
                            Class NumberFlux Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class NumberFlux
:
    public Flux<CloudType, NumberFlux<CloudType>>
{
public:

    //- Runtime type information
    TypeName("numberFlux");


    // Constructors

        //- Inherit constructors
        using Flux<CloudType, NumberFlux<CloudType>>::Flux;


    // Static Member Data

        //- Dimensions of the flux field
        static const dimensionSet dimensions;


    // Member Functions

        //- Return the flux contribution for a single particle
        static scalar dPhiDeltaT(const typename CloudType::particleType& p)
        {
            return p.nParticle();
        }
};


/*---------------------------------------------------------------------------*\
                            Class VolumeFlux Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class VolumeFlux
:
    public Flux<CloudType, VolumeFlux<CloudType>>
{
public:

    //- Runtime type information
    TypeName("volumeFlux");


    // Constructors

        //- Inherit constructors
        using Flux<CloudType, VolumeFlux<CloudType>>::Flux;


    // Static Member Data

        //- Dimensions of the flux field
        static const dimensionSet dimensions;


    // Member Functions

        //- Return the flux contribution for a single particle
        static scalar dPhiDeltaT(const typename CloudType::particleType& p)
        {
            return p.nParticle()*p.volume();
        }
};


/*---------------------------------------------------------------------------*\
                            Class MassFlux Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class MassFlux
:
    public Flux<CloudType, MassFlux<CloudType>>
{
public:

    //- Runtime type information
    TypeName("massFlux");


    // Constructors

        //- Inherit constructors
        using Flux<CloudType, MassFlux<CloudType>>::Flux;


    // Static Member Data

        //- Dimensions of the flux field
        static const dimensionSet dimensions;


    // Member Functions

        //- Return the flux contribution for a single particle
        static scalar dPhiDeltaT(const typename CloudType::particleType& p)
        {
            return p.nParticle()*p.mass();
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "Flux.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
