/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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
    Foam::CloudFilmTransfer

Description
    Thermo parcel<->film transfer model.

    Responsible for:
    - injecting parcels from the film model into the cloud, e.g. for dripping
    - parcel interaction with the film, e.g absorb, bounce, splash

    References:
    \verbatim
        Bai, C., & Gosman, A. D. (1996).
        Mathematical modelling of wall films formed by impinging sprays.
        SAE transactions, 782-796.

        Bai, C. X., Rusche, H., & Gosman, A. D. (2002).
        Modeling of gasoline spray impingement.
        Atomization and Sprays, 12(1-3).
    \endverbatim

SourceFiles
    CloudFilmTransfer.C

\*---------------------------------------------------------------------------*/

#ifndef CloudFilmTransfer_H
#define CloudFilmTransfer_H

#include "SurfaceFilmModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class Random;

namespace fv
{
    class filmCloudTransfer;
}


/*---------------------------------------------------------------------------*\
                       Class CloudFilmTransferBase Declaration
\*---------------------------------------------------------------------------*/

class CloudFilmTransferBase
{
public:

    //- Interaction types
    enum class interactionType
    {
        absorb,
        bounce,
        splashBai
    };

    //- Interaction type names
    static const NamedEnum<interactionType, 3> interactionTypeNames_;
};


/*---------------------------------------------------------------------------*\
                      Class CloudFilmTransfer Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class CloudFilmTransfer
:
    public SurfaceFilmModel<CloudType>,
    public CloudFilmTransferBase
{
    // Private Member Functions

        //- Return pointers to the film transfer fvModels
        UPtrList<fv::filmCloudTransfer>& filmTransferPtrs() const;


protected:

    // Protected Data

        //- Convenience typedef to the cloud's parcel type
        typedef typename CloudType::parcelType parcelType;

        //- Reference to the cloud random number generator
        Random& rndGen_;


        // Film models

            //- Pointers to the films
            mutable UPtrList<fv::filmCloudTransfer> filmTransfers_;

            //- List of film patches
            mutable labelList filmPatches_;


        // Cached injector fields per film patch

            //- Film velocity / patch face
            vectorField UFilmPatch_;

            //- Film density / patch face
            scalarField rhoFilmPatch_;

            //- Film temperature / patch face
            scalarField TFilmPatch_;

            //- Film specific heat capacity / patch face
            scalarField CpFilmPatch_;


        // Interaction model data

            //- Interaction type enumeration
            interactionType interactionType_;

            //- Film thickness beyond which patch is assumed to be wet
            scalar deltaWet_;

            //- Splash parcel type label - id assigned to identify parcel for
            //  post-processing. If not specified, defaults to originating cloud
            //  type
            label splashParcelType_;

            //- Number of new parcels resulting from splash event
            label parcelsPerSplash_;


            // Surface roughness coefficient typically in the range 1300 - 5200
            // and decreases with increasing surface roughness

                //- Dry surface roughness coefficient
                //  = 2630 for dry interaction (ref. Bai)
                scalar Adry_;

                //- Wet surface roughness coefficient
                //  = 1320 for wet interaction (ref. Bai)
                scalar Awet_;


            //- Skin friction typically in the range 0.6 < Cf < 0.8
            scalar Cf_;

            //- Counter for number of new splash parcels
            label nParcelsSplashed_;


    // Protected Member Functions

        //- Return splashed parcel direction
        vector splashDirection
        (
            const vector& tanVec1,
            const vector& tanVec2,
            const vector& nf
        ) const;


        // Interaction models

            //- Absorb parcel into film
            void absorbInteraction
            (
                fv::filmCloudTransfer&,
                const parcelType& p,
                const polyPatch& pp,
                const label facei,
                const scalar mass,
                bool& keepParticle
            );

            //- Bounce parcel (flip parcel normal velocity)
            void bounceInteraction
            (
                parcelType& p,
                const polyPatch& pp,
                const label facei,
                bool& keepParticle
            ) const;

            //- Parcel interaction with dry surface
            void drySplashInteraction
            (
                fv::filmCloudTransfer&,
                const parcelType& p,
                const polyPatch& pp,
                const label facei,
                bool& keepParticle
            );

            //- Parcel interaction with wetted surface
            void wetSplashInteraction
            (
                fv::filmCloudTransfer&,
                parcelType& p,
                const polyPatch& pp,
                const label facei,
                bool& keepParticle
            );

            //- Bai parcel splash interaction model
            void splashInteraction
            (
                fv::filmCloudTransfer&,
                const parcelType& p,
                const polyPatch& pp,
                const label facei,
                const scalar mRatio,
                const scalar We,
                const scalar Wec,
                const scalar sigma,
                bool& keepParticle
            );

        //- Return pointers to the films
        virtual const labelList& filmPatches() const;

        //- Cache the film fields in preparation for injection
        virtual void cacheFilmFields(const label filmi);

        //- Set the individual parcel properties
        virtual void setParcelProperties
        (
            parcelType& p,
            const label filmFacei
        ) const;


public:

    //- Runtime type information
    TypeName("cloudFilmTransfer");


    // Constructors

        //- Construct from components
        CloudFilmTransfer(const dictionary& dict, CloudType& owner);

        //- Construct copy
        CloudFilmTransfer(const CloudFilmTransfer<CloudType>& sfm);

        //- Construct and return a clone using supplied owner cloud
        virtual autoPtr<SurfaceFilmModel<CloudType>> clone() const
        {
            return autoPtr<SurfaceFilmModel<CloudType>>
            (
                new CloudFilmTransfer<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~CloudFilmTransfer();


    // Member Functions

        // Evaluation

            //- Transfer parcel from cloud to film
            //  Returns true if parcel is to be transferred
            virtual bool transferParcel
            (
                parcelType& p,
                const polyPatch& pp,
                bool& keepParticle
            );


        // I-O

            //- Write film info to stream
            virtual void info(Ostream& os);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "CloudFilmTransfer.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
