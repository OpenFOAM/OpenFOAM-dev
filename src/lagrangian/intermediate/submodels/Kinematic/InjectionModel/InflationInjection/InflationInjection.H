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
    Foam::InflationInjection

Description
    Inflation injection - creates new particles by splitting existing
    particles within in a set of generation cells, then inflating them
    to a target diameter within the generation cells and an additional
    set of inflation cells.

SourceFiles
    InflationInjection.C

\*---------------------------------------------------------------------------*/

#ifndef InflationInjection_H
#define InflationInjection_H

#include "InjectionModel.H"
#include "distributionModel.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Structure to hold:
// + position = vectorPairScalarPair::first().first()
// + velocity = vectorPairScalarPair::first().second()
// + diameter = vectorPairScalarPair::second().first()
// + target diameter = vectorPairScalarPair::second().second()
// One structure to allow single operation parallel comms
typedef Tuple2<Pair<vector>, Pair<scalar>> vectorPairScalarPair;


/*---------------------------------------------------------------------------*\
                      Class InflationInjection Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class InflationInjection
:
    public InjectionModel<CloudType>
{
    // Private data

        //- Name of cellSet for generating new particles
        word generationSetName_;

        //- Name of cellSet for inflating new particles
        word inflationSetName_;

        //- Set of cells to generate particles in
        labelList generationCells_;

        //- Set of cells to inflate particles in, includes all
        //  generation cells
        labelList inflationCells_;

        //- Injection duration [s]
        scalar duration_;

        //- Flow rate profile relative to SOI [m3/s]
        TimeFunction1<scalar> flowRateProfile_;

        //- Growth rate of particle diameters towards target [m/s]
        TimeFunction1<scalar> growthRate_;

        //- Positions, velocities, diameters and target diameters of
        //  new particles after splitting
        DynamicList<vectorPairScalarPair> newParticles_;

        //- Accumulation variable to carry over volume from one injection
        //  to the next
        scalar volumeAccumulator_;

        //- Fraction of injection controlled by this processor
        scalar fraction_;

        //- Switch to control whether or not the injector is allowed
        //  to create new particles in empty cells
        Switch selfSeed_;

        //- Diameter with which to create new seed particles
        scalar dSeed_;

        //- Parcel size distribution model
        const autoPtr<distributionModel> sizeDistribution_;


public:

    //- Runtime type information
    TypeName("inflationInjection");


    // Constructors

        //- Construct from dictionary
        InflationInjection
        (
            const dictionary& dict,
            CloudType& owner,
            const word& modelName
        );

        //- Construct copy
        InflationInjection(const InflationInjection<CloudType>& im);

        //- Construct and return a clone
        virtual autoPtr<InjectionModel<CloudType>> clone() const
        {
            return autoPtr<InjectionModel<CloudType>>
            (
                new InflationInjection<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~InflationInjection();


    // Member Functions

        //- Set injector locations when mesh is updated
        virtual void updateMesh();

        //- Return the end-of-injection time
        scalar timeEnd() const;

        //- Number of parcels to introduce relative to SOI
        virtual label parcelsToInject(const scalar time0, const scalar time1);

        //- Volume of parcels to introduce relative to SOI
        virtual scalar volumeToInject(const scalar time0, const scalar time1);


        // Injection geometry

            //- Set the injection position and owner cell, tetFace and tetPt
            virtual void setPositionAndCell
            (
                const label parcelI,
                const label nParcels,
                const scalar time,
                vector& position,
                label& cellOwner,
                label& tetFacei,
                label& tetPti
            );

            //- Set the parcel properties
            virtual void setProperties
            (
                const label parcelI,
                const label nParcels,
                const scalar time,
                typename CloudType::parcelType& parcel
            );

            //- Flag to identify whether model fully describes the parcel
            virtual bool fullyDescribed() const;

            //- Return flag to identify whether or not injection of parcelI is
            //  permitted
            virtual bool validInjection(const label parcelI);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "InflationInjection.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
