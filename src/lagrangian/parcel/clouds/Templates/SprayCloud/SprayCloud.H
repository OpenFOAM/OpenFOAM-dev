/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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
    Foam::SprayCloud

Description
    Templated base class for spray cloud

    - sub-models:
      - atomisation model
      - break-up model

\*---------------------------------------------------------------------------*/

#ifndef SprayCloud_H
#define SprayCloud_H

#include "volFieldsFwd.H"
#include "fvMatricesFwd.H"
#include "dimensionedTypes.H"
#include "fvMesh.H"
#include "fluidThermo.H"
#include "Cloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
template<class CloudType>
class AtomisationModel;

template<class CloudType>
class BreakupModel;

/*---------------------------------------------------------------------------*\
                       Class SprayCloudName Declaration
\*---------------------------------------------------------------------------*/

TemplateName(SprayCloud);


/*---------------------------------------------------------------------------*\
                         Class SprayCloud Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class SprayCloud
:
    public CloudType,
    public SprayCloudName
{
public:

    // Public Typedefs

        //- Type of cloud this cloud was instantiated for
        typedef CloudType cloudType;

        //- Type of parcel the cloud was instantiated for
        typedef typename CloudType::particleType parcelType;

        //- Convenience typedef for this cloud type
        typedef SprayCloud<CloudType> sprayCloudType;


private:

    // Private Data

        //- Cloud copy pointer
        autoPtr<SprayCloud<CloudType>> cloudCopyPtr_;

        //- Average parcel mass
        scalar averageParcelMass_;


protected:

    // Protected data

        // References to the cloud sub-models

            //- Atomisation model
            autoPtr<AtomisationModel<SprayCloud<CloudType>>>
                atomisationModel_;

            //- Break-up model
            autoPtr<BreakupModel<SprayCloud<CloudType>>> breakupModel_;


    // Protected Member Functions

        // Initialisation

            //- Set cloud sub-models
            void setModels();


        // Cloud evolution functions

            //- Reset state of cloud
            void cloudReset(SprayCloud<CloudType>& c);


public:

    // Constructors

        //- Construct given carrier gas fields
        SprayCloud
        (
            const word& cloudName,
            const volScalarField& rho,
            const volVectorField& U,
            const dimensionedVector& g,
            const fluidThermo& carrierThermo,
            bool readFields = true
        );

        //- Copy constructor with new name
        SprayCloud(SprayCloud<CloudType>& c, const word& name);

        //- Copy constructor with new name - creates bare cloud
        SprayCloud
        (
            const fvMesh& mesh,
            const word& name,
            const SprayCloud<CloudType>& c
        );

        //- Disallow default bitwise copy construction
        SprayCloud(const SprayCloud&) = delete;

        //- Construct and return clone based on (this) with new name
        virtual autoPtr<Cloud<parcelType>> clone(const word& name)
        {
            return autoPtr<Cloud<parcelType>>
            (
                new SprayCloud(*this, name)
            );
        }

        //- Construct and return bare clone based on (this) with new name
        virtual autoPtr<Cloud<parcelType>> cloneBare(const word& name) const
        {
            return autoPtr<Cloud<parcelType>>
            (
                new SprayCloud(this->mesh(), name, *this)
            );
        }


    //- Destructor
    virtual ~SprayCloud();


    // Member Functions

        // Access

            //- Return a reference to the cloud copy
            inline const SprayCloud& cloudCopy() const;

            //- Return const-access to the average parcel mass
            inline scalar averageParcelMass() const;


        // Check

            //- Penetration for fraction [0-1] of the current total mass
            inline scalar penetration(const scalar fraction) const;


            // Sub-models

                //- Return const-access to the atomisation model
                inline const AtomisationModel<SprayCloud<CloudType>>&
                    atomisation() const;

                //- Return reference to the atomisation model
                inline AtomisationModel<SprayCloud<CloudType>>& atomisation();

                //- Return const-access to the breakup model
                inline const BreakupModel<SprayCloud<CloudType>>&
                    breakup() const;

                //- Return reference to the breakup model
                inline BreakupModel<SprayCloud<CloudType>>& breakup();


        // Cloud evolution functions

            //- Set parcel thermo properties
            void setParcelThermoProperties
            (
                parcelType& parcel,
                const scalar lagrangianDt
            );

            //- Check parcel properties
            void checkParcelProperties
            (
                parcelType& parcel,
                const scalar lagrangianDt,
                const bool fullyDescribed
            );

            //- Store the current cloud state
            void storeState();

            //- Reset the current cloud to the previously stored state
            void restoreState();

            //- Evolve the spray (inject, move)
            void evolve();


        // I-O

            //- Print cloud information
            void info();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const SprayCloud&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "SprayCloudI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "SprayCloud.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
