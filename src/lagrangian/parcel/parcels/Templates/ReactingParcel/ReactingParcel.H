/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2025 OpenFOAM Foundation
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
    Foam::ReactingParcel

Description
    Reacting parcel class with one/two-way coupling with the continuous
    phase.

SourceFiles
    ReactingParcelI.H
    ReactingParcel.C
    ReactingParcelIO.C

\*---------------------------------------------------------------------------*/

#ifndef ReactingParcel_H
#define ReactingParcel_H

#include "particle.H"
#include "interpolation.H"
#include "fluidThermo.H"
#include "demandDrivenEntry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class ParcelType>
class ReactingParcel;

template<class ParcelType>
Ostream& operator<<
(
    Ostream&,
    const ReactingParcel<ParcelType>&
);


/*---------------------------------------------------------------------------*\
                     Class ReactingParcelName Declaration
\*---------------------------------------------------------------------------*/

TemplateName(ReactingParcel);


/*---------------------------------------------------------------------------*\
                       Class ReactingParcel Declaration
\*---------------------------------------------------------------------------*/

template<class ParcelType>
class ReactingParcel
:
    public ParcelType,
    public ReactingParcelName
{
    // Private Data

        //- Size in bytes of the fields
        static const std::size_t sizeofFields_;


public:

    //- Class to hold reacting parcel constant properties
    class constantProperties
    :
        public ParcelType::constantProperties
    {
        // Private Data

            //- Minimum pressure [Pa]
            demandDrivenEntry<scalar> pMin_;

            //- Constant volume flag - e.g. during mass transfer
            demandDrivenEntry<bool> constantVolume_;


    public:

        // Constructors

            //- Null constructor
            constantProperties();

            //- Copy constructor
            constantProperties(const constantProperties& cp);

            //- Construct from dictionary
            constantProperties(const dictionary& parentDict);


        // Access

            //- Return const access to the minimum pressure
            inline scalar pMin() const;

            //- Return const access to the constant volume flag
            inline bool constantVolume() const;
    };


    //- Use base tracking data
    typedef typename ParcelType::trackingData trackingData;


protected:

    // Protected data

        // Parcel properties

            //- Mass fractions of mixture []
            scalarField Y_;


    // Protected Member Functions

        //- Calculate Phase change
        template<class TrackCloudType>
        void calcPhaseChange
        (
            TrackCloudType& cloud,
            trackingData& td,
            const scalar dt,           // timestep
            const scalar Re,           // Reynolds number
            const scalar Pr,           // Prandtl number
            const scalar Ts,           // Surface temperature
            const scalar nus,          // Surface kinematic viscosity
            const scalar d,            // diameter
            const scalar T,            // temperature
            const scalar mass,         // mass
            const label idPhase,       // id of phase involved in phase change
            const scalar YPhase,       // total mass fraction
            const scalarField& YComponents, // component mass fractions
            scalarField& dMassPC,      // mass transfer - local to parcel
            scalar& Sh,                // explicit parcel enthalpy source
            scalar& N,                 // flux of species emitted from parcel
            scalar& NCpW,              // sum of N*Cp*W of emission species
            scalarField& Cs            // carrier conc. of emission species
        );

        //- Update mass fraction
        scalar updateMassFraction
        (
            const scalar mass0,
            const scalarField& dMass,
            scalarField& Y
        ) const;


public:

    // Static Data Members

        //- String representation of properties
        AddToPropertyList
        (
            ParcelType,
            " nPhases(Y1..YN)"
        );


    // Constructors

        //- Construct from mesh, coordinates and topology
        //  Other properties initialised as null
        inline ReactingParcel
        (
            const polyMesh& mesh,
            const barycentric& coordinates,
            const label celli,
            const label tetFacei,
            const label tetPti,
            const label facei
        );

        //- Construct from a position and a cell, searching for the rest of the
        //  required topology. Other properties are initialised as null.
        inline ReactingParcel
        (
            const meshSearch& searchEngine,
            const vector& position,
            const label celli,
            label& nLocateBoundaryHits
        );

        //- Construct from Istream
        ReactingParcel(Istream& is, bool readFields = true);

        //- Construct as a copy
        ReactingParcel(const ReactingParcel& p);

        //- Construct and return a clone
        virtual autoPtr<particle> clone() const
        {
            return autoPtr<particle>(new ReactingParcel<ParcelType>(*this));
        }

        //- Construct from Istream and return
        static autoPtr<ReactingParcel> New(Istream& is)
        {
            return autoPtr<ReactingParcel>(new ReactingParcel(is));
        }


    // Member Functions

        // Access

            //- Return const access to mass fractions of mixture []
            inline const scalarField& Y() const;


        // Edit

            //- Return access to mass fractions of mixture []
            inline scalarField& Y();


        // Main calculation loop

            //- Set cell values
            template<class TrackCloudType>
            void setCellValues(TrackCloudType& cloud, trackingData& td);

            //- Correct cell values using latest transfer information
            template<class TrackCloudType>
            void cellValueSourceCorrection
            (
                TrackCloudType& cloud,
                trackingData& td,
                const scalar dt
            );

            //- Correct surface values due to emitted species
            template<class TrackCloudType>
            void correctSurfaceValues
            (
                TrackCloudType& cloud,
                trackingData& td,
                const scalar T,
                const scalarField& Cs,
                scalar& rhos,
                scalar& mus,
                scalar& Prs,
                scalar& kappas
            );

            //- Update parcel properties over the time interval
            template<class TrackCloudType>
            void calc
            (
                TrackCloudType& cloud,
                trackingData& td,
                const scalar dt
            );


        // I-O

            //- Read
            template<class CloudType, class CompositionType>
            static void readFields
            (
                CloudType& c,
                const CompositionType& compModel
            );

            //- Read - no composition
            template<class CloudType>
            static void readFields(CloudType& c);

            //- Write
            template<class CloudType, class CompositionType>
            static void writeFields
            (
                const CloudType& c,
                const CompositionType& compModel
            );

            //- Write - composition supplied
            template<class CloudType>
            static void writeFields(const CloudType& c);


    // Ostream Operator

        friend Ostream& operator<< <ParcelType>
        (
            Ostream&,
            const ReactingParcel<ParcelType>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ReactingParcelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "ReactingParcel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
