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
    Foam::ThermoParcel

Description
    Thermodynamic parcel class with one/two-way coupling with the continuous
    phase

SourceFiles
    ThermoParcelI.H
    ThermoParcel.C
    ThermoParcelIO.C

\*---------------------------------------------------------------------------*/

#ifndef ThermoParcel_H
#define ThermoParcel_H

#include "particle.H"
#include "interpolation.H"
#include "fluidThermo.H"
#include "demandDrivenEntry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class ParcelType>
class ThermoParcel;

template<class ParcelType>
Ostream& operator<<
(
    Ostream&,
    const ThermoParcel<ParcelType>&
);


/*---------------------------------------------------------------------------*\
                      Class ThermoParcelName Declaration
\*---------------------------------------------------------------------------*/

TemplateName(ThermoParcel);


/*---------------------------------------------------------------------------*\
                        Class ThermoParcel Declaration
\*---------------------------------------------------------------------------*/

template<class ParcelType>
class ThermoParcel
:
    public ParcelType,
    public ThermoParcelName
{
    // Private Data

        //- Size in bytes of the fields
        static const std::size_t sizeofFields_;


public:

    //- Class to hold thermo particle constant properties
    class constantProperties
    :
        public ParcelType::constantProperties
    {

        // Private Data

            //- Particle initial temperature [K]
            demandDrivenEntry<scalar> T0_;

            //- Minimum temperature [K]
            demandDrivenEntry<scalar> TMin_;

            //- Maximum temperature [K]
            demandDrivenEntry<scalar> TMax_;

            //- Particle specific heat capacity [J/kg/K]
            demandDrivenEntry<scalar> Cp0_;

            //- Particle emissivity [] (radiation)
            demandDrivenEntry<scalar> epsilon0_;

            //- Particle scattering factor [] (radiation)
            demandDrivenEntry<scalar> f0_;


    public:

        // Constructors

            //- Null constructor
            constantProperties();

            //- Copy constructor
            constantProperties(const constantProperties& cp);

            //- Construct from dictionary
            constantProperties(const dictionary& parentDict);


        // Member Functions

            // Access

                //- Return const access to the particle initial temperature [K]
                inline scalar T0() const;

                //- Return const access to minimum temperature [K]
                inline scalar TMin() const;

                //- Return const access to maximum temperature [K]
                inline scalar TMax() const;

                //- Set the maximum temperature [K]
                inline void setTMax(const scalar TMax);

                //- Return const access to the particle specific heat capacity
                //  [J/kg/K]
                inline scalar Cp0() const;

                //- Return const access to the particle emissivity []
                //  Active for radiation only
                inline scalar epsilon0() const;

                //- Return const access to the particle scattering factor []
                //  Active for radiation only
                inline scalar f0() const;
    };


    class trackingData
    :
        public ParcelType::trackingData
    {
    private:

        // Private Data

            //- Local copy of carrier specific heat field
            //  Cp not stored on carrier thermo, but returned as tmp<...>
            const volScalarField Cp_;

            //- Local copy of carrier thermal conductivity field
            //  kappa not stored on carrier thermo, but returned as tmp<...>
            const volScalarField kappa_;


            // Interpolators for continuous phase fields

                //- Interpolator for continuous phase pressure field
                autoPtr<interpolation<scalar>> pInterp_;

                //- Temperature field interpolator
                autoPtr<interpolation<scalar>> TInterp_;

                //- Specific heat capacity field interpolator
                autoPtr<interpolation<scalar>> CpInterp_;

                //- Thermal conductivity field interpolator
                autoPtr<interpolation<scalar>> kappaInterp_;

                //- Radiation field interpolator
                autoPtr<interpolation<scalar>> GInterp_;


            // Cached continuous phase properties

                //- Pressure [Pa]
                scalar pc_;

                //- Temperature [K]
                scalar Tc_;

                //- Specific heat capacity [J/kg/K]
                scalar Cpc_;


    public:

        // Constructors

            //- Construct from components
            template <class TrackCloudType>
            inline trackingData(const TrackCloudType& cloud);


        // Member Functions

            //- Return access to the locally stored carrier Cp field
            inline const volScalarField& Cp() const;

            //- Return access to the locally stored carrier kappa field
            inline const volScalarField& kappa() const;

            //- Return const access to the interpolator for continuous phase
            //  pressure field
            inline const interpolation<scalar>& pInterp() const;

            //- Return const access to the interpolator for continuous
            //  phase temperature field
            inline const interpolation<scalar>& TInterp() const;

            //- Return const access to the interpolator for continuous
            //  phase specific heat capacity field
            inline const interpolation<scalar>& CpInterp() const;

            //- Return const access to the interpolator for continuous
            //  phase thermal conductivity field
            inline const interpolation<scalar>& kappaInterp() const;

            //- Return const access to the interpolator for continuous
            //  radiation field
            inline const interpolation<scalar>& GInterp() const;

            //- Return the continuous phase pressure
            inline scalar pc() const;

            //- Access the continuous phase pressure
            inline scalar& pc();

            //- Return the continuous phase temperature
            inline scalar Tc() const;

            //- Access the continuous phase temperature
            inline scalar& Tc();

            //- Return the continuous phase specific heat capacity
            inline scalar Cpc() const;

            //- Access the continuous phase specific heat capacity
            inline scalar& Cpc();
    };


protected:

    // Protected data

        // Parcel properties

            //- Temperature [K]
            scalar T_;

            //- Specific heat capacity [J/kg/K]
            scalar Cp_;


    // Protected Member Functions

        //- Calculate new particle temperature
        template<class TrackCloudType>
        scalar calcHeatTransfer
        (
            TrackCloudType& cloud,
            trackingData& td,
            const scalar dt,           // timestep
            const scalar Re,           // Reynolds number
            const scalar Pr,           // Prandtl number - surface
            const scalar kappa,        // Thermal conductivity - surface
            const scalar NCpW,         // Sum of N*Cp*W of emission species
            const scalar Sh,           // explicit particle enthalpy source
            scalar& dhsTrans,          // sensible enthalpy transfer to carrier
            scalar& Sph                // linearised heat transfer coefficient
        );


public:

    // Static Data Members

        //- String representation of properties
        AddToPropertyList
        (
            ParcelType,
            " T"
          + " Cp"
        );


    // Constructors

        //- Construct from mesh, coordinates and topology
        //  Other properties initialised as null
        inline ThermoParcel
        (
            const polyMesh& mesh,
            const barycentric& coordinates,
            const label celli,
            const label tetFacei,
            const label tetPti
        );

        //- Construct from a position and a cell, searching for the rest of the
        //  required topology. Other properties are initialised as null.
        inline ThermoParcel
        (
            const polyMesh& mesh,
            const vector& position,
            const label celli
        );

        //- Construct from Istream
        ThermoParcel
        (
            const polyMesh& mesh,
            Istream& is,
            bool readFields = true
        );

        //- Construct as a copy
        ThermoParcel(const ThermoParcel& p);

        //- Construct as a copy
        ThermoParcel(const ThermoParcel& p, const polyMesh& mesh);

        //- Construct and return a (basic particle) clone
        virtual autoPtr<particle> clone() const
        {
            return autoPtr<particle>(new ThermoParcel(*this));
        }

        //- Construct and return a (basic particle) clone
        virtual autoPtr<particle> clone(const polyMesh& mesh) const
        {
            return autoPtr<particle>(new ThermoParcel(*this, mesh));
        }

        //- Factory class to read-construct particles used for
        //  parallel transfer
        class iNew
        {
            const polyMesh& mesh_;

        public:

            iNew(const polyMesh& mesh)
            :
                mesh_(mesh)
            {}

            autoPtr<ThermoParcel<ParcelType>> operator()(Istream& is) const
            {
                return autoPtr<ThermoParcel<ParcelType>>
                (
                    new ThermoParcel<ParcelType>(mesh_, is, true)
                );
            }
        };


    // Member Functions

        // Access

            //- Return const access to temperature
            inline scalar T() const;

            //- Return const access to specific heat capacity
            inline scalar Cp() const;


        // Edit

            //- Return access to temperature
            inline scalar& T();

            //- Return access to specific heat capacity
            inline scalar& Cp();


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

            //- Calculate surface thermo properties
            template<class TrackCloudType>
            void calcSurfaceValues
            (
                TrackCloudType& cloud,
                trackingData& td,
                const scalar T,
                scalar& Ts,
                scalar& rhos,
                scalar& mus,
                scalar& Pr,
                scalar& kappas
            ) const;

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
            template<class CloudType>
            static void readFields(CloudType& c);

            //- Write
            template<class CloudType>
            static void writeFields(const CloudType& c);

            //- Write
            template<class CloudType, class CompositionType>
            static void writeFields
            (
                const CloudType& c,
                const CompositionType& compModel
            );


    // Ostream Operator

        friend Ostream& operator<< <ParcelType>
        (
            Ostream&,
            const ThermoParcel<ParcelType>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ThermoParcelI.H"
#include "ThermoParcelTrackingDataI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "ThermoParcel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
