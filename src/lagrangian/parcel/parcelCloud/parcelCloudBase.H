/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020-2021 OpenFOAM Foundation
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
    Foam::parcelCloudBase

Description
    Virtual abstract base class for parcel clouds. Inserted by ParcelCloudBase
    into the base of the cloud template hierarchy and adds virtualisation of
    most methods defined by the clouds.

    Note that the "evolve" method is not virtualised here. Due to the way in
    which TrackCloudType and trackingData templating work, it is not possible
    to virtualise this method directly. Instead it has to be wrapped. That is
    achieved by the parcelCloud and ParcelCloud classes.

SourceFiles
    parcelCloudBase.C

\*---------------------------------------------------------------------------*/

#ifndef parcelCloudBase_H
#define parcelCloudBase_H

#include "volFields.H"
#include "fvMatricesFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class parcelCloudBase
{
public:

    //- Runtime type information
    TypeName("parcelCloudBase");


    // Constructors

        //- Null constructor
        parcelCloudBase();

        //- Disallow default bitwise copy construction
        parcelCloudBase(const parcelCloudBase&) = delete;


    //- Destructor
    virtual ~parcelCloudBase();


    // Member Functions

        // Check

            //-  Number of parcels
            virtual label nParcels() const = 0;

            //- Total mass in system
            virtual scalar massInSystem() const = 0;


        // Fields

            //- Volume swept rate of parcels per cell
            virtual const tmp<volScalarField> vDotSweep() const = 0;

            //- Return the particle volume fraction field
            virtual const tmp<volScalarField> theta() const = 0;

            //- Return the particle mass fraction field
            virtual const tmp<volScalarField> alpha() const = 0;

            //- Return the particle effective density field
            virtual const tmp<volScalarField> rhoEff() const = 0;


        // Sources

            // Momentum

                //- Return momentum source term [kg m/s^2]
                virtual tmp<fvVectorMatrix> SU
                (
                    const volVectorField& U
                ) const = 0;

                //- Momentum transfer [kg m/s]
                virtual tmp<volVectorField::Internal> UTrans() const = 0;

                //- Momentum transfer coefficient [kg]
                virtual tmp<volScalarField::Internal> UCoeff() const = 0;


            // Energy

                //- Return sensible enthalpy source term [J/s]
                virtual tmp<fvScalarMatrix> Sh
                (
                    const volScalarField& hs
                ) const = 0;

                //- Sensible enthalpy transfer [J]
                virtual tmp<volScalarField::Internal> hsTrans() const = 0;

                //- Sensible enthalpy transfer coefficient [J/K]
                virtual tmp<volScalarField::Internal> hsCoeff() const = 0;

                //- Return equivalent particulate emission [kg/m/s^3]
                virtual tmp<volScalarField> Ep() const = 0;

                //- Return equivalent particulate absorption [1/m]
                virtual tmp<volScalarField> ap() const = 0;

                //- Return equivalent particulate scattering factor [1/m]
                virtual tmp<volScalarField> sigmap() const = 0;


            // Mass

                //- Return mass source term for specie [kg/s]
                virtual tmp<fvScalarMatrix> SYi
                (
                    const label i,
                    const volScalarField& Yi
                ) const = 0;

                //- Return mass source for specie [kg/m^3/s]
                //  Why is this not called SYi ???
                virtual tmp<volScalarField::Internal> Srho
                (
                    const label i
                ) const = 0;

                //- Return total mass source term [kg/s]
                virtual tmp<fvScalarMatrix> Srho(volScalarField& rho) const = 0;

                //- Return total mass source [kg/m^3/s]
                virtual tmp<volScalarField::Internal> Srho() const = 0;


        // I-O

            //- Print cloud information
            virtual void info() = 0;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const parcelCloudBase&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
