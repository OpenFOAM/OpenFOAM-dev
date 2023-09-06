/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020-2023 OpenFOAM Foundation
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
    Foam::parcelCloudList

Description
    List of parcel clouds, with the same interface as an individual parcel
    cloud. This is the object that should be constructed by an fvModel, or any
    system that can call this class' mesh change functions. A solver should
    *not* construct this object, as that would not provide a mechanism for the
    mesh change functions to be executed. A solver should construct a
    parcelClouds object instead.

SourceFiles
    parcelCloudList.C

\*---------------------------------------------------------------------------*/

#ifndef parcelCloudList_H
#define parcelCloudList_H

#include "parcelCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class parcelCloudList
:
    public PtrList<parcelCloud>
{
private:

    // Private Static Member Functions

        //- Get the cloud names for this case
        static wordList cloudNames(const objectRegistry& db);


    // Private data

        //- Reference to the mesh
        const fvMesh& mesh_;


public:

    // Static Data Members

        //- The default cloud name
        static const word defaultCloudName;

        //- The default cloud names (i.e., a list of length one with the
        //  defaultCloudName as its value)
        static const wordList defaultCloudNames;

        //- The name of the clouds file in which multiple cloud names are
        //  specified
        static const word cloudNamesName;


    // Constructors

        //- Construct specified clouds with given carrier fields
        parcelCloudList
        (
            const wordList& cloudNames,
            const volScalarField& rho,
            const volVectorField& U,
            const volScalarField& mu,
            const dimensionedVector& g
        );

        //- Construct specified clouds with given carrier fields and thermo
        parcelCloudList
        (
            const wordList& cloudNames,
            const volScalarField& rho,
            const volVectorField& U,
            const dimensionedVector& g,
            const fluidThermo& carrierThermo
        );

        //- Construct detected clouds with given carrier fields
        parcelCloudList
        (
            const volScalarField& rho,
            const volVectorField& U,
            const volScalarField& mu,
            const dimensionedVector& g
        );

        //- Construct detected clouds with given carrier fields and thermo
        parcelCloudList
        (
            const volScalarField& rho,
            const volVectorField& U,
            const dimensionedVector& g,
            const fluidThermo& carrierThermo
        );

        //- Disallow default bitwise copy construction
        parcelCloudList(const parcelCloudList&) = delete;


    //- Destructor
    ~parcelCloudList();


    // Member Functions

        // Fields

            //- Return the particle volume fraction field
            const tmp<volScalarField> alpha() const;


        // Sources

            // Momentum

                //- Return momentum source term [kg m/s^2]
                tmp<fvVectorMatrix> SU(const volVectorField& U) const;

                //- Momentum transfer [kg m/s]
                tmp<volVectorField::Internal> UTrans() const;

                //- Momentum transfer coefficient [kg]
                tmp<volScalarField::Internal> UCoeff() const;


            // Energy

                //- Return sensible enthalpy source term [J/s]
                tmp<fvScalarMatrix> Sh(const volScalarField& hs) const;

                //- Sensible enthalpy transfer [J]
                tmp<volScalarField::Internal> hsTrans() const;

                //- Sensible enthalpy transfer coefficient [J/K]
                tmp<volScalarField::Internal> hsCoeff() const;

                //- Return equivalent particulate emission [kg/m/s^3]
                tmp<volScalarField> Ep() const;

                //- Return equivalent particulate absorption [1/m]
                tmp<volScalarField> ap() const;

                //- Return equivalent particulate scattering factor [1/m]
                tmp<volScalarField> sigmap() const;


            // Mass

                //- Return mass source term for specie [kg/s]
                tmp<fvScalarMatrix> SYi
                (
                    const label speciei,
                    const volScalarField& Yi
                ) const;

                //- Return total mass source term [kg/s]
                tmp<fvScalarMatrix> Srho(const volScalarField& rho) const;

                //- Return total mass source [kg/m^3/s]
                tmp<volScalarField::Internal> Srho() const;


        // I-O

            //- Print cloud information
            void info();


        // Evolution

            //- Evolve the cloud
            void evolve();


        // Mesh changes

            //- Call this before a topology change. Stores the particles global
            //  positions in the database for use during mapping.
            void storeGlobalPositions();

            //- Update topology using the given map
            void topoChange(const polyTopoChangeMap&);

            //- Update from another mesh using the given map
            void mapMesh(const polyMeshMap&);

            //- Redistribute or update using the given distribution map
            void distribute(const polyDistributionMap&);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const parcelCloudList&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
