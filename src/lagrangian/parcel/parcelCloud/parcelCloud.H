/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020 OpenFOAM Foundation
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
    Foam::parcelCloud

Description
    Virtual abstract base class for parcel clouds. As parcelCloudBase but with
    additional virtualisation of the evolve method, plus some additional
    methods that are defined below the parcel-cloud layer (i.e., in Cloud).
    These methods are implemented by forwarding in the ParcelCloud class.

    This class contains the selection table and represents the high level
    interface used by a solver.

SourceFiles
    parcelCloud.C

\*---------------------------------------------------------------------------*/

#ifndef parcelCloud_H
#define parcelCloud_H

#include "parcelCloudBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class SLGThermo;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class parcelCloud
:
    virtual public parcelCloudBase
{
public:

    //- Runtime type information
    TypeName("parcelCloud");


    // Declare viscosity-based runtime construction
    declareRunTimeSelectionTable
    (
        autoPtr,
        parcelCloud,
        viscosity,
        (
            const word& name,
            const volScalarField& rho,
            const volVectorField& U,
            const volScalarField& mu,
            const dimensionedVector& g
        ),
        (name, rho, U, mu, g)
    );

    // Declare thermo-based runtime construction
    declareRunTimeSelectionTable
    (
        autoPtr,
        parcelCloud,
        thermo,
        (
            const word& name,
            const volScalarField& rho,
            const volVectorField& U,
            const dimensionedVector& g,
            const SLGThermo& thermo
        ),
        (name, rho, U, g, thermo)
    );


    // Constructors

        //- Null constructor
        parcelCloud();

        //- Disallow default bitwise copy construction
        parcelCloud(const parcelCloud&) = delete;


    //- Selectors

        //- Select and construct with given carrier fields
        static autoPtr<parcelCloud> New
        (
            const word& name,
            const volScalarField& rho,
            const volVectorField& U,
            const volScalarField& mu,
            const dimensionedVector& g
        );

        //- Select and construct with given carrier fields and thermo
        static autoPtr<parcelCloud> New
        (
            const word& name,
            const volScalarField& rho,
            const volVectorField& U,
            const dimensionedVector& g,
            const SLGThermo& thermo
        );


    //- Destructor
    virtual ~parcelCloud();


    // Member Functions

        //- Call this before a topology change. Stores the particles global
        //  positions in the database for use during mapping.
        virtual void storeGlobalPositions() = 0;


        // Evolution

            //- Evolve the cloud
            virtual void evolve() = 0;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const parcelCloud&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
