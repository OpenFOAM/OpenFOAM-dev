/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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
    Foam::carrierDensity

Description
    Class which provides the density of the carrier phase

SourceFiles
    carrierDensity.C

\*---------------------------------------------------------------------------*/

#ifndef carrierDensity_H
#define carrierDensity_H

#include "carrierProperties.H"
#include "carried.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class carrierDensity Declaration
\*---------------------------------------------------------------------------*/

class carrierDensity
:
    virtual public carrierProperties
{
protected:

    // Protected Static Member Functions

        //- Get the carrier density field
        static tmp<volScalarField> getRhoVf
        (
            const cloud& cloud,
            const word& phaseName,
            const bool allowEmpty
        );


public:

    // Public Classes

        //- Forward declare the implementation class
        class implementation;


    //- Runtime type information
    TypeName("carrierDensity");


    //- Select a carrier density
    static autoPtr<carrierDensity> New
    (
        const cloud& c,
        const clouds::carried& carriedCloud,
        const word& phaseName
    );


    //- Destructor
    virtual ~carrierDensity();


    // Member Functions

        //- Density [kg/m^3]
        virtual const CarrierField<scalar>& rho() const = 0;
};

/*---------------------------------------------------------------------------*\
             Class carrierDensity::implementation Declaration
\*---------------------------------------------------------------------------*/

class carrierDensity::implementation
:
    virtual public carrierDensity
{
    // Private Data

        //- Reference to the cloud
        const cloud& cloud_;

        //- The name of the phase
        const word& phaseName_;

        //- Cached carrier density field
        tmp<volScalarField> trhoVf_;

        //- Carrier density
        const CarrierField<scalar>& rho_;


public:

    // Constructors

        //- Construct from a reference to the cloud and its base classes
        implementation
        (
            const cloud& c,
            const clouds::carried& carriedCloud,
            const word& phaseName
        );


    //- Destructor
    virtual ~implementation();


    // Member Functions

        //- Correct
        virtual void correct();

        //- Density [kg/m^3]
        virtual const CarrierField<scalar>& rho() const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
