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
    Foam::carrierViscosity

Description
    Class which provides the kinematic viscosity of the carrier phase

SourceFiles
    carrierViscosity.C

\*---------------------------------------------------------------------------*/

#ifndef carrierViscosity_H
#define carrierViscosity_H

#include "carrierProperties.H"
#include "carried.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class carrierViscosity Declaration
\*---------------------------------------------------------------------------*/

class carrierViscosity
:
    virtual public carrierProperties
{
protected:

    // Protected Static Member Functions

        //- Get the carrier kinematic viscosity field
        static tmp<volScalarField> getNuVf
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
    TypeName("carrierViscosity");


    //- Select a carrier viscosity
    static autoPtr<carrierViscosity> New
    (
        const cloud& c,
        const clouds::carried& carriedCloud,
        const word& phaseName
    );


    //- Destructor
    virtual ~carrierViscosity();


    // Member Functions

        //- Kinematic viscosity [m^2/s]
        virtual const CloudDerivedField<scalar>& nu() const = 0;
};

/*---------------------------------------------------------------------------*\
             Class carrierViscosity::implementation Declaration
\*---------------------------------------------------------------------------*/

class carrierViscosity::implementation
:
    virtual public carrierViscosity
{
    // Private Data

        //- Reference to the cloud
        const cloud& cloud_;

        //- The name of the phase
        const word& phaseName_;

        //- Cached carrier kinematic viscosity field
        tmp<volScalarField> tnuVf_;

        //- Carrier kinematic viscosity
        const CloudDerivedField<scalar>& nu_;


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

        //- Kinematic viscosity [m^2/s]
        virtual const CloudDerivedField<scalar>& nu() const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
