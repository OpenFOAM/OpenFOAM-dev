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
    Foam::carrierDynamicViscosity

Description
    Class which provides the dynamic viscosity of the carrier phase

SourceFiles
    carrierDynamicViscosity.C

\*---------------------------------------------------------------------------*/

#ifndef carrierDynamicViscosity_H
#define carrierDynamicViscosity_H

#include "carrierViscosity.H"
#include "carrierDensity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class carrierDynamicViscosity Declaration
\*---------------------------------------------------------------------------*/

class carrierDynamicViscosity
:
    virtual public carrierViscosity,
    virtual public carrierDensity
{
protected:

    // Protected Static Member Functions

        //- Get the carrier dynamic viscosity field
        static const volScalarField& getMuVf
        (
            const cloud& cloud,
            const word& phaseName,
            const bool allowNull
        );


public:

    // Public Classes

        //- Forward declare the implementation class
        class implementation;


    //- Runtime type information
    TypeName("carrierDynamicViscosity");


    //- Select a carrier dynamic viscosity
    static autoPtr<carrierDynamicViscosity> New
    (
        const cloud& c,
        const clouds::carried& carriedCloud,
        const word& phaseName
    );


    //- Destructor
    virtual ~carrierDynamicViscosity();


    // Member Functions

        //- Dynamic viscosity [kg/m/s]
        virtual const CloudDerivedField<scalar>& mu() const = 0;
};

/*---------------------------------------------------------------------------*\
             Class carrierDynamicViscosity::implementation Declaration
\*---------------------------------------------------------------------------*/

class carrierDynamicViscosity::implementation
:
    virtual public carrierDynamicViscosity,
    public carrierDensity::implementation
{
    // Private Data

        //- Reference to the cloud
        const cloud& cloud_;

        //- The name of the phase
        const word& phaseName_;

        //- Carrier kinematic viscosity
        tmp<volScalarField> tnuVf_;

        //- Carrier dynamic viscosity
        const volScalarField& muVf_;

        //- Carrier kinematic viscosity
        const CloudDerivedField<scalar>& nu_;

        //- Carrier dynamic viscosity
        const CloudDerivedField<scalar>& mu_;


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

        //- Dynamic viscosity [kg/m/s]
        virtual const CloudDerivedField<scalar>& mu() const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
