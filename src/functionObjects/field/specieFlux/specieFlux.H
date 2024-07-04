/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024 OpenFOAM Foundation
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
    Foam::functionObjects::specieFlux

Description
    These functions calculate the specie-flux and write it as a
    surfaceScalarField called 'specie<Type>Flux(<specieName>)'. There are
    three such functions; specieAdvectiveFlux and specieDiffusiveFlux return
    the advective and diffusive parts of the flux, respectively, and
    specieFlux returns the total combined flux.

    Example of function object specification:
    \verbatim
    specieFlux
    {
        type    specieFlux; // specieAdvectiveFlux, specieDiffusiveFlux
        libs    ("libfieldFunctionObjects.so");
        field   NH3;
    }
    \endverbatim

    Or, using the standard configuration:
    \verbatim
    #includeFunc specieFlux(NH3)
    \endverbatim

Usage
    \table
        Property     | Description                | Required | Default value
        type         | Type name: specieFlux, \\
                       specieAdvectiveFlux, or \\
                       specieDiffusiveFlux        | yes      |
        field        | Name of the specie/mass \\
                       fraction field             | yes      |
        schemesField | Name of the field from \\
                       which schemes are taken    | no       | Yi
        region       | Region to be evaluated     | no       | default region
    \endtable

SourceFiles
    specieFlux.C

\*---------------------------------------------------------------------------*/

#ifndef specieFlux_H
#define specieFlux_H

#include "fieldExpression.H"
#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fluidThermophysicalTransportModel;

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                          Class specieFluxBase Declaration
\*---------------------------------------------------------------------------*/

class specieFluxBase
:
    public fieldExpression
{
    // Protected Data

        //- Name of field from which schemes are taken
        const word schemesField_;


    // Private Member Functions

        //- Calculate the specie flux field and return true if successful
        virtual bool calc();

        //- Calculate the specie flux field
        virtual tmp<surfaceScalarField> calc
        (
            const fluidThermophysicalTransportModel& ttm,
            const volScalarField& Yi
        ) = 0;


protected:

    // Protected Member Functions

        //- Return the advective flux
        tmp<surfaceScalarField> calcPhiYif
        (
            const fluidThermophysicalTransportModel& ttm,
            const volScalarField& Yi
        ) const;

        //- Return the diffusive flux
        tmp<surfaceScalarField> calcJ
        (
            const fluidThermophysicalTransportModel& ttm,
            const volScalarField& Yi
        ) const;

public:


    // Constructors

        //- Construct from Time and dictionary
        specieFluxBase
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict,
            const word& typeName
        );


    //- Destructor
    virtual ~specieFluxBase();
};


/*---------------------------------------------------------------------------*\
                       Class specieFlux Declaration
\*---------------------------------------------------------------------------*/

class specieFlux
:
    public specieFluxBase
{
private:

    // Private Member Functions

        //- Calculate the specie flux field
        virtual tmp<surfaceScalarField> calc
        (
            const fluidThermophysicalTransportModel& ttm,
            const volScalarField& Yi
        );


public:

    //- Runtime type information
    TypeName("specieFlux");


    // Constructors

        //- Construct from Time and dictionary
        specieFlux
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );
};


/*---------------------------------------------------------------------------*\
                       Class specieAdvectiveFlux Declaration
\*---------------------------------------------------------------------------*/

class specieAdvectiveFlux
:
    public specieFluxBase
{
private:

    // Private Member Functions

        //- Calculate the specie flux field
        virtual tmp<surfaceScalarField> calc
        (
            const fluidThermophysicalTransportModel& ttm,
            const volScalarField& Yi
        );


public:

    //- Runtime type information
    TypeName("specieAdvectiveFlux");


    // Constructors

        //- Construct from Time and dictionary
        specieAdvectiveFlux
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );
};


/*---------------------------------------------------------------------------*\
                       Class specieDiffusiveFlux Declaration
\*---------------------------------------------------------------------------*/

class specieDiffusiveFlux
:
    public specieFluxBase
{
private:

    // Private Member Functions

        //- Calculate the specie flux field
        virtual tmp<surfaceScalarField> calc
        (
            const fluidThermophysicalTransportModel& ttm,
            const volScalarField& Yi
        );


public:

    //- Runtime type information
    TypeName("specieDiffusiveFlux");


    // Constructors

        //- Construct from Time and dictionary
        specieDiffusiveFlux
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
