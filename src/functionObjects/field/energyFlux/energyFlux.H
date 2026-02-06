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
    Foam::functionObjects::energyFlux

Description
    These functions calculate the energy-flux and write it as a
    surfaceScalarField. There are three such functions; energyAdvectiveFlux and
    heatFlux return the advective and diffusive parts of the energy flux,
    respectively, and energyFlux returns the total combined flux.

    Example of function object specification:
    \verbatim
    energyFlux
    {
        type    energyFlux; // energyAdvectiveFlux, heatFlux
        libs    ("libfieldFunctionObjects.so");
    }
    \endverbatim

    Or, using the standard configuration:
    \verbatim
    #includeFunc energyFlux
    \endverbatim

Usage
    \table
        Property     | Description                | Required | Default value
        type         | Type name: energyFlux, \\
                       energyAdvectiveFlux, or \\
                       heatFlux                   | yes      |
        phase        | Name of the phase          | no       |
        region       | Region to be evaluated     | no       | default region
    \endtable

SourceFiles
    energyFlux.C

\*---------------------------------------------------------------------------*/

#ifndef energyFlux_functionObject_H
#define energyFlux_functionObject_H

#include "fieldExpression.H"
#include "surfaceFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class thermophysicalTransportModel;

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                         Class energyFluxBase Declaration
\*---------------------------------------------------------------------------*/

class energyFluxBase
:
    public fieldExpression
{
    // Protected Data

        //- Name of field from which schemes are taken
        const word schemesField_;

        //- The name of the phase
        word phaseName_;


    // Private Member Functions

        //- Calculate the specie flux field
        virtual bool calc();

        //- Calculate the specie flux field
        virtual tmp<surfaceScalarField> calc
        (
            const thermophysicalTransportModel& ttm
        ) const = 0;


protected:

    // Protected Member Functions

        //- Return the advective energy flux
        tmp<surfaceScalarField> calcPhihs
        (
            const thermophysicalTransportModel& ttm
        ) const;

        //- Return the heat flux
        tmp<surfaceScalarField> calcQ
        (
            const thermophysicalTransportModel& ttm
        ) const;


public:

    // Constructors

        //- Construct from Time and dictionary
        energyFluxBase
        (
            const word& name,
            const Time& runTime,
            const dictionary&,
            const word& typeName
        );


    //- Destructor
    virtual ~energyFluxBase();


    // Member Functions

        //- Read the data
        virtual bool read(const dictionary&);
};


/*---------------------------------------------------------------------------*\
                          Class energyFlux Declaration
\*---------------------------------------------------------------------------*/

class energyFlux
:
    public energyFluxBase
{
private:

    // Private Member Functions

        //- Calculate the energy flux field
        virtual tmp<surfaceScalarField> calc
        (
            const thermophysicalTransportModel& ttm
        ) const;


public:

    //- Runtime type information
    TypeName("energyFlux");


    // Constructors

        //- Construct from Time and dictionary
        energyFlux
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );
};


/*---------------------------------------------------------------------------*\
                      Class energyAdvectiveFlux Declaration
\*---------------------------------------------------------------------------*/

class energyAdvectiveFlux
:
    public energyFluxBase
{
private:

    // Private Member Functions

        //- Calculate the energy flux field
        virtual tmp<surfaceScalarField> calc
        (
            const thermophysicalTransportModel& ttm
        ) const;


public:

    //- Runtime type information
    TypeName("energyAdvectiveFlux");


    // Constructors

        //- Construct from Time and dictionary
        energyAdvectiveFlux
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );
};


/*---------------------------------------------------------------------------*\
                          Class heatFlux Declaration
\*---------------------------------------------------------------------------*/

class heatFlux
:
    public energyFluxBase
{
private:

    // Private Member Functions

        //- Calculate the heat flux field
        virtual tmp<surfaceScalarField> calc
        (
            const thermophysicalTransportModel& ttm
        ) const;


public:

    //- Runtime type information
    TypeName("heatFlux");


    // Constructors

        //- Construct from Time and dictionary
        heatFlux
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
