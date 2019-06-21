/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2019 OpenFOAM Foundation
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
    Foam::regionModels::surfaceFilmModels::liquidFilmThermo

Description
    Liquid thermo model

SourceFiles
    liquidFilmThermo.C

\*---------------------------------------------------------------------------*/

#ifndef liquidFilmThermo_H
#define liquidFilmThermo_H

#include "filmThermoModel.H"
#include "liquidProperties.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

// Forward class declarations
class thermoSingleLayer;

/*---------------------------------------------------------------------------*\
                      Class liquidFilmThermo Declaration
\*---------------------------------------------------------------------------*/

class liquidFilmThermo
:
    public filmThermoModel
{
protected:

    // Protected data

        //- Liquid name
        word name_;

        //- Pointer to the liquid properties
        const liquidProperties* liquidPtr_;

        //- Flag to indicate that model owns the liquid object
        bool ownLiquid_;

        //- Flag to indicate that reference values of p and T should be used
        bool useReferenceValues_;

        //- Reference pressure [pa]
        scalar pRef_;

        //- Reference temperature [K]
        scalar TRef_;


    // Protected member functions

        //- Return a reference to a thermo film
        const thermoSingleLayer& thermoFilm() const;

        //- Initialise the liquid pointer
        void initLiquid(const dictionary& dict);


public:

    //- Runtime type information
    TypeName("liquid");


    // Constructors

        //- Construct from surface film model and dictionary
        liquidFilmThermo
        (
            surfaceFilmRegionModel& film,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        liquidFilmThermo(const liquidFilmThermo&) = delete;


    //- Destructor
    virtual ~liquidFilmThermo();


    // Member Functions

        //- Return the specie name
        virtual const word& name() const;


        // Elemental access

            //- Return density [kg/m^3]
            virtual scalar rho(const scalar p, const scalar T) const;

            //- Return dynamic viscosity [Pa.s]
            virtual scalar mu(const scalar p, const scalar T) const;

            //- Return surface tension [kg/s^2]
            virtual scalar sigma(const scalar p, const scalar T) const;

            //- Return specific heat capacity [J/kg/K]
            virtual scalar Cp(const scalar p, const scalar T) const;

            //- Return thermal conductivity [W/m/K]
            virtual scalar kappa(const scalar p, const scalar T) const;

            //- Return diffusivity [m2/s]
            virtual scalar D(const scalar p, const scalar T) const;

            //- Return latent heat [J/kg]
            virtual scalar hl(const scalar p, const scalar T) const;

            //- Return vapour pressure [Pa]
            virtual scalar pv(const scalar p, const scalar T) const;

            //- Return molecular weight [kg/kmol]
            virtual scalar W() const;

            //- Return boiling temperature [K]
            virtual scalar Tb(const scalar p) const;


        // Field access

            //- Return density [kg/m^3]
            virtual tmp<volScalarField> rho() const;

            //- Return dynamic viscosity [Pa.s]
            virtual tmp<volScalarField> mu() const;

            //- Return surface tension [kg/s^2]
            virtual tmp<volScalarField> sigma() const;

            //- Return specific heat capacity [J/kg/K]
            virtual tmp<volScalarField> Cp() const;

            //- Return thermal conductivity [W/m/K]
            virtual tmp<volScalarField> kappa() const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const liquidFilmThermo&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
