/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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
    Foam::constantFilmThermo

Description
    Constant thermo model

SourceFiles
    constantFilmThermo.C

\*---------------------------------------------------------------------------*/

#ifndef constantFilmThermo_H
#define constantFilmThermo_H

#include "filmThermoModel.H"
#include "dimensionSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

/*---------------------------------------------------------------------------*\
                     Class constantFilmThermo Declaration
\*---------------------------------------------------------------------------*/

class constantFilmThermo
:
    public filmThermoModel
{
public:

    struct thermoData
    {
//    private:
        word name_;
        scalar value_;
        bool set_;

//    public:
        thermoData()
        :
            name_("unknown"),
            value_(0.0),
            set_(false)
        {}
        thermoData(const word& n)
        :
            name_(n),
            value_(0.0),
            set_(false)
        {}

//        virtual ~thermoData()
//        {}
    };


private:

    // Private data

        //- Specie name
        word name_;

        //- Density [kg/m3]
        mutable thermoData rho0_;

        //- Dynamic viscosity [Pa.s]
        mutable thermoData mu0_;

        //- Surface tension [kg/s2]
        mutable thermoData sigma0_;

        //- Specific heat capacity [J/kg/K]
        mutable thermoData Cp0_;

        //- Thermal conductivity [W/m/K]
        mutable thermoData kappa0_;

        //- Diffusivity [m2/s]
        mutable thermoData D0_;

        //- Latent heat [J/kg]
        mutable thermoData hl0_;

        //- Vapour pressure [Pa]
        mutable thermoData pv0_;

        //- Molecular weight [kg/kmol]
        mutable thermoData W0_;

        //- Boiling temperature [K]
        mutable thermoData Tb0_;


    // Private member functions

        //- Initialise thermoData object
        void init(thermoData& td);

        //- Disallow default bitwise copy construct
        constantFilmThermo(const constantFilmThermo&);

        //- Disallow default bitwise assignment
        void operator=(const constantFilmThermo&);


public:

    //- Runtime type information
    TypeName("constant");


    // Constructors

        //- Construct from surface film model and dictionary
        constantFilmThermo
        (
            surfaceFilmModel& owner,
            const dictionary& dict
        );


    //- Destructor
    virtual ~constantFilmThermo();


    // Member Functions

        //- Return the specie name
        virtual const word& name() const;


        // Elemental access

            //- Return density [kg/m3]
            virtual scalar rho(const scalar p, const scalar T) const;

            //- Return dynamic viscosity [Pa.s]
            virtual scalar mu(const scalar p, const scalar T) const;

            //- Return surface tension [kg/s2]
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

            //- Return density [kg/m3]
            virtual tmp<volScalarField> rho() const;

            //- Return dynamic viscosity [Pa.s]
            virtual tmp<volScalarField> mu() const;

            //- Return surface tension [kg/s2]
            virtual tmp<volScalarField> sigma() const;

            //- Return specific heat capacity [J/kg/K]
            virtual tmp<volScalarField> Cp() const;

            //- Return thermal conductivity [W/m/K]
            virtual tmp<volScalarField> kappa() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
