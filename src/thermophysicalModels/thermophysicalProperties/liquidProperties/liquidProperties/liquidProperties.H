/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    Foam::liquidProperties

Description
    The thermophysical properties of a liquid

SourceFiles
    liquidProperties.C

\*---------------------------------------------------------------------------*/

#ifndef liquidProperties_H
#define liquidProperties_H

#include "thermophysicalProperties.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class liquidProperties Declaration
\*---------------------------------------------------------------------------*/

class liquidProperties
:
    public thermophysicalProperties
{
    // Private data

        //- Critical temperature [K]
        scalar Tc_;

        //- Critical pressure [Pa]
        scalar Pc_;

        //- Critical volume [m^3/kmol]
        scalar Vc_;

        //- Critical compressibility factor []
        scalar Zc_;

        //- Triple point temperature [K]
        scalar Tt_;

        //- Triple point pressure [Pa]
        scalar Pt_;

        //- Normal boiling temperature [K]
        scalar Tb_;

        //- Dipole moment []
        scalar dipm_;

        //- Pitzer's accentric factor []
        scalar omega_;

        //- Solubility parameter [(J/m^3)^0.5]
        scalar delta_;


public:

    TypeName("liquid");


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            autoPtr,
            liquidProperties,
            ,
            (),
            ()
        );

        declareRunTimeSelectionTable
        (
            autoPtr,
            liquidProperties,
            dictionary,
            (const dictionary& dict),
            (dict)
        );


    // Constructors

        //- Construct from components
        liquidProperties
        (
            scalar W,
            scalar Tc,
            scalar Pc,
            scalar Vc,
            scalar Zc,
            scalar Tt,
            scalar Pt,
            scalar Tb,
            scalar dipm,
            scalar omega,
            scalar delta
        );

        //- Construct from dictionary
        liquidProperties(const dictionary& dict);

        //- Construct and return clone
        virtual autoPtr<liquidProperties> clone() const = 0;


    // Selectors

        //- Return a pointer to a new liquidProperties created from name
        static autoPtr<liquidProperties> New(const word& name);

        //- Return a pointer to a new liquidProperties created from dictionary
        static autoPtr<liquidProperties> New(const dictionary& dict);


    //- Destructor
    virtual ~liquidProperties()
    {}


    // Static data

        //- Is the equation of state is incompressible i.e. rho != f(p)
        static const bool incompressible = true;

        //- Is the equation of state is isochoric i.e. rho = const
        static const bool isochoric = false;


    // Member Functions

        // Physical constants which define the specie

            //- No of moles of this species in mixture
            //  Note Mixing of liquidProperties is not currently supported
            //  so Y = 1
            inline scalar Y() const;

            //- Critical temperature [K]
            inline scalar Tc() const;

            //- Critical pressure [Pa]
            inline scalar Pc() const;

            //- Critical volume [m^3/kmol]
            inline scalar Vc() const;

            //- Critical compressibilty factor
            inline scalar Zc() const;

            //- Triple point temperature [K]
            inline scalar Tt() const;

            //- Triple point pressure [Pa]
            inline scalar Pt() const;

            //- Normal boiling temperature [K]
            inline scalar Tb() const;

            //- Dipole moment []
            inline scalar dipm() const;

            //- Pitzer's ascentric factor []
            inline scalar omega() const;

            //- Solubility parameter [(J/m^3)^(1/2)]
            inline scalar delta() const;

            //- Limit the temperature to be in the range Tlow_ to Thigh_
            inline scalar limit(const scalar T) const;


        // Fundamental equation of state properties

            //- Liquid compressibility rho/p [s^2/m^2]
            //  Note: currently it is assumed the liquid is incompressible
            inline scalar psi(scalar p, scalar T) const;

            //- Return (Cp - Cv) [J/(kg K]
            //  Note: currently it is assumed the liquid is incompressible
            //  so CpMCv = 0
            inline scalar CpMCv(scalar p, scalar T) const;


        // Fundamental thermodynamic properties

            //- Absolute Enthalpy [J/kg]
            inline scalar Ha(const scalar p, const scalar T) const;

            //- Sensible enthalpy [J/kg]
            inline scalar Hs(const scalar p, const scalar T) const;

            //- Chemical enthalpy [J/kg]
            inline scalar Hc() const;

            // Entropy [J/(kg K)]
            scalar S(const scalar p, const scalar T) const;


        // Physical properties

            //- Vapour pressure [Pa]
            virtual scalar pv(scalar p, scalar T) const = 0;

            //- Heat of vapourisation [J/kg]
            virtual scalar hl(scalar p, scalar T) const = 0;

            //- Liquid enthalpy [J/kg] - reference to 298.15 K
            virtual scalar h(scalar p, scalar T) const = 0;

            //- Vapour heat capacity [J/(kg K)]
            virtual scalar Cpg(scalar p, scalar T) const = 0;

            //- Liquid viscosity [Pa s]
            virtual scalar mu(scalar p, scalar T) const = 0;

            //- Vapour viscosity [Pa s]
            virtual scalar mug(scalar p, scalar T) const = 0;

            //- Liquid thermal conductivity  [W/(m K)]
            virtual scalar kappa(scalar p, scalar T) const = 0;

            //- Liquid thermal diffusivity of enthalpy [kg/ms]
            inline scalar alphah(const scalar p, const scalar T) const;

            //- Vapour thermal conductivity  [W/(m K)]
            virtual scalar kappag(scalar p, scalar T) const = 0;

            //- Surface tension [N/m]
            virtual scalar sigma(scalar p, scalar T) const = 0;

            //- Vapour diffussivity [m2/s]
            virtual scalar D(scalar p, scalar T) const = 0;

            //- Vapour diffussivity [m2/s] with specified binary pair
            virtual scalar D(scalar p, scalar T, scalar Wb) const = 0;

            //- Invert the vapour pressure relationship to retrieve the
            //  boiling temperuture as a function of pressure
            virtual scalar pvInvert(scalar p) const;


    // I-O

        //- Read and set the properties present it the given dictionary
        void readIfPresent(const dictionary& dict);

        //- Read and set the function coefficients
        //  if present it the given dictionary
        template<class Func>
        inline void readIfPresent
        (
            Func& f,
            const word& name,
            const dictionary& dict
        );

        //- Read and set the function coefficients
        //  if present it the given dictionary
        template<class Liquid>
        inline void readIfPresent(Liquid& l, const dictionary& dict);

        //- Write the function coefficients
        virtual void writeData(Ostream& os) const = 0;

        //- Write dictionary to Ostream
        void write(Ostream& os) const;

        //- Write the data for each of the property functions
        template<class Liquid>
        inline void writeData(const Liquid& l, Ostream& os) const;

        //- Ostream Operator
        friend Ostream& operator<<(Ostream& os, const liquidProperties& l);
};


Ostream& operator<<(Ostream& os, const liquidProperties& l);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "liquidPropertiesI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
