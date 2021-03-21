/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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
    Foam::C7H8

Description
    toluene

SourceFiles
    C7H8.C

\*---------------------------------------------------------------------------*/

#ifndef C7H8_H
#define C7H8_H

#include "liquidProperties.H"
#include "NSRDS0.H"
#include "NSRDS1.H"
#include "NSRDS2.H"
#include "NSRDS3.H"
#include "NSRDS4.H"
#include "NSRDS5.H"
#include "NSRDS6.H"
#include "NSRDS7.H"
#include "NSRDS14.H"
#include "APIdiffCoef.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                            Class C7H8 Declaration
\*---------------------------------------------------------------------------*/

class C7H8
:
    public liquidProperties
{
    // Private Data

        Function1s::NSRDS5 rho_;
        Function1s::NSRDS1 pv_;
        Function1s::NSRDS6 hl_;
        Function1s::NSRDS0 Cp_;
        Function1s::NSRDS0 h_;
        Function1s::NSRDS7 Cpg_;
        Function1s::NSRDS4 B_;
        Function1s::NSRDS1 mu_;
        Function1s::NSRDS2 mug_;
        Function1s::NSRDS0 kappa_;
        Function1s::NSRDS2 kappag_;
        Function1s::NSRDS6 sigma_;
        Function2s::APIdiffCoef D_;

        //- Liquid heat of formation [J/kg]
        scalar Hf_;


public:

    friend class liquidProperties;

    //- Runtime type information
    TypeName("C7H8");


    // Constructors

        //- Construct null
        C7H8();

        //- Construct from components
        C7H8
        (
            const liquidProperties& l,
            const Function1s::NSRDS5& density,
            const Function1s::NSRDS1& vapourPressure,
            const Function1s::NSRDS6& heatOfVapourisation,
            const Function1s::NSRDS0& heatCapacity,
            const Function1s::NSRDS0& enthalpy,
            const Function1s::NSRDS7& idealGasHeatCapacity,
            const Function1s::NSRDS4& secondVirialCoeff,
            const Function1s::NSRDS1& dynamicViscosity,
            const Function1s::NSRDS2& vapourDynamicViscosity,
            const Function1s::NSRDS0& thermalConductivity,
            const Function1s::NSRDS2& vapourThermalConductivity,
            const Function1s::NSRDS6& surfaceTension,
            const Function2s::APIdiffCoef& vapourDiffusivity
        );

        //- Construct from dictionary
        C7H8(const dictionary& dict);

        //- Construct and return clone
        virtual autoPtr<liquidProperties> clone() const
        {
            return autoPtr<liquidProperties>(new C7H8(*this));
        }


    // Member Functions

        //- Liquid density [kg/m^3]
        inline scalar rho(scalar p, scalar T) const;

        //- Vapour pressure [Pa]
        inline scalar pv(scalar p, scalar T) const;

        //- Heat of vapourisation [J/kg]
        inline scalar hl(scalar p, scalar T) const;

        //- Liquid heat capacity [J/kg/K]
        inline scalar Cp(scalar p, scalar T) const;

        //- Liquid sensible enthalpy [J/kg]
        inline scalar Hs(scalar p, scalar T) const;

        //- Liquid heat of formation [J/kg]
        inline scalar Hf() const;

        //- Liquid absolute enthalpy [J/kg]
        inline scalar Ha(scalar p, scalar T) const;

        //- Ideal gas heat capacity [J/kg/K]
        inline scalar Cpg(scalar p, scalar T) const;

        //- Second Virial Coefficient [m^3/kg]
        inline scalar B(scalar p, scalar T) const;

        //- Liquid viscosity [Pa s]
        inline scalar mu(scalar p, scalar T) const;

        //- Vapour viscosity [Pa s]
        inline scalar mug(scalar p, scalar T) const;

        //- Liquid thermal conductivity  [W/m/K]
        inline scalar kappa(scalar p, scalar T) const;

        //- Vapour thermal conductivity  [W/m/K]
        inline scalar kappag(scalar p, scalar T) const;

        //- Surface tension [N/m]
        inline scalar sigma(scalar p, scalar T) const;

        //- Vapour diffusivity [m^2/s]
        inline scalar D(scalar p, scalar T) const;

        //- Vapour diffusivity [m^2/s] with specified binary pair
        inline scalar D(scalar p, scalar T, scalar Wb) const;


    // I-O

        //- Write the function coefficients
        void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "C7H8I.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
