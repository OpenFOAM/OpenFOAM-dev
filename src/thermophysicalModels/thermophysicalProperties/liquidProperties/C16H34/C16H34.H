/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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
    Foam::C16H34

Description
    nHexaDecane

SourceFiles
    C16H34.C

\*---------------------------------------------------------------------------*/

#ifndef C16H34_H
#define C16H34_H

#include "liquidProperties.H"
#include "NSRDS0ThermophysicalFunction.H"
#include "NSRDS1ThermophysicalFunction.H"
#include "NSRDS2ThermophysicalFunction.H"
#include "NSRDS3ThermophysicalFunction.H"
#include "NSRDS4ThermophysicalFunction.H"
#include "NSRDS5ThermophysicalFunction.H"
#include "NSRDS6ThermophysicalFunction.H"
#include "NSRDS7ThermophysicalFunction.H"
#include "NSRDS14ThermophysicalFunction.H"
#include "APIdiffCoefThermophysicalFunction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class C16H34 Declaration
\*---------------------------------------------------------------------------*/

class C16H34
:
    public liquidProperties
{
    // Private Data

        thermophysicalFunctions::NSRDS5 rho_;
        thermophysicalFunctions::NSRDS1 pv_;
        thermophysicalFunctions::NSRDS6 hl_;
        thermophysicalFunctions::NSRDS0 Cp_;
        thermophysicalFunctions::NSRDS0 h_;
        thermophysicalFunctions::NSRDS7 Cpg_;
        thermophysicalFunctions::NSRDS4 B_;
        thermophysicalFunctions::NSRDS1 mu_;
        thermophysicalFunctions::NSRDS2 mug_;
        thermophysicalFunctions::NSRDS0 kappa_;
        thermophysicalFunctions::NSRDS2 kappag_;
        thermophysicalFunctions::NSRDS6 sigma_;
        thermophysicalFunctions::APIdiffCoef D_;


public:

    friend class liquidProperties;

    //- Runtime type information
    TypeName("C16H34");


    // Constructors

        //- Construct null
        C16H34();

        //- Construct from components
        C16H34
        (
            const liquidProperties& l,
            const thermophysicalFunctions::NSRDS5& density,
            const thermophysicalFunctions::NSRDS1& vapourPressure,
            const thermophysicalFunctions::NSRDS6& heatOfVapourisation,
            const thermophysicalFunctions::NSRDS0& heatCapacity,
            const thermophysicalFunctions::NSRDS0& enthalpy,
            const thermophysicalFunctions::NSRDS7& idealGasHeatCapacity,
            const thermophysicalFunctions::NSRDS4& secondVirialCoeff,
            const thermophysicalFunctions::NSRDS1& dynamicViscosity,
            const thermophysicalFunctions::NSRDS2& vapourDynamicViscosity,
            const thermophysicalFunctions::NSRDS0& thermalConductivity,
            const thermophysicalFunctions::NSRDS2& vapourThermalConductivity,
            const thermophysicalFunctions::NSRDS6& surfaceTension,
            const thermophysicalFunctions::APIdiffCoef& vapourDiffusivity
        );

        //- Construct from dictionary
        C16H34(const dictionary& dict);

        //- Construct and return clone
        virtual autoPtr<liquidProperties> clone() const
        {
            return autoPtr<liquidProperties>(new C16H34(*this));
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

        //- Liquid enthalpy [J/kg]
        inline scalar h(scalar p, scalar T) const;

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

#include "C16H34I.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
