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
    Foam::IDEA

Description
    The IDEA fuel is constructed by adding 30% alphaMethylNaphthalene
    with 70% n-decane.

    The new properties have been calculated by adding the values in these
    proportions and making a least square fit, using the same NSRDS-eq.
    so that Y = 0.3*Y_naphthalene + 0.7*Y_decane

    The valid Temperature range for n-decane is normally 243.51 - 617.70 K
    and for the naphthalene it is 242.67 - 772.04 K
    The least square fit was done in the interval 244 - 617 K

    The critical temperature was taken to be 618.074 K, since this
    is the 'c'-value in the rho-equation, which corresponds to Tcrit,
    This value was then used in the fit for the NSRDS6-eq, which uses Tcrit.
    (important for the latent heat and surface tension)

    The molecular weights are 142.20 and 142.285 and for the IDEA fuel
    it is thus 142.26 ( approximately 0.3*142.2 + 0.7*142.285 )

    Critical pressure was set to the lowest one (n-Decane)

    Critical volume... also the lowest one (naphthalene) 0.523 m^3/kmol

    Second Virial Coefficient is n-Decane

SourceFiles
    IDEA.C

\*---------------------------------------------------------------------------*/

#ifndef IDEA_H
#define IDEA_H

#include "liquidProperties.H"
#include "NSRDS0.H"
#include "NSRDS1.H"
#include "NSRDS2.H"
#include "NSRDS3.H"
#include "NSRDS4.H"
#include "NSRDS5.H"
#include "NSRDS6.H"
#include "NSRDS7.H"
#include "APIdiffCoef.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                            Class IDEA Declaration
\*---------------------------------------------------------------------------*/

class IDEA
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
    TypeName("IDEA");


    // Constructors

        //- Construct null
        IDEA();

        // Construct from components
        IDEA
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
        IDEA(const dictionary& dict);

        //- Construct and return clone
        virtual autoPtr<liquidProperties> clone() const
        {
            return autoPtr<liquidProperties>(new IDEA(*this));
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

#include "IDEAI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
