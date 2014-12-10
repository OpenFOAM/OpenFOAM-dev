/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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
    Foam::thermo

Description
    Basic thermodynamics type based on the use of fitting functions for
    cp, h, s obtained from the template argument type thermo.  All other
    properties are derived from these primitive functions.

SourceFiles
    thermoI.H
    thermo.C

\*---------------------------------------------------------------------------*/

#ifndef thermo_H
#define thermo_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace species
{

// Forward declaration of friend functions and operators

template<class Thermo, template<class> class Type> class thermo;

template<class Thermo, template<class> class Type>
inline thermo<Thermo, Type> operator+
(
    const thermo<Thermo, Type>&,
    const thermo<Thermo, Type>&
);

template<class Thermo, template<class> class Type>
inline thermo<Thermo, Type> operator-
(
    const thermo<Thermo, Type>&,
    const thermo<Thermo, Type>&
);

template<class Thermo, template<class> class Type>
inline thermo<Thermo, Type> operator*
(
    const scalar,
    const thermo<Thermo, Type>&
);

template<class Thermo, template<class> class Type>
inline thermo<Thermo, Type> operator==
(
    const thermo<Thermo, Type>&,
    const thermo<Thermo, Type>&
);

template<class Thermo, template<class> class Type>
Ostream& operator<<
(
    Ostream&,
    const thermo<Thermo, Type>&
);


/*---------------------------------------------------------------------------*\
                           Class thermo Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo, template<class> class Type>
class thermo
:
    public Thermo,
    public Type<thermo<Thermo, Type> >
{
    // Private data

        //- Convergence tolerance of energy -> temperature inversion functions
        static const scalar tol_;

        //- Max number of iterations in energy->temperature inversion functions
        static const int maxIter_;


    // Private Member Functions

        //- Return the temperature corresponding to the value of the
        //  thermodynamic property f, given the function f = F(p, T)
        //  and dF(p, T)/dT
        inline scalar T
        (
            scalar f,
            scalar p,
            scalar T0,
            scalar (thermo::*F)(const scalar, const scalar) const,
            scalar (thermo::*dFdT)(const scalar, const scalar) const,
            scalar (thermo::*limit)(const scalar) const
        ) const;


public:

    //- The thermodynamics of the individual species'
    typedef thermo<Thermo, Type> thermoType;


    // Constructors

        //- construct from components
        inline thermo(const Thermo& sp);

        //- Construct from Istream
        thermo(Istream&);

        //- Construct from dictionary
        thermo(const dictionary& dict);

        //- Construct as named copy
        inline thermo(const word& name, const thermo&);


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return
                  Thermo::typeName() + ','
                + Type<thermo<Thermo, Type> >::typeName();
        }

        // Fundamental properties
        // (These functions must be provided in derived types)

            // Heat capacity at constant pressure [J/(kmol K)]
            // scalar cp(const scalar) const;

            // Absolute Enthalpy [J/kmol]
            // scalar ha(const scalar) const;

            // Sensible enthalpy [J/kmol]
            // scalar hs(const scalar) const;

            // Chemical enthalpy [J/kmol]
            // scalar hc(const scalar) const;

            // Entropy [J/(kmol K)]
            // scalar s(const scalar) const;


        // Calculate and return derived properties
        // (These functions need not provided in derived types)

            // Mole specific properties

                //- Name of Enthalpy/Internal energy
                static inline word heName();

                //- Enthalpy/Internal energy [J/kmol]
                inline scalar he(const scalar p, const scalar T) const;

                //- Heat capacity at constant volume [J/(kmol K)]
                inline scalar cv(const scalar p, const scalar T) const;

                //- Heat capacity at constant pressure/volume [J/(kmol K)]
                inline scalar cpv(const scalar p, const scalar T) const;

                //- gamma = cp/cv []
                inline scalar gamma(const scalar p, const scalar T) const;

                //- Ratio of heat capacity at constant pressure to that at
                //  constant pressure/volume []
                inline scalar cpBycpv(const scalar p, const scalar T) const;

                //- Sensible internal energy [J/kmol]
                inline scalar es(const scalar p, const scalar T) const;

                //- Absolute internal energy [J/kmol]
                inline scalar ea(const scalar p, const scalar T) const;

                //- Gibbs free energy [J/kmol]
                inline scalar g(const scalar p, const scalar T) const;

                //- Helmholtz free energy [J/kmol]
                inline scalar a(const scalar p, const scalar T) const;


            // Mass specific properties

                //- Heat capacity at constant pressure [J/(kg K)]
                inline scalar Cp(const scalar p, const scalar T) const;

                //- Heat capacity at constant volume [J/(kg K)]
                inline scalar Cv(const scalar p, const scalar T) const;

                //- Heat capacity at constant pressure/volume [J/(kg K)]
                inline scalar Cpv(const scalar p, const scalar T) const;

                //- Enthalpy/Internal energy [J/kg]
                inline scalar HE(const scalar p, const scalar T) const;

                //- Sensible enthalpy [J/kg]
                inline scalar Hs(const scalar p, const scalar T) const;

                //- Chemical enthalpy [J/kg]
                inline scalar Hc() const;

                //- Absolute Enthalpy [J/kg]
                inline scalar Ha(const scalar p, const scalar T) const;

                //- Entropy [J/(kg K)]
                inline scalar S(const scalar p, const scalar T) const;

                //- Internal energy [J/kg]
                inline scalar E(const scalar p, const scalar T) const;

                //- Sensible internal energy [J/kg]
                inline scalar Es(const scalar p, const scalar T) const;

                //- Absolute internal energy [J/kg]
                inline scalar Ea(const scalar p, const scalar T) const;

                //- Gibbs free energy [J/kg]
                inline scalar G(const scalar p, const scalar T) const;

                //- Helmholtz free energy [J/kg]
                inline scalar A(const scalar p, const scalar T) const;


        // Equilibrium reaction thermodynamics

            //- Equilibrium constant [] i.t.o fugacities
            //  = PIi(fi/Pstd)^nui
            inline scalar K(const scalar p, const scalar T) const;

            //- Equilibrium constant [] i.t.o. partial pressures
            //  = PIi(pi/Pstd)^nui
            //  For low pressures (where the gas mixture is near perfect) Kp = K
            inline scalar Kp(const scalar p, const scalar T) const;

            //- Equilibrium constant i.t.o. molar concentration
            //  = PIi(ci/cstd)^nui
            //  For low pressures (where the gas mixture is near perfect)
            //  Kc = Kp(pstd/(RR*T))^nu
            inline scalar Kc(const scalar p, const scalar T) const;

            //- Equilibrium constant [] i.t.o. mole-fractions
            //  For low pressures (where the gas mixture is near perfect)
            //  Kx = Kp(pstd/p)^nui
            inline scalar Kx
            (
                const scalar p,
                const scalar T
            ) const;

            //- Equilibrium constant [] i.t.o. number of moles
            //  For low pressures (where the gas mixture is near perfect)
            //  Kn = Kp(n*pstd/p)^nui where n = number of moles in mixture
            inline scalar Kn
            (
                const scalar p,
                const scalar T,
                const scalar n
            ) const;


        // Energy->temperature  inversion functions

            //- Temperature from enthalpy or internal energy
            //  given an initial temperature T0
            inline scalar THE
            (
                const scalar H,
                const scalar p,
                const scalar T0
            ) const;

            //- Temperature from sensible enthalpy given an initial T0
            inline scalar THs
            (
                const scalar Hs,
                const scalar p,
                const scalar T0
            ) const;

            //- Temperature from absolute enthalpy
            //  given an initial temperature T0
            inline scalar THa
            (
                const scalar H,
                const scalar p,
                const scalar T0
            ) const;

            //- Temperature from sensible internal energy
            //  given an initial temperature T0
            inline scalar TEs
            (
                const scalar E,
                const scalar p,
                const scalar T0
            ) const;

            //- Temperature from absolute internal energy
            //  given an initial temperature T0
            inline scalar TEa
            (
                const scalar E,
                const scalar p,
                const scalar T0
            ) const;


        // I-O

            //- Write to Ostream
            void write(Ostream& os) const;


    // Member operators

        inline void operator+=(const thermo&);
        inline void operator-=(const thermo&);

        inline void operator*=(const scalar);


    // Friend operators

        friend thermo operator+ <Thermo, Type>
        (
            const thermo&,
            const thermo&
        );

        friend thermo operator- <Thermo, Type>
        (
            const thermo&,
            const thermo&
        );

        friend thermo operator* <Thermo, Type>
        (
            const scalar s,
            const thermo&
        );

        friend thermo operator== <Thermo, Type>
        (
            const thermo&,
            const thermo&
        );


    // Ostream Operator

        friend Ostream& operator<< <Thermo, Type>
        (
            Ostream&,
            const thermo&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace species
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "thermoI.H"

#ifdef NoRepository
#   include "thermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
