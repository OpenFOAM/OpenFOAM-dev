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
    Foam::NSRDSThermo

Description
    Enthalpy based thermodynamics package using NSRDS function number 100 to
    calculate specific heat capacity at constant pressure as a function of the
    temperature.

    Reference:
    \verbatim
    Daubert, T. E., & Danner, R. P. (1985).
    Data compilation tables of properties of pure compounds.
    \endverbatim

Usage
    \table
        Property     | Description
        hf           | Heat of formation
        sf           | Standard entropy
        a            | Coefficient
        b            | Coefficient
        c            | Coefficient
        d            | Coefficient
        e            | Coefficient
    \endtable

    Example specification of NSRDS thermodynamics for water, with coefficients
    as given by Daubert & Danner:
    \verbatim
    thermodynamics
    {
        hf              -15858419.8191934;
        sf              0; // ?
        a               2.92167638e+03;
        b               1.33882875e+01;
        c               -4.72300860e-02;
        d               5.55092978e-05;
        e               0;
    }
    \endverbatim

SourceFiles
    NSRDSThermoI.H
    NSRDSThermo.C

See also
    Foam::FixedPolynomial

\*---------------------------------------------------------------------------*/

#ifndef NSRDSThermo_H
#define NSRDSThermo_H

#include "scalar.H"
#include "FixedLaurentPolynomial.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class EquationOfState>
class NSRDSThermo;

template<class EquationOfState>
inline NSRDSThermo<EquationOfState> operator+
(
    const NSRDSThermo<EquationOfState>&,
    const NSRDSThermo<EquationOfState>&
);

template<class EquationOfState>
inline NSRDSThermo<EquationOfState> operator*
(
    const scalar,
    const NSRDSThermo<EquationOfState>&
);

template<class EquationOfState>
inline NSRDSThermo<EquationOfState> operator==
(
    const NSRDSThermo<EquationOfState>&,
    const NSRDSThermo<EquationOfState>&
);

template<class EquationOfState>
Ostream& operator<<(Ostream&, const NSRDSThermo<EquationOfState>&);


/*---------------------------------------------------------------------------*\
                      Class NSRDSThermo Declaration
\*---------------------------------------------------------------------------*/

template<class EquationOfState>
class NSRDSThermo
:
    public EquationOfState
{
    // Private Data

        //- Heat of formation
        scalar hf_;

        //- Standard entropy
        scalar sf_;

        //- Specific heat at constant pressure polynomial
        FixedLaurentPolynomial<scalar, 0, 5> CpCoeffs_;

        //- Reference sensible enthalpy
        scalar hsRef_;

        //- Reference entropy
        scalar sRef_;


    // Private Constructors

        //- Construct from name and parent dictionary and dictionary
        NSRDSThermo
        (
            const word& name,
            const dictionary& dict,
            const dictionary& subDict
        );


public:

    // Constructors

        //- Construct from components
        inline NSRDSThermo
        (
            const EquationOfState& ns,
            const scalar hf,
            const scalar sf,
            const FixedLaurentPolynomial<scalar, 0, 5>& CpCoeffs,
            const scalar hsRef,
            const scalar sRef
        );

        //- Construct from name and dictionary
        NSRDSThermo(const word& name, const dictionary& dict);

        //- Construct as a named copy
        inline NSRDSThermo(const word&, const NSRDSThermo&);


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "NSRDS<" + EquationOfState::typeName() + '>';
        }

        //- Limit the temperature to be in the range Tlow_ to Thigh_
        inline scalar limit(const scalar) const;


        // Fundamental properties

            //- Heat capacity at constant pressure [J/kg/K]
            inline scalar Cp(const scalar p, const scalar T) const;

            //- Absolute enthalpy [J/kg]
            inline scalar ha(const scalar p, const scalar T) const;

            //- Sensible enthalpy [J/kg]
            inline scalar hs(const scalar p, const scalar T) const;

            //- Enthalpy of formation [J/kg]
            inline scalar hf() const;

            //- Entropy [J/kg/K]
            inline scalar s(const scalar p, const scalar T) const;

            //- Gibbs free energy of the mixture in the standard state [J/kg]
            inline scalar gStd(const scalar T) const;

            #include "HtoEthermo.H"


        // Derivative term used for Jacobian

            //- Temperature derivative of heat capacity at constant pressure
            inline scalar dCpdT(const scalar p, const scalar T) const;


        // I-O

            //- Write to Ostream
            void write(Ostream& os) const;


    // Member Operators

        inline void operator+=(const NSRDSThermo&);
        inline void operator*=(const scalar);


    // Friend operators

        friend NSRDSThermo operator+ <EquationOfState>
        (
            const NSRDSThermo&,
            const NSRDSThermo&
        );

        friend NSRDSThermo operator* <EquationOfState>
        (
            const scalar,
            const NSRDSThermo&
        );

        friend NSRDSThermo operator== <EquationOfState>
        (
            const NSRDSThermo&,
            const NSRDSThermo&
        );


    // Ostream Operator

        friend Ostream& operator<< <EquationOfState>
        (
            Ostream&,
            const NSRDSThermo&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "NSRDSThermoI.H"

#ifdef NoRepository
    #include "NSRDSThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
