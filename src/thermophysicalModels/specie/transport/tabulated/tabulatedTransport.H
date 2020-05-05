/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020 OpenFOAM Foundation
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
    Foam::tabulatedTransport

Description
    Transport properties package using non-uniform tabulated data for
    viscosity and thermal conductivity vs temperature.

Usage

    \table
        Property        | Description
        mu              | Dynamic viscosity vs temperature table
        kappa           | Thermal conductivity vs temperature table
    \endtable

    Example of the specification of the transport properties:
    \verbatim
    transport
    {
        mu
        (
            (200 1.82e-05)
            (350 2.61e-05)
            (400 3.93e-05)
        );

        kappa
        (
            (200 2.56e-5)
            (350 3.33e-5)
            (400 4.72e-5)
        );
    }
    \endverbatim

SourceFiles
    tabulatedTransportI.H
    tabulatedTransport.C

See also
    Foam::thermophysicalFunctions::nonUniformTable

\*---------------------------------------------------------------------------*/

#ifndef tabulatedTransport_H
#define tabulatedTransport_H

#include "nonUniformTableThermophysicalFunction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Thermo> class tabulatedTransport;

template<class Thermo>
Ostream& operator<<(Ostream&, const tabulatedTransport<Thermo>&);


/*---------------------------------------------------------------------------*\
                     Class tabulatedTransport Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo>
class tabulatedTransport
:
    public Thermo
{
    typedef thermophysicalFunctions::nonUniformTable nonUniformTable;


    // Private Data

        //- Dynamic viscosity table [Pa.s]
        nonUniformTable mu_;

        //- Thermal conductivity table [W/m/K]
        nonUniformTable kappa_;


    // Private Member Functions

        //- Construct from components
        inline tabulatedTransport
        (
            const Thermo& t,
            const nonUniformTable& muPoly,
            const nonUniformTable& kappaPoly
        );


public:

    // Constructors

        //- Construct as named copy
        inline tabulatedTransport(const word&, const tabulatedTransport&);

        //- Construct from dictionary
        tabulatedTransport(const dictionary& dict);

        //- Construct and return a clone
        inline autoPtr<tabulatedTransport> clone() const;

        // Selector from dictionary
        inline static autoPtr<tabulatedTransport> New(const dictionary& dict);


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "tabulated<" + Thermo::typeName() + '>';
        }

        //- Dynamic viscosity [kg/m/s]
        inline scalar mu(const scalar p, const scalar T) const;

        //- Thermal conductivity [W/m/K]
        inline scalar kappa(const scalar p, const scalar T) const;

        //- Thermal diffusivity of enthalpy [kg/m/s]
        inline scalar alphah(const scalar p, const scalar T) const;

        // Species diffusivity
        // inline scalar D(const scalar p, const scalar T) const;

        //- Write to Ostream
        void write(Ostream& os) const;


    // Ostream Operator

        friend Ostream& operator<< <Thermo>
        (
            Ostream&,
            const tabulatedTransport&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "tabulatedTransportI.H"

#ifdef NoRepository
    #include "tabulatedTransport.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
