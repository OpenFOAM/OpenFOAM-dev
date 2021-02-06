/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020-2021 OpenFOAM Foundation
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
    Transport properties package using uniform tabulated data for
    viscosity and thermal conductivity vs pressure and temperature.

Usage

    \table
        Property        | Description
        mu              | Dynamic viscosity vs pressure and temperature table
        kappa           | Thermal conductivity vs pressure and temperature table
    \endtable

    Example of the specification of the transport properties:
    \verbatim
    transport
    {
        mu
        {
            pLow    1e3;
            pHigh   1e6;

            Tlow    200;
            Thigh   1500;

            values
            <m> <n>
            (
                (..........)
                .
                .
                .
                (..........)
            );
        }

        kappa
        {
            pLow    1e3;
            pHigh   1e6;

            Tlow    200;
            Thigh   1500;

            values
            <m> <n>
            (
                (..........)
                .
                .
                .
                (..........)
            );
        }
    }
    \endverbatim

SourceFiles
    tabulatedTransportI.H
    tabulatedTransport.C

See also
    Foam::Function2s::uniformTable

\*---------------------------------------------------------------------------*/

#ifndef tabulatedTransport_H
#define tabulatedTransport_H

#include "UniformTable2.H"

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
    typedef Function2s::UniformTable<scalar> table2D;


    // Private Data

        //- Dynamic viscosity table [Pa.s]
        table2D mu_;

        //- Thermal conductivity table [W/m/K]
        table2D kappa_;


    // Private Member Functions

        //- Construct from components
        inline tabulatedTransport
        (
            const Thermo& t,
            const table2D& muPoly,
            const table2D& kappaPoly
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
