/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2023 OpenFOAM Foundation
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
    Foam::WLFTransport

Description
    Transport package using the Williams-Landel-Ferry model for viscosity of
    polymer melts:

    \verbatim
        mu = mu0*exp(-C1*(T - Tr)/(C2 + T - Tr))
    \endverbatim

    References:
    \verbatim
        Williams, M. L., Landel, R. F., & Ferry, J. D. (1955).
        The temperature dependence of relaxation mechanisms
        in amorphous polymers and other glass-forming liquids.
        Journal of the American Chemical society, 77(14), 3701-3707.
    \endverbatim

    The thermal conductivity is obtained using a constant Prandtl number.

Usage
    \table
        Property     | Description
        mu0          | Reference dynamic viscosity [Pa.s]
        Tref         | Reference temperature [K]
        C1           | WLF constant []
        C2           | WLF constant [K]
        Pr           | Constant Prandtl number []
    \endtable

    Example specification of WLFTransport for a polymer:
    \verbatim
    transport
    {
        mu0         50000;
        Tr          416;
        C1          20.4;
        C2          101.6;
        Pr          10000;
    }
    \endverbatim

SourceFiles
    WLFTransportI.H
    WLFTransport.C

\*---------------------------------------------------------------------------*/

#ifndef WLFTransport_H
#define WLFTransport_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Thermo> class WLFTransport;

template<class Thermo>
inline WLFTransport<Thermo> operator+
(
    const WLFTransport<Thermo>&,
    const WLFTransport<Thermo>&
);

template<class Thermo>
inline WLFTransport<Thermo> operator*
(
    const scalar,
    const WLFTransport<Thermo>&
);

template<class Thermo>
Ostream& operator<<
(
    Ostream&,
    const WLFTransport<Thermo>&
);


/*---------------------------------------------------------------------------*\
                        Class WLFTransport Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo>
class WLFTransport
:
    public Thermo
{
    // Private Data

        //- Dynamic viscosity at the reference temperature [Pa.s]
        scalar mu0_;

        //- Reference temperature [T]
        scalar Tr_;

        //- WLF coefficient 1 []
        scalar C1_;

        //- WLF coefficient 2 [T]
        scalar C2_;

        //- Reciprocal Prandtl Number []
        scalar rPr_;


    // Private Member Functions

        //- Read coefficient from dictionary
        scalar readCoeff(const word& coeffName, const dictionary& dict);


public:

    // Constructors

        //- Construct from components
        inline WLFTransport
        (
            const Thermo& t,
            const scalar mu0,
            const scalar Tr,
            const scalar C1,
            const scalar C2,
            const scalar Pr
        );

        //- Construct as named copy
        inline WLFTransport(const word&, const WLFTransport&);

        //- Construct from name and dictionary
        WLFTransport(const word& name, const dictionary& dict);

        //- Construct and return a clone
        inline autoPtr<WLFTransport> clone() const;


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "WLF<" + Thermo::typeName() + '>';
        }

        //- Dynamic viscosity [kg/m/s]
        inline scalar mu(const scalar p, const scalar T) const;

        //- Thermal conductivity [W/m/K]
        inline scalar kappa(const scalar p, const scalar T) const;

        //- Write to Ostream
        void write(Ostream& os) const;


    // Member Operators

        inline void operator+=(const WLFTransport&);
        inline void operator*=(const scalar);


    // Friend operators

        friend WLFTransport operator+ <Thermo>
        (
            const WLFTransport&,
            const WLFTransport&
        );

        friend WLFTransport operator* <Thermo>
        (
            const scalar,
            const WLFTransport&
        );


    // Ostream Operator

        friend Ostream& operator<< <Thermo>
        (
            Ostream&,
            const WLFTransport&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "WLFTransportI.H"

#ifdef NoRepository
    #include "WLFTransport.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
