/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2023 OpenFOAM Foundation
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
    Foam::tabulatedSolidTransport

Description
    Transport properties package using non-uniformly-spaced tabulated data for
    thermal conductivity vs temperature.

Usage
    \table
        Property        | Description
        kappa           | Thermal conductivity vs temperature table
    \endtable

    Example of the specification of the transport properties:
    \verbatim
    transport
    {
        kappa
        {
            values
            (
                (200 380)
                (350 400)
                (400 450)
            );
        }
    }
    \endverbatim

SourceFiles
    tabulatedSolidTransportI.H
    tabulatedSolidTransport.C

\*---------------------------------------------------------------------------*/

#ifndef tabulatedSolidTransport_H
#define tabulatedSolidTransport_H

#include "NonUniformTable1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Thermo> class tabulatedSolidTransport;

template<class Thermo>
inline tabulatedSolidTransport<Thermo> operator*
(
    const scalar,
    const tabulatedSolidTransport<Thermo>&
);

template<class Thermo>
Ostream& operator<<
(
    Ostream&,
    const tabulatedSolidTransport<Thermo>&
);


/*---------------------------------------------------------------------------*\
                  Class tabulatedSolidTransport Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo>
class tabulatedSolidTransport
:
    public Thermo
{
    typedef Function1s::NonUniformTable<scalar> nonUniformTable;


    // Private Data

        //- Thermal conductivity table [W/m/K]
        nonUniformTable kappa_;


public:

    // Constructors

        //- Construct from components
        inline tabulatedSolidTransport
        (
            const Thermo& t,
            const nonUniformTable& kappa
        );

        //- Construct as named copy
        inline tabulatedSolidTransport
        (
            const word&,
            const tabulatedSolidTransport&
        );

        //- Construct from name and dictionary
        tabulatedSolidTransport(const word& name, const dictionary& dict);

        //- Construct and return a clone
        inline autoPtr<tabulatedSolidTransport> clone() const;


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "tabulatedSolid<" + Thermo::typeName() + '>';
        }

        //- Is the thermal conductivity isotropic
        static const bool isotropic = true;

        //- Isotropic thermal conductivity [W/m/K]
        inline scalar kappa(const scalar p, const scalar T) const;

        //- Anisotropic thermal conductivity [W/m/K]
        //  Not implemented
        inline vector Kappa(const scalar p, const scalar T) const;

        //- Write to Ostream
        void write(Ostream& os) const;


    // Ostream Operator

        friend Ostream& operator<< <Thermo>
        (
            Ostream&,
            const tabulatedSolidTransport&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "tabulatedSolidTransportI.H"

#ifdef NoRepository
    #include "tabulatedSolidTransport.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
