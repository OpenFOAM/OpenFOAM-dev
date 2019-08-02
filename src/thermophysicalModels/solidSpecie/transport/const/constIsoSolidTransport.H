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
    Foam::constIsoSolidTransport

Description
    Constant properties Transport package.
    Templated into a given thermodynamics package (needed for thermal
    conductivity).

SourceFiles
    constIsoSolidTransportI.H
    constIsoSolidTransport.C

\*---------------------------------------------------------------------------*/

#ifndef constIsoSolidTransport_H
#define constIsoSolidTransport_H

#include "vector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Thermo> class constIsoSolidTransport;

template<class Thermo>
inline constIsoSolidTransport<Thermo> operator*
(
    const scalar,
    const constIsoSolidTransport<Thermo>&
);

template<class Thermo>
Ostream& operator<<
(
    Ostream&,
    const constIsoSolidTransport<Thermo>&
);


/*---------------------------------------------------------------------------*\
                   Class constIsoSolidTransport Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo>
class constIsoSolidTransport
:
    public Thermo
{
    // Private Data

        //- Constant isotropic thermal conductivity
        scalar kappa_;


    // Private Member Functions

        //- Construct from components
        inline constIsoSolidTransport(const Thermo& t, const scalar kappa);


public:

    // Constructors

        //- Construct as named copy
        inline constIsoSolidTransport
        (
            const word&,
            const constIsoSolidTransport&
        );

        //- Construct from dictionary
        constIsoSolidTransport(const dictionary& dict);

        //- Construct and return a clone
        inline autoPtr<constIsoSolidTransport> clone() const;

        // Selector from dictionary
        inline static autoPtr<constIsoSolidTransport> New
        (
            const dictionary& dict
        );


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "constIso<" + Thermo::typeName() + '>';
        }

        //- Is the thermal conductivity isotropic
        static const bool isotropic = true;

        //- Isotropic thermal conductivity [W/m/K]
        inline scalar kappa(const scalar p, const scalar T) const;

        //- Un-isotropic thermal conductivity [W/m/K]
        inline vector Kappa(const scalar p, const scalar T) const;

         //- Dynamic viscosity [kg/m/s]
        inline scalar mu(const scalar p, const scalar T) const;

        //- Thermal diffusivity of enthalpy [kg/m/s]
        inline scalar alphah(const scalar p, const scalar T) const;


        //- Write to Ostream
        void write(Ostream& os) const;


    // Member Operators

        inline void operator+=(const constIsoSolidTransport&);


    // Friend operators

        friend constIsoSolidTransport operator* <Thermo>
        (
            const scalar,
            const constIsoSolidTransport&
        );


    // Ostream Operator

        friend Ostream& operator<< <Thermo>
        (
            Ostream&,
            const constIsoSolidTransport&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "constIsoSolidTransportI.H"

#ifdef NoRepository
    #include "constIsoSolidTransport.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
