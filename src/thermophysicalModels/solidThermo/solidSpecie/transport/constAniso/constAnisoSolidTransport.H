/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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
    Foam::constAnisoSolidTransport

Description
    Constant properties Transport package.
    Templated into a given Thermodynamics package (needed for thermal
    conductivity).

SourceFiles
    constAnisoSolidTransportI.H
    constAnisoSolidTransport.C

\*---------------------------------------------------------------------------*/

#ifndef constAnisoSolidTransport_H
#define constAnisoSolidTransport_H
#include "vector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Thermo> class constAnisoSolidTransport;

template<class Thermo>
inline constAnisoSolidTransport<Thermo> operator*
(
    const scalar,
    const constAnisoSolidTransport<Thermo>&
);

template<class Thermo>
Ostream& operator<<
(
    Ostream&,
    const constAnisoSolidTransport<Thermo>&
);


/*---------------------------------------------------------------------------*\
                  Class constAnisoSolidTransport Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo>
class constAnisoSolidTransport
:
    public Thermo
{
    // Private Data

        //- Constant anisotropic thermal conductivity.
        vector kappa_;


public:

    // Constructors

        //- Construct from components
        inline constAnisoSolidTransport(const Thermo& t, const vector kappa);

        //- Construct as named copy
        inline constAnisoSolidTransport
        (
            const word&,
            const constAnisoSolidTransport&
        );

        //- Construct from name and dictionary
        constAnisoSolidTransport(const word& name, const dictionary& dict);

        //- Construct and return a clone
        inline autoPtr<constAnisoSolidTransport> clone() const;


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "constAnisoSolid<" + Thermo::typeName() + '>';
        }

        //- Is the thermal conductivity isotropic
        static const bool isotropic = false;

        //- Isotropic thermal conductivity [W/m/K]
        //  Not implemented
        inline scalar kappa(const scalar p, const scalar T) const;

        //- Anisotropic thermal conductivity [W/m/K]
        inline vector Kappa(const scalar p, const scalar T) const;

        //- Write to Ostream
        void write(Ostream& os) const;


    // Member Operators

        inline void operator+=(const constAnisoSolidTransport&);


    // Friend operators

        friend constAnisoSolidTransport operator* <Thermo>
        (
            const scalar,
            const constAnisoSolidTransport&
        );


    // Ostream Operator

        friend Ostream& operator<< <Thermo>
        (
            Ostream&,
            const constAnisoSolidTransport&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "constAnisoSolidTransportI.H"

#ifdef NoRepository
    #include "constAnisoSolidTransport.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
