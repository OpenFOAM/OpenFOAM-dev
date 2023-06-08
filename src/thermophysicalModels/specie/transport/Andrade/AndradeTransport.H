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
    Foam::AndradeTransport

Description
    Transport package using the Andrade function for the natural logarithm of
    dynamic viscosity and thermal conductivity of liquids:

    \verbatim
        log(mu) = muCoeffs[0] + muCoeffs[1]*T + muCoeffs[2]*sqr(T)
          + muCoeffs_[3]/(muCoeffs_[4] + T)

        log(kappa) = kappaCoeffs[0] + kappaCoeffs[1]*T + kappaCoeffs[2]*sqr(T)
          + kappaCoeffs_[3]/(kappaCoeffs_[4] + T)
    );
    \endverbatim

    References:
    \verbatim
        Andrade, E. D. C. (1934).
        XLI. A theory of the viscosity of liquids.—Part I.
        The London, Edinburgh, and Dublin Philosophical Magazine
        and Journal of Science, 17(112), 497-511.

        Andrade, E. D. C. (1934).
        LVIII. A theory of the viscosity of liquids.—Part II.
        The London, Edinburgh, and Dublin Philosophical Magazine
        and Journal of Science, 17(113), 698-732.
    \endverbatim

Usage
    \table
        Property        | Description
        muCoeffs        | Dynamic viscosity polynomial coefficients
        kappaCoeffs     | Thermal conductivity polynomial coefficients
    \endtable

    Example of the specification of the transport properties for water@200bar:
    \verbatim
    transport
    {
        muCoeffs    (-25.8542 0.031256 -2.2e-05 3289.918 -11.4784);
        kappaCoeffs (-2.56543 0.008794 -9.8e-06 100.368 0);
    }
    \endverbatim

SourceFiles
    AndradeTransportI.H
    AndradeTransport.C

\*---------------------------------------------------------------------------*/

#ifndef AndradeTransport_H
#define AndradeTransport_H

#include "FixedList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Thermo> class AndradeTransport;

template<class Thermo>
inline AndradeTransport<Thermo> operator+
(
    const AndradeTransport<Thermo>&,
    const AndradeTransport<Thermo>&
);

template<class Thermo>
inline AndradeTransport<Thermo> operator*
(
    const scalar,
    const AndradeTransport<Thermo>&
);

template<class Thermo>
Ostream& operator<<
(
    Ostream&,
    const AndradeTransport<Thermo>&
);


/*---------------------------------------------------------------------------*\
                   Class AndradeTransport Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo>
class AndradeTransport
:
    public Thermo
{
    // Private Typedefs

        //- Coefficient list type
        typedef FixedList<scalar, 5> coeffList;


    // Private Data

        //- Dynamic viscosity coefficients
        coeffList muCoeffs_;

        //- Thermal conductivity coefficients
        coeffList kappaCoeffs_;


public:

    // Constructors

        //- Construct from components
        inline AndradeTransport
        (
            const Thermo& t,
            const coeffList& muCoeffs,
            const coeffList& kappaCoeffs
        );

        //- Construct as named copy
        inline AndradeTransport

        (
            const word&,
            const AndradeTransport&
        );

        //- Construct from name and dictionary
        AndradeTransport(const word& name, const dictionary& dict);

        //- Construct and return a clone
        inline autoPtr<AndradeTransport> clone() const;


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "Andrade<" + Thermo::typeName() + '>';
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
            const AndradeTransport&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "AndradeTransportI.H"

#ifdef NoRepository
    #include "AndradeTransport.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
