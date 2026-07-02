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
    Foam::NSRDSTransport

Description
    Transport properties package using NSRDS function 101 for viscosity and
    NSRDS function 100 for thermal conductivity, both as a function of
    temperature.

    Reference:
    \verbatim
    Daubert, T. E., & Danner, R. P. (1985).
    Data compilation tables of properties of pure compounds.
    \endverbatim

Usage
    \table
        Property        | Description
        mu              | Dynamic viscosity coefficients
        kappa           | Thermal conductivity coefficients
    \endtable

    Example specification of NSRDS transport for water, with coefficients as
    given by Daubert & Danner:
    \verbatim
    transport
    {
        mu
        {
            a               -51.964;
            b               3670.6;
            c               5.7331;
            d               -5.3495e-29;
            e               10;
        }
        kappa
        {
            a               -0.4267;
            b               0.0056903;
            c               -8.0065e-06;
            d               1.815e-09;
            e               0;
        }
    }
    \endverbatim

SourceFiles
    NSRDSTransportI.H
    NSRDSTransport.C

See also
    Foam::Function1s::NonUniformTable

\*---------------------------------------------------------------------------*/

#ifndef NSRDSTransport_H
#define NSRDSTransport_H

#include "scalar.H"
#include "FixedPolynomial.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Thermo> class NSRDSTransport;

template<class Thermo>
Ostream& operator<<(Ostream&, const NSRDSTransport<Thermo>&);


/*---------------------------------------------------------------------------*\
                     Class NSRDSTransport Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo>
class NSRDSTransport
:
    public Thermo
{
    // Private Data

        //- Dynamic viscosity [Pa.s]
        scalar muA_, muB_, muC_, muD_, muE_;

        //- Thermal conductivity [W/m/K]
        FixedPolynomial<scalar, 5> kappaCoeffs_;


    // Private Constructors

        //- Construct from name and parent dictionary and dictionary
        NSRDSTransport
        (
            const word& name,
            const dictionary& dict,
            const dictionary& muDict,
            const dictionary& kappaDict
        );


public:

    // Constructors

        //- Construct as named copy
        inline NSRDSTransport(const word&, const NSRDSTransport&);

        //- Construct from name and dictionary
        NSRDSTransport(const word& name, const dictionary& dict);

        //- Construct and return a clone
        inline autoPtr<NSRDSTransport> clone() const;


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "NSRDS<" + Thermo::typeName() + '>';
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
            const NSRDSTransport&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "NSRDSTransportI.H"

#ifdef NoRepository
    #include "NSRDSTransport.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
