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
    Foam::reactionRates::pLog

Description
    Pressure-dependent reaction rate coefficient evaluated by logarithmic
    interpolation of Arrhenius rate coefficients evaluated from a pressure
    dependent table of parameters.

SourceFiles
    pLog_reactionRateI.H
    pLog_reactionRates.C

\*---------------------------------------------------------------------------*/

#ifndef pLog_reactionRate_H
#define pLog_reactionRate_H

#include "reaction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace reactionRates
{

/*---------------------------------------------------------------------------*\
                             Class pLog Declaration
\*---------------------------------------------------------------------------*/

class pLog
{
public:

    // Public Classes

        //- Coefficients for the specified pressure
        struct pCoeffs
        {
            //- Type name needed for list IO
            inline static const word typeName = "pLog::pCoeffs";

            //- Pressure
            scalar p;

            //- Pre-exponential factor
            scalar A;

            //- Temperature exponent
            scalar beta;

            //- Activation temperature
            scalar Ta;

            //- Log of the pressure
            scalar lnp;

            //- Log of the pre-exponential factor
            scalar lnA;
        };


private:

    // Private Data

        //- List of the coefficients for the specified pressures
        List<pCoeffs> coeffs_;


    // Private Member Functions

        //- Return the Arrhenius rate coefficient
        //  for given pressure coefficients
        inline scalar ki
        (
            const pCoeffs& coeffs,
            const reaction::state& s
        ) const;

        //- Return the log of the Arrhenius rate coefficient
        // for given pressure coefficients
        inline scalar lnki
        (
            const pCoeffs& coeffs,
            const reaction::state& s
        ) const;

        //- Return the temperature derivative of the Arrhenius rate coefficient
        // for given pressure coefficients
        inline scalar ddTi
        (
            const pCoeffs& coeffs,
            const reaction::state& s
        ) const;

        //- Return the pressure interpolation fraction for the given pressure
        inline scalar fp(const scalar lnp, const label pi) const;


public:

    // Constructors

        //- Construct from coefficients
        inline pLog(const List<pCoeffs>& coeffs);

        //- Construct from dictionary
        inline pLog
        (
            const speciesTable& species,
            const dimensionSet& dims,
            const dictionary& dict
        );

        //- Construct for objectRegistry
        inline pLog
        (
            const speciesTable& species,
            const objectRegistry& ob,
            const dimensionSet& dims,
            const dictionary& dict
        );


    // Member Functions

        //- Return the type name
        static word type()
        {
            return "pLog";
        }

        //- Pre-evaluation hook
        inline void preEvaluate() const;

        //- Post-evaluation hook
        inline void postEvaluate() const;

        //- Return the rate constant
        inline scalar k(const reaction::state& s) const;

        //- The derivative of the rate w.r.t. temperature
        inline scalar ddT(const reaction::state& s) const;

        //- Is the rate a function of concentration?
        static const bool hasDdc = true;

        //- The derivative of the rate w.r.t. concentration
        inline void ddc(const reaction::state& s, scalarField& ddc) const;

        //- Write to stream
        inline void write(Ostream& os) const;
};


//- Compare coeffs
inline bool operator==
(
    const pLog::pCoeffs& ce1,
    const pLog::pCoeffs& ce2
);

//- Compare coeffs
inline bool operator!=
(
    const pLog::pCoeffs& ce1,
    const pLog::pCoeffs& ce2
);

//- Read coeffs from stream
inline Istream& operator>>(Istream& is, pLog::pCoeffs&);

//- Write coeffs to stream
inline Ostream& operator<<(Ostream& os, const pLog::pCoeffs&);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace reactionRates

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Write to stream
inline Ostream& operator<<(Ostream&, const reactionRates::pLog&);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "pLog_reactionRateI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
