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
    Foam::pLogReactionRate

Description
    Pressure-dependent reaction rate coefficient evaluated by logarithmic
    interpolation of Arrhenius rate coefficients evaluated from a pressure
    dependent table of parameters.

SourceFiles
    pLogReactionRateI.H

\*---------------------------------------------------------------------------*/

#ifndef pLogReactionRate_H
#define pLogReactionRate_H

#include "speciesTable.H"
#include "scalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class pLogReactionRate Declaration
\*---------------------------------------------------------------------------*/

class pLogReactionRate
{
public:

    //- Coefficients for the specified pressure
    struct pCoeffs
    {
        scalar p;
        scalar A;
        scalar beta;
        scalar Ta;

        scalar lnp;
        scalar lnA;

        inline static const word typeName = "coeffs";
    };


private:

    // Private Data

        //- List of the coefficients for the specified pressures
        List<pCoeffs> coeffs_;


    // Private Member Functions

        //- Return the Arrhenius rate coefficient
        //  for given pressure coefficients
        inline scalar ki(const pCoeffs& coeffs, const scalar T) const;

        //- Return the log of the Arrhenius rate coefficient
        // for given pressure coefficients
        inline scalar lnki(const pCoeffs& coeffs, const scalar T) const;

        //- Return the temperature derivative of the Arrhenius rate coefficient
        // for given pressure coefficients
        inline scalar ddTi(const pCoeffs& coeffs, const scalar T) const;

        //- Return the pressure interpolation fraction for the given pressure
        inline scalar fp(const scalar p, const label pi) const;


public:

    // Constructors

        //- Construct from coefficients
        inline pLogReactionRate(const List<pCoeffs>& coeffs);

        //- Construct from dictionary
        inline pLogReactionRate
        (
            const speciesTable& species,
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

        //- Return the rate
        inline scalar operator()
        (
            const scalar p,
            const scalar T,
            const scalarField& c,
            const label li
        ) const;

        //- The derivative of the rate w.r.t. temperature
        inline scalar ddT
        (
            const scalar p,
            const scalar T,
            const scalarField& c,
            const label li
        ) const;

        //- Is the rate a function of concentration?
        inline bool hasDdc() const;

        //- The derivative of the rate w.r.t. concentration
        inline void ddc
        (
            const scalar p,
            const scalar T,
            const scalarField& c,
            const label li,
            scalarField& ddc
        ) const;

        //- Write to stream
        inline void write(Ostream& os) const;
};


inline bool operator==
(
    const pLogReactionRate::pCoeffs& ce1,
    const pLogReactionRate::pCoeffs& ce2
);

inline bool operator!=
(
    const pLogReactionRate::pCoeffs& ce1,
    const pLogReactionRate::pCoeffs& ce2
);

inline Istream& operator>>(Istream& is, pLogReactionRate::pCoeffs&);
inline Ostream& operator<<(Ostream& os, const pLogReactionRate::pCoeffs&);

inline Ostream& operator<<(Ostream&, const pLogReactionRate&);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "pLogReactionRateI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
