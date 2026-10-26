/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2026 OpenFOAM Foundation
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
    Foam::reactionRates::Arrhenius

Description
    Arrhenius reaction rate given by:

        k = A * T^beta * exp(-Ta/T)

SourceFiles
    Arrhenius_reactionRateI.H
    Arrhenius_reactionRates.C

\*---------------------------------------------------------------------------*/

#ifndef Arrhenius_reactionRate_H
#define Arrhenius_reactionRate_H

#include "reaction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace reactionRates
{

/*---------------------------------------------------------------------------*\
                          Class Arrhenius Declaration
\*---------------------------------------------------------------------------*/

class Arrhenius
{
    // Private Data

        //- Temperature exponent
        scalar beta_;

        //- Pre-exponential factor
        scalar A_;

        //- Activation temperature
        scalar Ta_;


public:

    // Constructors

        //- Construct from components
        inline Arrhenius
        (
            const scalar A,
            const scalar beta,
            const scalar Ta
        );

        //- Construct from dictionary
        inline Arrhenius
        (
            const speciesTable& species,
            const dimensionSet& dims,
            const dictionary& dict
        );

        //- Construct for objectRegistry
        inline Arrhenius
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
            return "Arrhenius";
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
        static const bool hasDdc = false;

        //- The derivative of the rate w.r.t. concentration
        inline void ddc(const reaction::state& s, scalarField& ddc) const;

        //- Write to stream
        inline void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace reactionRates

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Write to stream
Ostream& operator<<(Ostream&, const reactionRates::Arrhenius&);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "Arrhenius_reactionRateI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
