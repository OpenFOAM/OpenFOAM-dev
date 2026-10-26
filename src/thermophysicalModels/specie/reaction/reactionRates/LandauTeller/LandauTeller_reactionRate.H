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
    Foam::reactionRates::LandauTeller

Description
    Landau-Teller reaction rate.

SourceFiles
    LandauTeller_reactionRateI.H
    LandauTeller_reactionRates.C

\*---------------------------------------------------------------------------*/

#ifndef LandauTeller_reactionRate_H
#define LandauTeller_reactionRate_H

#include "reaction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace reactionRates
{

/*---------------------------------------------------------------------------*\
                  Class LandauTeller Declaration
\*---------------------------------------------------------------------------*/

class LandauTeller
{
    // Private Data

        scalar beta_;

        scalar A_;

        scalar Ta_;

        scalar B_;

        scalar C_;


public:

    // Constructors

        //- Construct from components
        inline LandauTeller
        (
            const scalar A,
            const scalar beta,
            const scalar Ta,
            const scalar B,
            const scalar C
        );

        //- Construct from dictionary
        inline LandauTeller
        (
            const speciesTable& species,
            const dimensionSet& dims,
            const dictionary& dict
        );

        //- Construct for objectRegistry
        inline LandauTeller
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
            return "LandauTeller";
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
Ostream& operator<<(Ostream&, const reactionRates::LandauTeller&);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "LandauTeller_reactionRateI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
