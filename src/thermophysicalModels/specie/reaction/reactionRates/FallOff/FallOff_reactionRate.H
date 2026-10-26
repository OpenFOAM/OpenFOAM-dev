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
    Foam::reactionRates::FallOff

Description
    General class for handling unimolecular/recombination fall-off reactions.

SourceFiles
    FallOff_reactionRateI.H
    FallOff_reactionRates.C

\*---------------------------------------------------------------------------*/

#ifndef FallOff_reactionRate_H
#define FallOff_reactionRate_H

#include "reaction.H"
#include "thirdBodyEfficiencies.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define FALL_OFF_FUNCTION_JACOBIAN false

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace reactionRates
{

/*---------------------------------------------------------------------------*\
                           Class FallOff Declaration
\*---------------------------------------------------------------------------*/

template<class Rate, class FallOffFunction>
class FallOff
{
    // Private Data

        //- Initial rate constant
        Rate k0_;

        //- Limiting rate constant
        Rate kInf_;

        //- Fall-off function
        FallOffFunction F_;

        //- Third body efficiencies
        thirdBodyEfficiencies thirdBodyEfficiencies_;


public:

    // Constructors

        //- Construct from components
        inline FallOff
        (
            const Rate& k0,
            const Rate& kInf,
            const FallOffFunction& F,
            const thirdBodyEfficiencies& tbes
        );

        //- Construct from dictionary
        inline FallOff
        (
            const speciesTable& species,
            const dimensionSet& dims,
            const dictionary& dict
        );

        //- Construct for objectRegistry
        inline FallOff
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
            return word
            (
                "FallOff<"
              + Rate::type()
              + ','
              + FallOffFunction::type()
              + '>'
            );
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
        static const bool hasDdc = FALL_OFF_FUNCTION_JACOBIAN && Rate::hasDdc;

        //- The derivative of the rate w.r.t. concentration
        inline void ddc(const reaction::state& s, scalarField& ddc) const;

        //- Write to stream
        inline void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace reactionRates

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Write to stream
template<class Rate, class FallOffFunction>
inline Ostream& operator<<
(
    Ostream&,
    const reactionRates::FallOff<Rate, FallOffFunction>&
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "FallOff_reactionRateI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
