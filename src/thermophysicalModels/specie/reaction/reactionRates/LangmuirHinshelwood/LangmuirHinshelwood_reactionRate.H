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
    Foam::reactionRates::LangmuirHinshelwood

Description
    Langmuir-Hinshelwood reaction rate for gaseous reactions on surfaces.

    Reference:
    \verbatim
        Hinshelwood, C.N. (1940).
        The Kinetics of Chemical Change.
        Oxford Clarendon Press
    \endverbatim

SourceFiles
    LangmuirHinshelwood_reactionRateI.H
    LangmuirHinshelwood_reactionRates.C

\*---------------------------------------------------------------------------*/

#ifndef LangmuirHinshelwood_reactionRate_H
#define LangmuirHinshelwood_reactionRate_H

#include "reaction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace reactionRates
{

/*---------------------------------------------------------------------------*\
                    Class LangmuirHinshelwood Declaration
\*---------------------------------------------------------------------------*/

class LangmuirHinshelwood
{
    // Private Data

        //- Reactant names
        FixedList<word, 2> reactantNames_;

        //- Reactant indices
        FixedList<label, 3> r_;

        scalar a_;

        FixedList<scalar, 3> A_;

        FixedList<scalar, 3> Ta_;

        FixedList<scalar, 3> beta_;

        FixedList<scalar, 3> m_;


public:

    // Constructors

        //- Construct from dictionary
        inline LangmuirHinshelwood
        (
            const speciesTable& species,
            const dimensionSet& dims,
            const dictionary& dict
        );

        //- Construct for objectRegistry
        inline LangmuirHinshelwood
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
            return "LangmuirHinshelwood";
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


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace reactionRates

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Write to stream
Ostream& operator<<(Ostream&, const reactionRates::LangmuirHinshelwood&);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "LangmuirHinshelwood_reactionRateI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
