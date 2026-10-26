/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2026 OpenFOAM Foundation
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
    Foam::reactionRates::PhaseSurface

Description
    Reaction rate modified by multiplying by a phase's surface area per unit
    volume, with the name of the phase specified by the user.

SourceFiles
    PhaseSurface_reactionRateI.H
    PhaseSurface_reactionRates.C

\*---------------------------------------------------------------------------*/

#ifndef PhaseSurface_reactionRate_H
#define PhaseSurface_reactionRate_H

#include "Arrhenius_reactionRate.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace reactionRates
{

/*---------------------------------------------------------------------------*\
                      Class PhaseSurface Declaration
\*---------------------------------------------------------------------------*/

template<class Rate>
class PhaseSurface
{
    // Private Data

        //- Unmodified rate constant
        const Rate k_;

        //- Name of the phase
        const word phaseName_;

        //- Reference to the object registry
        const objectRegistry& ob_;

        //- Surface area per unit volume field
        mutable tmp<volScalarField> tAv_;


public:

    // Constructors

        //- Construct for objectRegistry
        inline PhaseSurface
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
            return word("PhaseSurface<", Rate::type(), '>');
        }

        //- Pre-evaluation hook
        inline void preEvaluate() const;

        //- Post-evaluation hook
        inline void postEvaluate() const;

        //- Evaluate the rate constant
        inline scalar k(const reaction::state& s) const;

        //- Evaluate the derivative
        inline scalar ddT(const reaction::state& s) const;

        //- Is the rate a function of concentration?
        static const bool hasDdc = Rate::hasDdc;

        //- The derivative of the rate w.r.t. concentration
        inline void ddc(const reaction::state& s, scalarField& ddc) const;

        //- Write to stream
        inline void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace reactionRates

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Write to stream
template<class Rate>
Ostream& operator<<(Ostream&, const reactionRates::PhaseSurface<Rate>&);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "PhaseSurface_reactionRateI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
