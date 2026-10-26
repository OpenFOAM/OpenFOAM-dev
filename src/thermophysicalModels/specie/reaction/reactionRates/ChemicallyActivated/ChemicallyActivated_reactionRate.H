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
    Foam::reactionRates::ChemicallyActivated

Description
    General class for handling chemically-activated bimolecular reactions.

SourceFiles
    ChemicallyActivated_reactionRateI.H
    ChemicallyActivated_reactionRates.C

\*---------------------------------------------------------------------------*/

#ifndef ChemicallyActivated_reactionRate_H
#define ChemicallyActivated_reactionRate_H

#include "reaction.H"
#include "thirdBodyEfficiencies.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define CHEMICAL_ACTIVATION_FUNCTION_JACOBIAN false

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace reactionRates
{

/*---------------------------------------------------------------------------*\
                      Class ChemicallyActivated Declaration
\*---------------------------------------------------------------------------*/

template<class Rate, class ChemicalActivationFunction>
class ChemicallyActivated
{
    // Private Data

        //- Initial rate constant
        Rate k0_;

        //- Limiting rate constant
        Rate kInf_;

        //- Activation function
        ChemicalActivationFunction F_;

        //- Third body efficiencies
        thirdBodyEfficiencies thirdBodyEfficiencies_;


public:

    // Constructors

        //- Construct from components
        inline ChemicallyActivated
        (
            const Rate& k0,
            const Rate& kInf,
            const ChemicalActivationFunction& F,
            const thirdBodyEfficiencies& tbes
        );

        //- Construct from dictionary
        inline ChemicallyActivated
        (
            const speciesTable& species,
            const dimensionSet& dims,
            const dictionary& dict
        );

        //- Construct for objectRegistry
        inline ChemicallyActivated
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
                "ChemicallyActivated<",
                Rate::type(),
                ',',
                ChemicalActivationFunction::type(),
                '>'
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
        static const bool hasDdc =
            CHEMICAL_ACTIVATION_FUNCTION_JACOBIAN && Rate::hasDdc;

        //- The derivative of the rate w.r.t. concentration
        inline void ddc(const reaction::state& s, scalarField& ddc) const;

        //- Write to stream
        inline void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace reactionRates

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Write to stream
template<class Rate, class ChemicalActivationFunction>
inline Ostream& operator<<
(
    Ostream&,
    const reactionRates::ChemicallyActivated<Rate, ChemicalActivationFunction>&
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ChemicallyActivated_reactionRateI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
