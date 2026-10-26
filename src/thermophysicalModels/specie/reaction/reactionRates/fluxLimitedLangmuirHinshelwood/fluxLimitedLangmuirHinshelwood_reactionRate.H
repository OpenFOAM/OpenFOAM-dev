/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019-2026 OpenFOAM Foundation
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
    Foam::reactionRates::fluxLimitedLangmuirHinshelwood

Description
    Langmuir-Hinshelwood reaction rate for gaseous reactions on surfaces
    including the optional flux limiter of Waletzko and Schmidt.

    References:
    \verbatim
        Hinshelwood, C.N. (1940).
        The Kinetics of Chemical Change.
        Oxford Clarendon Press

        Waletzko, N., & Schmidt, L. D. (1988).
        Modeling catalytic gauze reactors: HCN synthesis.
        AIChE journal, 34(7), 1146-1156.
    \endverbatim

SourceFiles
    fluxLimitedLangmuirHinshelwood_reactionRateI.H
    fluxLimitedLangmuirHinshelwood_reactionRates.C

\*---------------------------------------------------------------------------*/

#ifndef fluxLimitedLangmuirHinshelwood_reactionRate_H
#define fluxLimitedLangmuirHinshelwood_reactionRate_H

#include "reaction.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace reactionRates
{

/*---------------------------------------------------------------------------*\
              Class fluxLimitedLangmuirHinshelwood Declaration
\*---------------------------------------------------------------------------*/

class fluxLimitedLangmuirHinshelwood
{
    // Private data

        //- Number reactants (1 or 2)
        label nReactants_;

        //- Additional adsorbable specie names
        wordList additionalAdsorbableSpecieNames_;

        //- Reactant and additional adsorbable specie indices
        labelList ra_;

        //- ?
        scalar a_;

        //- Pre-exponential factor
        scalarList A_;

        //- Activation temperature
        scalarList Ta_;

        //- Temperature exponent
        scalarList beta_;

        //- ?
        scalarList m_;

        //- Is the reaction limited?
        bool limited_;

        //- Is the surface area per unit volume a uniform value?
        bool AvUniform_;

        //- Surface area per unit volume uniform value
        scalar Av_;

        //- Surface area per unit volume field name
        const word AvName_;

        //- ?
        scalarList nu_;

        //- ?
        scalarList exp_;

        //- ?
        scalarList s_;

        //- ?
        scalarList W_;

        //- Reference to the object registry
        const objectRegistry& ob_;

        //- Surface area per unit volume field
        mutable tmp<volInternalScalarField> tAv_;


    // Private Member Functions

        //- Return either the uniform Av value or the field element
        inline scalar Av(const label li) const;


public:

    // Constructors

        //- Construct for objectRegistry
        inline fluxLimitedLangmuirHinshelwood
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
            return "fluxLimitedLangmuirHinshelwood";
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
Ostream& operator<<
(
    Ostream&,
    const reactionRates::fluxLimitedLangmuirHinshelwood&
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "fluxLimitedLangmuirHinshelwood_reactionRateI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
