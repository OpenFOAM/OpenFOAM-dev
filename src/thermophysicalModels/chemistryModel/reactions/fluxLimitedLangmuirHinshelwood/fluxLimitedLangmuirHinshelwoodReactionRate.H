/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019-2021 OpenFOAM Foundation
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
    Foam::fluxLimitedLangmuirHinshelwoodReactionRate

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
    fluxLimitedLangmuirHinshelwoodReactionRateI.H

\*---------------------------------------------------------------------------*/

#ifndef fluxLimitedLangmuirHinshelwoodReactionRate_H
#define fluxLimitedLangmuirHinshelwoodReactionRate_H

#include "scalarField.H"
#include "typeInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

class fluxLimitedLangmuirHinshelwoodReactionRate;

Ostream& operator<<
(
    Ostream&,
    const fluxLimitedLangmuirHinshelwoodReactionRate&
);


/*---------------------------------------------------------------------------*\
         Class fluxLimitedLangmuirHinshelwoodReactionRate Declaration
\*---------------------------------------------------------------------------*/

class fluxLimitedLangmuirHinshelwoodReactionRate
{
    // Private data

        //- Number reactants (1 or 2)
        label nReactants_;

        //- Additional adsorbable specie names
        wordList additionalAdsorbableSpecieNames_;

        //- Reactant and additional adsorbable specie indices
        labelList ra_;

        scalar a_;
        scalarList A_;
        scalarList Ta_;
        scalarList beta_;
        scalarList m_;

        bool limited_;
        bool AvUniform_;
        scalar Av_;
        const word AvName_;
        const scalarField& AvField_;
        scalarList nu_;
        scalarList exp_;
        scalarList s_;
        scalarList W_;


    // Private Member Functions

        //- Return either the uniform Av value or the field element
        inline scalar Av(const label li) const;


public:

    // Constructors

        //- Construct from dictionary
        inline fluxLimitedLangmuirHinshelwoodReactionRate
        (
            const speciesTable& species,
            const objectRegistry& ob,
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

        inline scalar operator()
        (
            const scalar p,
            const scalar T,
            const scalarField& c,
            const label li
        ) const;

        inline scalar ddT
        (
            const scalar p,
            const scalar T,
            const scalarField& c,
            const label li
        ) const;

        inline bool hasDdc() const;

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


    // Ostream Operator

        inline friend Ostream& operator<<
        (
            Ostream&,
            const fluxLimitedLangmuirHinshelwoodReactionRate&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "fluxLimitedLangmuirHinshelwoodReactionRateI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
