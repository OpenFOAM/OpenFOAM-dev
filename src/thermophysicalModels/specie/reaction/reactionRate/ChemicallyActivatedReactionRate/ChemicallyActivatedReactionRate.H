/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    Foam::ChemicallyActivatedReactionRate

Description
    General class for handling chemically-activated bimolecular reactions.

SourceFiles
    ChemicallyActivatedReactionRateI.H

\*---------------------------------------------------------------------------*/

#ifndef ChemicallyActivatedReactionRate_H
#define ChemicallyActivatedReactionRate_H

#include "thirdBodyEfficiencies.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class ReactionRate, class ChemicallyActivationFunction>
class ChemicallyActivatedReactionRate;

template<class ReactionRate, class ChemicallyActivationFunction>
inline Ostream& operator<<
(
    Ostream&,
    const ChemicallyActivatedReactionRate
        <ReactionRate, ChemicallyActivationFunction>&
);


/*---------------------------------------------------------------------------*\
              Class ChemicallyActivatedReactionRate Declaration
\*---------------------------------------------------------------------------*/

template<class ReactionRate, class ChemicallyActivationFunction>
class ChemicallyActivatedReactionRate
{
    // Private data

        ReactionRate k0_;
        ReactionRate kInf_;
        ChemicallyActivationFunction F_;
        thirdBodyEfficiencies thirdBodyEfficiencies_;
        List<Tuple2<label, scalar>> beta_;


public:

    // Constructors

        //- Construct from components
        inline ChemicallyActivatedReactionRate
        (
            const ReactionRate& k0,
            const ReactionRate& kInf,
            const ChemicallyActivationFunction& F,
            const thirdBodyEfficiencies& tbes
        );

        //- Construct from dictionary
        inline ChemicallyActivatedReactionRate
        (
            const speciesTable& species,
            const dictionary& dict
        );


    // Member Functions

        //- Return the type name
        static word type()
        {
            return ReactionRate::type()
                + ChemicallyActivationFunction::type()
                + "ChemicallyActivated";
        }

        inline scalar operator()
        (
            const scalar p,
            const scalar T,
            const scalarField& c
        ) const;

        inline scalar ddT
        (
            const scalar p,
            const scalar T,
            const scalarField& c
        ) const;

        //- Third-body efficiencies (beta = 1-alpha)
        //  non-empty only for third-body reactions
        //  with enhanced molecularity (alpha != 1)
        inline const List<Tuple2<label, scalar>>& beta() const
        {
            return beta_;
        }

        //- Species concentration derivative of the pressure dependent term
        //  By default this value is 1 as it multiplies the third-body term
        inline void dcidc
        (
            const scalar p,
            const scalar T,
            const scalarField& c,
            scalarField& dcidc
        ) const;

        //- Temperature derivative of the pressure dependent term
        //  By default this value is 0 since ddT of molecularity is approx.0
        inline scalar dcidT
        (
            const scalar p,
            const scalar T,
            const scalarField& c
        ) const;

        //- Write to stream
        inline void write(Ostream& os) const;


    // Ostream Operator

        friend Ostream& operator<< <ReactionRate, ChemicallyActivationFunction>
        (
            Ostream&,
            const ChemicallyActivatedReactionRate
                <ReactionRate, ChemicallyActivationFunction>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ChemicallyActivatedReactionRateI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
