/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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
    Foam::NonEquilibriumReversibleReaction

Description
    Simple extension of Reaction to handle reversible reactions using
    equilibrium thermodynamics.

SourceFiles
    NonEquilibriumReversibleReaction.C

\*---------------------------------------------------------------------------*/

#ifndef NonEquilibriumReversibleReaction_H
#define NonEquilibriumReversibleReaction_H

#include "Reaction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class NonEquilibriumReversibleReaction Declaration
\*---------------------------------------------------------------------------*/

template
<
    template<class> class ReactionType,
    class ReactionThermo,
    class ReactionRate
>
class NonEquilibriumReversibleReaction
:
    public ReactionType<ReactionThermo>
{
    // Private data

        ReactionRate fk_;
        ReactionRate rk_;


    // Private Member Functions

        //- Disallow default bitwise assignment
        void operator=
        (
            const NonEquilibriumReversibleReaction
                <ReactionType, ReactionThermo, ReactionRate>&
        );


public:

    //- Runtime type information
    TypeName("nonEquilibriumReversible");


    // Constructors

        //- Construct from components
        NonEquilibriumReversibleReaction
        (
            const ReactionType<ReactionThermo>& reaction,
            const ReactionRate& forwardReactionRate,
            const ReactionRate& reverseReactionRate
        );

        //- Construct as copy given new speciesTable
        NonEquilibriumReversibleReaction
        (
            const NonEquilibriumReversibleReaction
                <ReactionType, ReactionThermo, ReactionRate>&,
            const speciesTable& species
        );

        //- Construct from Istream
        NonEquilibriumReversibleReaction
        (
            const speciesTable& species,
            const HashPtrTable<ReactionThermo>& thermoDatabase,
            Istream& is
        );

        //- Construct from dictionary
        NonEquilibriumReversibleReaction
        (
            const speciesTable& species,
            const HashPtrTable<ReactionThermo>& thermoDatabase,
            const dictionary& dict
        );

        //- Construct and return a clone
        virtual autoPtr<ReactionType<ReactionThermo> > clone() const
        {
            return autoPtr<ReactionType<ReactionThermo> >
            (
                new NonEquilibriumReversibleReaction
                    <ReactionType, ReactionThermo, ReactionRate>(*this)
            );
        }

        //- Construct and return a clone with new speciesTable
        virtual autoPtr<ReactionType<ReactionThermo> > clone
        (
            const speciesTable& species
        ) const
        {
            return autoPtr<ReactionType<ReactionThermo> >
            (
                new NonEquilibriumReversibleReaction
                    <ReactionType, ReactionThermo, ReactionRate>
                    (*this, species)
            );
        }


    //- Destructor
    virtual ~NonEquilibriumReversibleReaction()
    {}


    // Member Functions

        // NonEquilibriumReversibleReaction rate coefficients

            //- Forward rate constant
            virtual scalar kf
            (
                const scalar p,
                const scalar T,
                const scalarField& c
            ) const;

            //- Reverse rate constant from the given formard rate constant
            virtual scalar kr
            (
                const scalar kfwd,
                const scalar p,
                const scalar T,
                const scalarField& c
            ) const;

            //- Reverse rate constant.
            //  Note this evaluates the forward rate constant and divides by the
            //  equilibrium constant
            virtual scalar kr
            (
                const scalar p,
                const scalar T,
                const scalarField& c
            ) const;


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "NonEquilibriumReversibleReaction.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
