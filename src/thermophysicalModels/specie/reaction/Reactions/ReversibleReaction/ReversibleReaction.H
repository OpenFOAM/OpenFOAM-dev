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
    Foam::ReversibleReaction

Description
    Simple extension of Reaction to handle reversible reactions using
    equilibrium thermodynamics.

SourceFiles
    ReversibleReaction.C

\*---------------------------------------------------------------------------*/

#ifndef ReversibleReaction_H
#define ReversibleReaction_H

#include "Reaction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class ReversibleReaction Declaration
\*---------------------------------------------------------------------------*/

template
<
    template<class> class ReactionType,
    class ReactionThermo,
    class ReactionRate
>
class ReversibleReaction
:
    public ReactionType<ReactionThermo>
{
    // Private data

        ReactionRate k_;


    // Private Member Functions

        //- Disallow default bitwise assignment
        void operator=
        (
            const ReversibleReaction
            <
                ReactionType,
                ReactionThermo,
                ReactionRate
            >&
        );


public:

    //- Runtime type information
    TypeName("reversible");


    // Constructors

        //- Construct from components
        ReversibleReaction
        (
            const ReactionType<ReactionThermo>& reaction,
            const ReactionRate& k
        );

        //- Construct as copy given new speciesTable
        ReversibleReaction
        (
            const ReversibleReaction
            <
                ReactionType,
                ReactionThermo,
                ReactionRate
            >&,
            const speciesTable& species
        );

        //- Construct from dictionary
        ReversibleReaction
        (
            const speciesTable& species,
            const HashPtrTable<ReactionThermo>& thermoDatabase,
            const dictionary& dict
        );

        //- Construct and return a clone
        virtual autoPtr<ReactionType<ReactionThermo>> clone() const
        {
            return autoPtr<ReactionType<ReactionThermo>>
            (
                new ReversibleReaction
                <
                    ReactionType,
                    ReactionThermo,
                    ReactionRate
                >(*this)
            );
        }

        //- Construct and return a clone with new speciesTable
        virtual autoPtr<ReactionType<ReactionThermo>> clone
        (
            const speciesTable& species
        ) const
        {
            return autoPtr<ReactionType<ReactionThermo>>
            (
                new ReversibleReaction
                <ReactionType, ReactionThermo, ReactionRate>
                (
                    *this,
                    species
                )
            );
        }


    //- Destructor
    virtual ~ReversibleReaction()
    {}


    // Member Functions

        // ReversibleReaction rate coefficients

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
            //  Note this evaluates the forward rate constant and divides by
            //  the equilibrium constant
            virtual scalar kr
            (
                const scalar p,
                const scalar T,
                const scalarField& c
            ) const;


        // ReversibleReaction Jacobian functions

            //- Temperature derivative of forward rate
            virtual scalar dkfdT
            (
                const scalar p,
                const scalar T,
                const scalarField& c
            ) const;

            //- Temperature derivative of backward rate
            virtual scalar dkrdT
            (
                const scalar p,
                const scalar T,
                const scalarField& c,
                const scalar dkfdT,
                const scalar kr
            ) const;

            //- Third-body efficiencies (beta = 1-alpha)
            //  non-empty only for third-body reactions
            //  with enhanced molecularity (alpha != 1)
            virtual const List<Tuple2<label, scalar>>& beta() const;

            //- Species concentration derivative of the pressure dependent term
            //  By default this value is 1 as it multiplies the third-body term
            virtual void dcidc
            (
                const scalar p,
                const scalar T,
                const scalarField& c,
                scalarField& dcidc
            ) const;

            //- Temperature derivative of the pressure dependent term
            //  By default this value is 0 since ddT of molecularity is approx.0
            virtual scalar dcidT
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
    #include "ReversibleReaction.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
