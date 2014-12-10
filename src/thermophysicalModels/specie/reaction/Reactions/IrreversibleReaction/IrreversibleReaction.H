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
    Foam::IrreversibleReaction

Description
    Simple extension of Reaction to handle reversible reactions using
    equilibrium thermodynamics.

SourceFiles
    IrreversibleReaction.C

\*---------------------------------------------------------------------------*/

#ifndef IrreversibleReaction_H
#define IrreversibleReaction_H

#include "Reaction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class IrreversibleReaction Declaration
\*---------------------------------------------------------------------------*/

template
<
    template<class> class ReactionType,
    class ReactionThermo,
    class ReactionRate
>
class IrreversibleReaction
:
    public ReactionType<ReactionThermo>
{
    // Private data

        ReactionRate k_;


    // Private Member Functions

        //- Disallow default bitwise assignment
        void operator=
        (
            const IrreversibleReaction
            <
                ReactionType,
                ReactionThermo,
                ReactionRate
            >&
        );


public:

    //- Runtime type information
    TypeName("irreversible");


    // Constructors

        //- Construct from components
        IrreversibleReaction
        (
            const ReactionType<ReactionThermo>& reaction,
            const ReactionRate& reactionRate
        );

        //- Construct as copy given new speciesTable
        IrreversibleReaction
        (
            const IrreversibleReaction
            <
                ReactionType,
                ReactionThermo,
                ReactionRate
            >&,
            const speciesTable& species
        );

        //- Construct from Istream
        IrreversibleReaction
        (
            const speciesTable& species,
            const HashPtrTable<ReactionThermo>& thermoDatabase,
            Istream& is
        );

        //- Construct from dictionary
        IrreversibleReaction
        (
            const speciesTable& species,
            const HashPtrTable<ReactionThermo>& thermoDatabase,
            const dictionary& dict
        );

        //- Construct and return a clone
        virtual autoPtr<Reaction<ReactionThermo> > clone() const
        {
            return autoPtr<Reaction<ReactionThermo> >
            (
                new IrreversibleReaction
                <
                    ReactionType,
                    ReactionThermo,
                    ReactionRate
                >(*this)
            );
        }

        //- Construct and return a clone with new speciesTable
        virtual autoPtr<Reaction<ReactionThermo> > clone
        (
            const speciesTable& species
        ) const
        {
            return autoPtr<Reaction<ReactionThermo> >
            (
                new IrreversibleReaction
                <
                    ReactionType,
                    ReactionThermo,
                    ReactionRate
                >
                (
                    *this,
                    species
                )
            );
        }


    //- Destructor
    virtual ~IrreversibleReaction()
    {}


    // Member Functions

        // IrreversibleReaction rate coefficients

            //- Forward rate constant
            virtual scalar kf
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
#   include "IrreversibleReaction.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
