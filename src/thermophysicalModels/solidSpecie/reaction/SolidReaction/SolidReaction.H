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
    Foam::SolidReaction

Description

    Read solid reactions of the type S1 = S2 + G1

SourceFiles
    SolidReactionI.H
    SolidReaction.C

\*---------------------------------------------------------------------------*/

#ifndef SolidReaction_H
#define SolidReaction_H

#include "speciesTable.H"
#include "Reaction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators
template<class ReactionThermo>
class SolidReaction;

template<class ReactionThermo>
inline Ostream& operator<<(Ostream&, const SolidReaction<ReactionThermo>&);

/*---------------------------------------------------------------------------*\
                        Class SolidReaction Declaration
\*---------------------------------------------------------------------------*/

template<class ReactionThermo>
class SolidReaction
:
    public Reaction<ReactionThermo>
{

private:

    // Private data

        typedef typename Reaction<ReactionThermo>::specieCoeffs specieCoeffs;

        //- List of gas species present in reaction system
        speciesTable pyrolisisGases_;

        //- Gas specie index for the left-hand-side of the reaction
        List<specieCoeffs> glhs_;

        //- Gas specie index for the right-hand-side of the reaction
        List<specieCoeffs> grhs_;


    // Private Member Functions


        //- Return string representation of reaction
        string solidReactionStr(OStringStream&) const;

        //- Return string representation of the left of the reaction
        void solidReactionStrLeft(OStringStream&) const;

        //- Return string representation of the right of the reaction
        void solidReactionStrRight(OStringStream&) const;

        //- Disallow default bitwise assignment
        void operator=(const SolidReaction&);


public:

    //- Runtime type information
    TypeName("SolidReaction");


    // Constructors

        //- Construct from components
        SolidReaction
        (
            const Reaction<ReactionThermo>& reaction,
            const speciesTable& pyrolisisGases,
            const List<specieCoeffs>& glhs,
            const List<specieCoeffs>& grhs
        );


        //- Construct as copy given new speciesTable
        SolidReaction
        (
            const SolidReaction<ReactionThermo>&,
            const speciesTable& pyrolisisGases
        );


        //- Construct from dictionary
        SolidReaction
        (
            const speciesTable& species,
            const HashPtrTable<ReactionThermo>& thermoDatabase,
            const dictionary& dict
        );


    //- Destructor
    virtual ~SolidReaction()
    {}


    // Member Functions

        // Access

            //- Access to the gas components of the left hand side
            virtual const List<specieCoeffs>& grhs() const;

            //- Access to the gas components of the right hand side
            virtual const List<specieCoeffs>& glhs() const;

            //- Access to gas specie list
            virtual const speciesTable& gasSpecies() const;


        //- Write
        virtual void write(Ostream&) const;


    // Ostream Operator

        friend Ostream& operator<< <ReactionThermo>
        (
            Ostream&,
            const SolidReaction<ReactionThermo>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "SolidReactionI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "SolidReaction.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
