/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020-2023 OpenFOAM Foundation
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
    Foam::reaction

Description
    Reaction base-class holding the specie names and coefficients

SourceFiles
    reactionI.H
    reaction.C

\*---------------------------------------------------------------------------*/

#ifndef reaction_H
#define reaction_H

#include "speciesTable.H"
#include "specieCoeffs.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

class reaction;

inline Ostream& operator<<(Ostream&, const reaction&);

/*---------------------------------------------------------------------------*\
                          Class reaction Declaration
\*---------------------------------------------------------------------------*/

class reaction
{

public:

    // Static data

        //- Number of un-named reactions
        static label nUnNamedReactions;


private:

    // Private Data

        //- Name of reaction
        const word name_;

        //- List of specie names present in reaction system
        const speciesTable& species_;

        //- Specie info for the left-hand-side of the reaction
        List<specieCoeffs> lhs_;

        //- Specie info for the right-hand-side of the reaction
        List<specieCoeffs> rhs_;


    // Private Member Functions

        //- Return new reaction ID for un-named reactions
        label getNewReactionIndex();


public:

    // Constructors

        //- Construct from components
        reaction
        (
            const speciesTable& species,
            const List<specieCoeffs>& lhs,
            const List<specieCoeffs>& rhs
        );

        //- Construct as copy given new speciesTable
        reaction(const reaction&, const speciesTable& species);

        //- Construct from dictionary
        reaction
        (
            const speciesTable& species,
            const dictionary& dict
        );


    //- Destructor
    ~reaction()
    {}


    // Member Functions

        // Access

            //- Return the name of the reaction
            inline const word& name() const;

            //- Return the components of the left hand side
            inline const List<specieCoeffs>& lhs() const;

            //- Return the components of the right hand side
            inline const List<specieCoeffs>& rhs() const;

            //- Return the specie list
            inline const speciesTable& species() const;


        //- Write
        void write(Ostream&) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const reaction&) = delete;


    // Ostream Operator

        friend Ostream& operator<<(Ostream&, const reaction&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "reactionI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
