/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020-2026 OpenFOAM Foundation
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
    Foam::reactionCoeffs

Description
    Reaction base-class holding the specie names and coefficients

SourceFiles
    reactionCoeffsI.H
    reactionCoeffs.C

\*---------------------------------------------------------------------------*/

#ifndef reactionCoeffs_H
#define reactionCoeffs_H

#include "speciesTable.H"
#include "specieCoeffs.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class reactionCoeffs Declaration
\*---------------------------------------------------------------------------*/

class reactionCoeffs
{
    // Private Data

        //- Name of reaction
        const word name_;

        //- List of specie names present in reaction system
        const speciesTable& species_;

        //- Specie info for the left-hand-side of the reaction
        List<specieCoeffs> lhs_;

        //- Specie info for the right-hand-side of the reaction
        List<specieCoeffs> rhs_;

        //- Flat list of the indices of all the species in the reaction. This
        //  is useful sometimes when optimising what would otherwise be two
        //  loops over the left and right hand sides. Sometimes it is also
        //  necessary to avoid duplication in the case when the left and right
        //  hand sides contain terms referring to the same specie.
        List<label> specieIndices_;


    // Private Member Data

        //- Set the specie indices
        void setSpecieIndices();


public:

    // Constructors

        //- Construct from components
        reactionCoeffs
        (
            const word& name,
            const speciesTable& species,
            const List<specieCoeffs>& lhs,
            const List<specieCoeffs>& rhs
        );

        //- Construct from dictionary
        reactionCoeffs
        (
            const speciesTable& species,
            const dictionary& dict
        );


    //- Destructor
    ~reactionCoeffs()
    {}


    // Member Functions

        // Access

            //- Return the name of the reaction
            inline const word& name() const;

            //- Return the specie list
            inline const speciesTable& species() const;

            //- Return the components of the left hand side
            inline const List<specieCoeffs>& lhs() const;

            //- Return the components of the right hand side
            inline const List<specieCoeffs>& rhs() const;

            //- Return the species' indices
            inline const List<label>& specieIndices() const;


        //- Write
        void write(Ostream&) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const reactionCoeffs&) = delete;
};


//- Write to a stream
Ostream& operator<<(Ostream&, const reactionCoeffs&);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "reactionCoeffsI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
