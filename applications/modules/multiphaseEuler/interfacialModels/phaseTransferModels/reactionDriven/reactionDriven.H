/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019-2024 OpenFOAM Foundation
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
    Foam::phaseTransferModels::reactionDriven

Description
    Phase transfer model representing change from one phase to another due
    to reactions. Intended for irreversible reactions.

SourceFiles
    reactionDriven.C

\*---------------------------------------------------------------------------*/

#ifndef reactionDriven_H
#define reactionDriven_H

#include "phaseTransferModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace phaseTransferModels
{

/*---------------------------------------------------------------------------*\
                       Class reactionDriven Declaration
\*---------------------------------------------------------------------------*/

class reactionDriven
:
    public phaseTransferModel
{
private:

    // Private Data

        //- Interface
        const phaseInterface interface_;

        //- List of species transferring out of phase 1
        const hashedWordList species1_;

        //- List of species transferring out of phase 2
        const hashedWordList species2_;

        //- List of all transferring species
        hashedWordList species_;


    // Private Member Functions

        //- Keyword for the list of species
        template<bool Index>
        word speciesKey() const;


public:

    //- Runtime type information
    TypeName("reactionDriven");


    // Constructors

        //- Construct from a dictionary and an interface
        reactionDriven
        (
            const dictionary& dict,
            const phaseInterface& interface
        );


    //- Destructor
    virtual ~reactionDriven();


    // Member Functions

        //- The list of individual species that are transferred
        virtual const hashedWordList& species() const;

        //- The mass transfer rate for individual species
        virtual HashPtrTable<volScalarField> dmidtf() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace phaseTransferModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
