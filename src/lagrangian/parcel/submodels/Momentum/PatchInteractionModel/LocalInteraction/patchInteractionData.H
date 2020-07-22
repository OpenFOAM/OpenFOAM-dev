/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
    Foam::patchInteractionData

Description
    Helper class for the LocalInteraction patch interaction model

\*---------------------------------------------------------------------------*/

#ifndef patchInteractionData_H
#define patchInteractionData_H

#include "Istream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                    Class patchInteractionData Declaration
\*---------------------------------------------------------------------------*/

// Forward declaration of classes
class patchInteractionData;

// Forward declaration of friend functions
Istream& operator>>
(
    Istream& is,
    patchInteractionData& pid
);


class patchInteractionData
{
    // Private Data

        //- Interaction type name
        word interactionTypeName_;

        //- Patch name
        word patchName_;

        //- Elasticity coefficient
        scalar e_;

        //- Restitution coefficient
        scalar mu_;


public:

    // Constructor

        //- Construct null
        patchInteractionData();


    // Member Functions

        // Access

            //- Return const access to the interaction type name
            const word& interactionTypeName() const;

            //- Return const access to the patch name
            const word& patchName() const;

            //- Return const access to the elasticity coefficient
            scalar e() const;

            //- Return const access to the restitution coefficient
            scalar mu() const;


        // I-O

            //- Istream operator
            friend Istream& operator>>
            (
                Istream& is,
                patchInteractionData& pid
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
