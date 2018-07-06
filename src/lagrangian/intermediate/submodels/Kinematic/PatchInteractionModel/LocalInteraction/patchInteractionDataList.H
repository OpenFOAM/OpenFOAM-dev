/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
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
    Foam::patchInteractionDataList

Description
    List container for patchInteractionData class

\*---------------------------------------------------------------------------*/

#ifndef patchInteractionDataList_H
#define patchInteractionDataList_H

#include "patchInteractionData.H"
#include "polyMesh.H"
#include "dictionary.H"
#include "labelList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                  Class patchInteractionDataList Declaration
\*---------------------------------------------------------------------------*/

class patchInteractionDataList
:
    public List<patchInteractionData>
{
private:

    // Private data

        //- List of patch IDs for each patch group
        labelListList patchGroupIDs_;


public:

    // Constructor

        //- Construct null
        patchInteractionDataList();

        //- Construct copy
        patchInteractionDataList(const patchInteractionDataList& pidl);

        //- Construct from Istream
        patchInteractionDataList(const polyMesh& mesh, const dictionary& dict);


    // Member functions

        //- Return label of group containing patch id
        //  Returns -1 if patch id is not present
        label applyToPatch(const label id) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
