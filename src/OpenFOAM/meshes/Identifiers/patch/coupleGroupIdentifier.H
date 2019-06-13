/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2019 OpenFOAM Foundation
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
    Foam::coupleGroupIdentifier

Description
    Encapsulates using patchGroups to specify coupled patch

SourceFiles
    coupleGroupIdentifierI.H
    coupleGroupIdentifier.C
    coupleGroupIdentifierIO.C

\*---------------------------------------------------------------------------*/

#ifndef coupleGroupIdentifier_H
#define coupleGroupIdentifier_H

#include "word.H"
#include "label.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class dictionary;
class polyMesh;
class polyPatch;
class Ostream;

// Forward declaration of friend functions and operators
class coupleGroupIdentifier;
Ostream& operator<<(Ostream&, const coupleGroupIdentifier&);


/*---------------------------------------------------------------------------*\
                         Class coupleGroupIdentifier Declaration
\*---------------------------------------------------------------------------*/

class coupleGroupIdentifier
{
    // Private Data

        //- Name of patchGroup
        word name_;


    // Private Member Functions

        //- Find other patch in specified mesh. Returns index of patch or -1.
        label findOtherPatchID(const polyMesh&, const polyPatch&) const;


public:

    // Constructors

        //- Construct null
        coupleGroupIdentifier();

        //- Construct from components
        coupleGroupIdentifier(const word& patchGroupName);

        //- Construct from dictionary
        coupleGroupIdentifier(const dictionary&);


    // Member Functions

        //- Name of patchGroup
        inline const word& name() const;

        //- Is a valid patchGroup
        inline bool valid() const;

        //- Find other patch in same region. Returns index of patch or -1.
        label findOtherPatchID(const polyPatch&) const;

        //- Find other patch and region. Returns index of patch and sets
        //  otherRegion to name of region. Fatal error if patch not found
        label findOtherPatchID(const polyPatch&, word&) const;

        //- Write the data as a dictionary
        void write(Ostream&) const;


    // IOstream Operators

        friend Ostream& operator<<(Ostream&, const coupleGroupIdentifier&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "coupleGroupIdentifierI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
