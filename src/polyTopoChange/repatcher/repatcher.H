/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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
    Foam::repatcher

Description
    A mesh which allows changes in the patch distribution of the
    boundary faces.  The change in patching is set using changePatchID. For a
    boundary face, a new patch ID is given.

SourceFiles
    repatcher.C

\*---------------------------------------------------------------------------*/

#ifndef repatcher_H
#define repatcher_H

#include "polyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class polyTopoChange;

/*---------------------------------------------------------------------------*\
                   Class repatcher Declaration
\*---------------------------------------------------------------------------*/

class repatcher
{
    // Private Data

        //- The polyMesh to be repatched
        polyMesh& mesh_;

        //- Topological change to accumulated all mesh changes
        autoPtr<polyTopoChange> meshModPtr_;


    // Private Member Functions

        //- Demand driven access to polyTopoChange
        polyTopoChange& meshMod();


public:

    // Constructors

        //- Construct for given mesh
        explicit repatcher(polyMesh& mesh);

        //- Disallow default bitwise copy construction
        repatcher(const repatcher&) = delete;


    //- Destructor
    virtual ~repatcher();


    // Member Functions

        //- Change patches.
        void changePatches(const List<polyPatch*>& patches);

        //- Change patch ID for a boundary face. Note: patchID should be in new
        //  numbering.
        void changePatchID
        (
            const label faceID,
            const label patchID
        );

        //- Change anchor point (zero'th point of face) for a boundary face.
        void changeAnchorPoint
        (
            const label faceID,
            const label fp
        );

        //- Re-patch the mesh
        void repatch();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const repatcher&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
