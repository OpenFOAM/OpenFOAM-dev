/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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
    Foam::polyTopoChanger

Description
    List of mesh modifiers

SourceFiles
    polyTopoChanger.C

\*---------------------------------------------------------------------------*/

#ifndef polyTopoChanger_H
#define polyTopoChanger_H

#include "polyMeshModifier.H"
#include "PtrList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class polyMesh;
class polyTopoChangeMap;
class polyBoundaryMesh;


/*---------------------------------------------------------------------------*\
                       Class polyTopoChanger Declaration
\*---------------------------------------------------------------------------*/

class polyTopoChanger
:
    public PtrList<polyMeshModifier>
{
    // Private member functions

        //- Is topology change required
        bool changeTopology() const;

        //- Return topology change request
        autoPtr<polyTopoChange> topoChangeRequest() const;

        //- Force recalculation of locally stored data on topological change
        void update(const polyTopoChangeMap& m);


protected:

    // Protected member data

        //- Reference to mesh
        polyMesh& mesh_;


public:

    //- Runtime type information
    TypeName("polyTopoChanger");


    // Constructors

        //- Read constructor for given polyMesh
        explicit polyTopoChanger(polyMesh&);

        //- Disallow default bitwise copy construction
        polyTopoChanger(const polyTopoChanger&) = delete;


    //- Destructor
    virtual ~polyTopoChanger()
    {}


    // Member Functions

        //- Return the mesh reference
        const polyMesh& mesh() const
        {
            return mesh_;
        }

        autoPtr<polyTopoChangeMap> changeMesh
        (
            const bool inflate,
            const bool syncParallel = true,
            const bool orderCells = false,
            const bool orderPoints = false
        );


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const polyTopoChanger&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
