/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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
    Foam::LagrangianBoundaryMesh

Description
    Boundary part of a Lagrangian mesh. Just a list of Lagrangian patches with
    some added convenience functions.

SourceFiles
    LagrangianBoundaryMesh.C

\*---------------------------------------------------------------------------*/

#ifndef LagrangianBoundaryMesh_H
#define LagrangianBoundaryMesh_H

#include "PtrList.H"
#include "LagrangianPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class LagrangianMesh;
class polyBoundaryMesh;

/*---------------------------------------------------------------------------*\
                         Class LagrangianBoundaryMesh Declaration
\*---------------------------------------------------------------------------*/

class LagrangianBoundaryMesh
:
    public PtrList<LagrangianPatch>
{
    // Private Data

        //- Reference to mesh
        const LagrangianMesh& mesh_;


public:

    // Constructors

        //- Construct from LagrangianMesh and polyBoundaryMesh
        LagrangianBoundaryMesh(const LagrangianMesh&, const polyBoundaryMesh&);

        //- Construct from LagrangianMesh and polyBoundaryMesh and a list of
        //  patch types
        LagrangianBoundaryMesh
        (
            const LagrangianMesh&,
            const polyBoundaryMesh&,
            const wordList&
        );

        //- Disallow default bitwise copy construction
        LagrangianBoundaryMesh(const LagrangianBoundaryMesh&) = delete;


    //- Destructor
    ~LagrangianBoundaryMesh();


    // Member Functions

        //- Return the mesh reference
        const LagrangianMesh& mesh() const
        {
            return mesh_;
        }

        //- Find patch index given a name
        label findIndex(const word& patchName) const;

        //- Find patch indices given a name
        labelList findIndices(const wordRe&, const bool useGroups) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const LagrangianBoundaryMesh&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
