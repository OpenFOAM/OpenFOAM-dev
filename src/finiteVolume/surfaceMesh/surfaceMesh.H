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
    Foam::surfaceMesh

Description
    Mesh data needed to do the Finite Volume discretisation.

\*---------------------------------------------------------------------------*/

#ifndef surfaceMesh_H
#define surfaceMesh_H

#include "GeoMesh.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class surfaceMesh Declaration
\*---------------------------------------------------------------------------*/

class surfaceMesh
:
    public GeoMesh<fvMesh>
{
public:

    // Public typedefs

        //- Patch field type
        template<class Type>
        using PatchField = fvsPatchField<Type>;

        //- Field source type
        template<class Type>
        using FieldSource = NoFieldSource<Type>;


    // Constructors

        //- Construct from fvMesh
        explicit surfaceMesh(const fvMesh& mesh)
        :
            GeoMesh<fvMesh>(mesh)
        {}


    // Member Functions

        //- Return size
        label size() const
        {
            return size(mesh_);
        }

        //- Return size
        static label size(const Mesh& mesh)
        {
            return mesh.nInternalFaces();
        }

        //- Return face centres
        const surfaceVectorField& C()
        {
            return mesh_.Cf();
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
