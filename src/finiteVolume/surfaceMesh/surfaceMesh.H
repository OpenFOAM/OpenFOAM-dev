/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2026 OpenFOAM Foundation
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
    Mesh data for surface fields

\*---------------------------------------------------------------------------*/

#ifndef surfaceMesh_H
#define surfaceMesh_H

#include "fvBoundaryMesh.H"
#include "GeometricFieldSources.H"
#include "fvsPatchFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;
    // class fvBoundaryMesh;
class faceZone;

/*---------------------------------------------------------------------------*\
                         Class surfaceMesh Declaration
\*---------------------------------------------------------------------------*/

class surfaceMesh
{
    // Private Data

        //- Reference to fvMesh
        const fvMesh& mesh_;


protected:

    // Constructors

        //- Construct from fvMesh
        explicit surfaceMesh(const fvMesh& mesh)
        :
            mesh_(mesh)
        {}


public:

    // Destructor
    virtual ~surfaceMesh()
    {}


    // Public typedefs

        //- Mesh type
        typedef fvMesh Mesh;

        //- Boundary mesh type
        typedef fvBoundaryMesh BoundaryMesh;

        //- Zone type
        typedef faceZone Zone;

        //- Patch field type
        template<class Type>
        using PatchField = fvsPatchField<Type>;

        //- Field source type
        template<class Type>
        using FieldSource = NoFieldSource<Type>;


    // Member Functions

        //- Return size (number of internal faces)
        inline label size() const;

        //- Return fvMesh database
        inline const objectRegistry& thisDb() const;

        //- Return time
        inline const Time& time() const;

        //- Return reference to fvBoundaryMesh
        inline const fvBoundaryMesh& boundary() const;


    // Member Operators

        //- Return reference to fvMesh
        inline const fvMesh& operator()() const;

        bool operator!=(const surfaceMesh&) const;

        bool operator==(const surfaceMesh&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
