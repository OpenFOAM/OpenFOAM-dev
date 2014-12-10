/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
    Foam::meshSearchFACECENTRETETSMeshObject

Description
    MeshObject wrapper around meshSearch(mesh,  polyMesh::FACECENTRETETS).

SourceFiles

\*---------------------------------------------------------------------------*/

#ifndef meshSearchFACECENTRETETSMeshObject_H
#define meshSearchFACECENTRETETSMeshObject_H

#include "MeshObject.H"
#include "meshSearch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class meshSearchFACECENTRETETSMeshObject Declaration
\*---------------------------------------------------------------------------*/

class meshSearchFACECENTRETETSMeshObject
:
    public MeshObject
    <
        polyMesh,
        GeometricMeshObject,
        meshSearchFACECENTRETETSMeshObject
    >,
    public meshSearch
{

public:

    // Declare name of the class and its debug switch
    TypeName("meshSearchFACECENTRETETSMeshObject");


    // Constructors

        //- Constructor given polyMesh
        explicit meshSearchFACECENTRETETSMeshObject(const polyMesh& mesh);

    //- Destructor
    virtual ~meshSearchFACECENTRETETSMeshObject()
    {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
