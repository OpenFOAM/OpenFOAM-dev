/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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
    Foam::centredCFCFaceToCellStencilObject

Description

SourceFiles

\*---------------------------------------------------------------------------*/

#ifndef centredCFCFaceToCellStencilObject_H
#define centredCFCFaceToCellStencilObject_H

#include "extendedCentredFaceToCellStencil.H"
#include "CFCFaceToCellStencil.H"
#include "DemandDrivenMeshObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class centredCFCFaceToCellStencilObject Declaration
\*---------------------------------------------------------------------------*/

class centredCFCFaceToCellStencilObject
:
    public DemandDrivenMeshObject
    <
        fvMesh,
        TopologicalMeshObject,
        centredCFCFaceToCellStencilObject
    >,
    public extendedCentredFaceToCellStencil
{

protected:

    friend class DemandDrivenMeshObject
    <
        fvMesh,
        TopologicalMeshObject,
        centredCFCFaceToCellStencilObject
    >;

    // Protected Constructors

        //- Construct from uncompacted face stencil
        explicit centredCFCFaceToCellStencilObject
        (
            const fvMesh& mesh
        )
        :
            DemandDrivenMeshObject
            <
                fvMesh,
                TopologicalMeshObject,
                centredCFCFaceToCellStencilObject
            >(mesh),
            extendedCentredFaceToCellStencil(CFCFaceToCellStencil(mesh))
        {}


public:

    TypeName("centredCFCFaceToCellStencil");

    //- Destructor
    virtual ~centredCFCFaceToCellStencilObject()
    {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
