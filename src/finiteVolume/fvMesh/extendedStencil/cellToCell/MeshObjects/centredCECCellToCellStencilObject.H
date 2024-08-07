/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2023 OpenFOAM Foundation
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
    Foam::centredCECCellToCellStencilObject

Description

SourceFiles

\*---------------------------------------------------------------------------*/

#ifndef centredCECCellToCellStencilObject_H
#define centredCECCellToCellStencilObject_H

#include "extendedCentredCellToCellStencil.H"
#include "CECCellToCellStencil.H"
#include "DemandDrivenMeshObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class centredCECCellToCellStencilObject Declaration
\*---------------------------------------------------------------------------*/

class centredCECCellToCellStencilObject
:
    public DemandDrivenMeshObject
    <
        fvMesh,
        MoveableMeshObject,
        centredCECCellToCellStencilObject
    >,
    public extendedCentredCellToCellStencil
{

protected:

    friend class DemandDrivenMeshObject
    <
        fvMesh,
        MoveableMeshObject,
        centredCECCellToCellStencilObject
    >;

    // Protected Constructors

        //- Construct from uncompacted cell stencil
        explicit centredCECCellToCellStencilObject
        (
            const fvMesh& mesh
        )
        :
            DemandDrivenMeshObject
            <
                fvMesh,
                MoveableMeshObject,
                centredCECCellToCellStencilObject
            >(mesh),
            extendedCentredCellToCellStencil(CECCellToCellStencil(mesh))
        {}


public:

    TypeName("centredCECCellToCellStencil");

    //- Destructor
    virtual ~centredCECCellToCellStencilObject()
    {}

    //- Update for mesh motion
    virtual bool movePoints()
    {
        return true;
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
