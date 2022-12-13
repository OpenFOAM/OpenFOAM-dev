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
    Foam::pureUpwindCFCCellToFaceStencilObject

Description

SourceFiles

\*---------------------------------------------------------------------------*/

#ifndef pureUpwindCFCCellToFaceStencilObject_H
#define pureUpwindCFCCellToFaceStencilObject_H

#include "extendedUpwindCellToFaceStencil.H"
#include "CFCCellToFaceStencil.H"
#include "DemandDrivenMeshObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class pureUpwindCFCCellToFaceStencilObject Declaration
\*---------------------------------------------------------------------------*/

class pureUpwindCFCCellToFaceStencilObject
:
    public DemandDrivenMeshObject
    <
        fvMesh,
        TopologicalMeshObject,
        pureUpwindCFCCellToFaceStencilObject
    >,
    public extendedUpwindCellToFaceStencil
{

protected:

    friend class DemandDrivenMeshObject
    <
        fvMesh,
        TopologicalMeshObject,
        pureUpwindCFCCellToFaceStencilObject
    >;

    // Protected Constructors

        //- Construct from uncompacted face stencil
        explicit pureUpwindCFCCellToFaceStencilObject
        (
            const fvMesh& mesh
        )
        :
            DemandDrivenMeshObject
            <
                fvMesh,
                TopologicalMeshObject,
                pureUpwindCFCCellToFaceStencilObject
            >(mesh),
            extendedUpwindCellToFaceStencil(CFCCellToFaceStencil(mesh))
        {
            if (extendedCellToFaceStencil::debug)
            {
                Info<< "Generated pure upwind stencil " << type()
                    << nl << endl;
                writeStencilStats(Info, ownStencil(), ownMap());
            }
        }


public:

    TypeName("pureUpwindCFCCellToFaceStencil");

    //- Destructor
    virtual ~pureUpwindCFCCellToFaceStencilObject()
    {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
