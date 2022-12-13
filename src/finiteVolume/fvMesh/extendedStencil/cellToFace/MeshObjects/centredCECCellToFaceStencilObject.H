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
    Foam::centredCECCellToFaceStencilObject

Description

SourceFiles

\*---------------------------------------------------------------------------*/

#ifndef centredCECCellToFaceStencilObject_H
#define centredCECCellToFaceStencilObject_H

#include "extendedCentredCellToFaceStencil.H"
#include "CECCellToFaceStencil.H"
#include "DemandDrivenMeshObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class centredCECCellToFaceStencilObject Declaration
\*---------------------------------------------------------------------------*/

class centredCECCellToFaceStencilObject
:
    public DemandDrivenMeshObject
    <
        fvMesh,
        TopologicalMeshObject,
        centredCECCellToFaceStencilObject
    >,
    public extendedCentredCellToFaceStencil
{

protected:

    friend class DemandDrivenMeshObject
    <
        fvMesh,
        TopologicalMeshObject,
        centredCECCellToFaceStencilObject
    >;

    // Protected Constructors

        //- Construct from uncompacted face stencil
        explicit centredCECCellToFaceStencilObject
        (
            const fvMesh& mesh
        )
        :
            DemandDrivenMeshObject
            <
                fvMesh,
                TopologicalMeshObject,
                centredCECCellToFaceStencilObject
            >(mesh),
            extendedCentredCellToFaceStencil(CECCellToFaceStencil(mesh))
        {
            if (extendedCellToFaceStencil::debug)
            {
                Info<< "Generated centred stencil " << type()
                    << nl << endl;
                writeStencilStats(Info, stencil(), map());
            }
        }


public:

    TypeName("centredCECCellToFaceStencil");

    //- Destructor
    virtual ~centredCECCellToFaceStencilObject()
    {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
