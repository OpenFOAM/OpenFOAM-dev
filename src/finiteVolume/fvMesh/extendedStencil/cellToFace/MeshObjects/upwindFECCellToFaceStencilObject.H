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
    Foam::upwindFECCellToFaceStencilObject

Description

SourceFiles

\*---------------------------------------------------------------------------*/

#ifndef upwindFECCellToFaceStencilObject_H
#define upwindFECCellToFaceStencilObject_H

#include "extendedUpwindCellToFaceStencil.H"
#include "FECCellToFaceStencil.H"
#include "DemandDrivenMeshObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class upwindFECCellToFaceStencilObject Declaration
\*---------------------------------------------------------------------------*/

class upwindFECCellToFaceStencilObject
:
    public DemandDrivenMeshObject
    <
        fvMesh,
        TopologicalMeshObject,
        upwindFECCellToFaceStencilObject
    >,
    public extendedUpwindCellToFaceStencil
{

protected:

    friend class DemandDrivenMeshObject
    <
        fvMesh,
        TopologicalMeshObject,
        upwindFECCellToFaceStencilObject
    >;

    // Protected Constructors

        //- Construct from uncompacted face stencil
        upwindFECCellToFaceStencilObject
        (
            const fvMesh& mesh,
            const bool pureUpwind,
            const scalar minOpposedness
        )
        :
            DemandDrivenMeshObject
            <
                fvMesh,
                TopologicalMeshObject,
                upwindFECCellToFaceStencilObject
            >(mesh),
            extendedUpwindCellToFaceStencil
            (
                FECCellToFaceStencil(mesh),
                pureUpwind,
                minOpposedness
            )
        {
            if (extendedCellToFaceStencil::debug)
            {
                Info<< "Generated off-centred stencil " << type()
                    << nl << endl;
                writeStencilStats(Info, ownStencil(), ownMap());
            }
        }

public:

    TypeName("upwindFECCellToFaceStencil");

    //- Destructor
    virtual ~upwindFECCellToFaceStencilObject()
    {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
