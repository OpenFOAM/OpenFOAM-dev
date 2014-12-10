/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    Foam::fvBoundaryMeshMapper

Description
    Foam::fvBoundaryMeshMapper

\*---------------------------------------------------------------------------*/

#ifndef fvBoundaryMeshMapper_H
#define fvBoundaryMeshMapper_H

#include "PtrList.H"
#include "fvPatchMapper.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class fvBoundaryMeshMapper Declaration
\*---------------------------------------------------------------------------*/

class fvBoundaryMeshMapper
:
    public PtrList<fvPatchMapper>
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        fvBoundaryMeshMapper(const fvBoundaryMeshMapper&);

        //- Disallow default bitwise assignment
        void operator=(const fvBoundaryMeshMapper&);


public:

    // Constructors

        //- Construct from components
        fvBoundaryMeshMapper
        (
            const fvMesh& mesh,
            const faceMapper& faceMap
        )
        :
            PtrList<fvPatchMapper>(mesh.boundary().size())
        {
            const fvBoundaryMesh& patches = mesh.boundary();

            forAll(patches, patchI)
            {
                set
                (
                    patchI,
                    new fvPatchMapper
                    (
                        patches[patchI],
                        faceMap
                    )
                );
            }
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
