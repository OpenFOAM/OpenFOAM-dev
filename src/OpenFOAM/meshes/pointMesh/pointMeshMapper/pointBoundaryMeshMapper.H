/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
    Foam::pointBoundaryMeshMapper

Description
    Foam::pointBoundaryMeshMapper

\*---------------------------------------------------------------------------*/

#ifndef pointBoundaryMeshMapper_H
#define pointBoundaryMeshMapper_H

#include "PtrList.H"
#include "pointPatchMapper.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class pointBoundaryMeshMapper Declaration
\*---------------------------------------------------------------------------*/

class pointBoundaryMeshMapper
:
    public PtrList<pointPatchMapper>
{

public:

    // Constructors

        //- Construct from components
        pointBoundaryMeshMapper
        (
            const pointMesh& mesh,
            const pointMapper& pointMap,
            const mapPolyMesh& mpm
        )
        :
            PtrList<pointPatchMapper>(mesh.boundary().size())
        {
            const pointBoundaryMesh& patches = mesh.boundary();

            forAll(patches, patchi)
            {
                set
                (
                    patchi,
                    new pointPatchMapper
                    (
                        patches[patchi],
                        pointMap,
                        mpm
                    )
                );
            }
        }

        //- Disallow default bitwise copy construction
        pointBoundaryMeshMapper(const pointBoundaryMeshMapper&) = delete;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const pointBoundaryMeshMapper&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
