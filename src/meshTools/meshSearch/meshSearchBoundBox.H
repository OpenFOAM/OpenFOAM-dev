/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2025 OpenFOAM Foundation
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
    Foam::meshSearchBoundBox

Description
    Mesh object to hold a bound box for use by mesh search engines

SourceFiles
    meshSearchBoundBox.C

\*---------------------------------------------------------------------------*/

#ifndef meshSearchBoundBox_H
#define meshSearchBoundBox_H

#include "DemandDrivenMeshObject.H"
#include "polyMesh.H"
#include "treeBoundBox.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class treeBoundBox;

/*---------------------------------------------------------------------------*\
                         Class meshSearchBoundBox Declaration
\*---------------------------------------------------------------------------*/

class meshSearchBoundBox
:
    public DemandDrivenMeshObject
    <
        polyMesh,
        DeletableMeshObject,
        meshSearchBoundBox
    >
{
    // Private Data

        //- Bound box
        const treeBoundBox bb_;


public:

    //- Run-time type information
    TypeName("meshSearchBoundBox");


    // Constructors

        //- Construct from mesh
        meshSearchBoundBox(const polyMesh& mesh);

        //- Disallow default bitwise copy construction
        meshSearchBoundBox(const meshSearchBoundBox&) = delete;


    //- Destructor
    virtual ~meshSearchBoundBox();


    // Member Functions

        // Access

            //- Access the bound box
            inline const treeBoundBox& bb() const
            {
                return bb_;
            }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const meshSearchBoundBox&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
