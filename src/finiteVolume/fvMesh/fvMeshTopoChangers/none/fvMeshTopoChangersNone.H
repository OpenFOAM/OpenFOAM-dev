/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021 OpenFOAM Foundation
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
    Foam::fvMeshTopoChangers::none

Description
    Dummy fvMeshTopoChanger which does not move the mesh points

SourceFiles
    none.C

\*---------------------------------------------------------------------------*/

#ifndef fvMeshTopoChangersNone_H
#define fvMeshTopoChangersNone_H

#include "fvMeshTopoChanger.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fvMeshTopoChangers
{

/*---------------------------------------------------------------------------*\
                        Class none Declaration
\*---------------------------------------------------------------------------*/

class none
:
    public fvMeshTopoChanger
{
public:

    //- Runtime type information
    TypeName("none");


    // Constructors

        //- Construct from fvMesh
        none(fvMesh& mesh);

        //- Disallow default bitwise copy construction
        none(const none&) = delete;


    //- Destructor
    ~none();


    // Member Functions

        //- The mesh is static and does not change so return false
        virtual bool dynamic() const
        {
            return false;
        }

        //- Dummy update function which does not change the mesh
        virtual bool update();

        //- Update corresponding to the given map
        virtual void updateMesh(const mapPolyMesh&)
        {}

        //- Update corresponding to the given distribution map
        virtual void distribute(const mapDistributePolyMesh&)
        {}


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const none&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvMeshTopoChangers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
