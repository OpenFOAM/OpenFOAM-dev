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
    Foam::fvMeshMoversLayeredEngine

Description
    Foam::fvMeshMoversLayeredEngine

SourceFiles
    fvMeshMoversLayeredEngine.C

\*---------------------------------------------------------------------------*/

#ifndef fvMeshMoversLayeredEngine_H
#define fvMeshMoversLayeredEngine_H

#include "fvMeshMoversEngine.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fvMeshMovers
{

/*---------------------------------------------------------------------------*\
                Class fvMeshMovers::layeredEngine Declaration
\*---------------------------------------------------------------------------*/

class layeredEngine
:
    public engine
{
    // Private Data

        dimensionedScalar pistonLayers_;


public:

    //- Runtime type information
    TypeName("layeredEngine");


    // Constructors

        //- Construct from fvMesh
        layeredEngine(fvMesh& mesh);

        //- Disallow default bitwise copy construction
        layeredEngine(const layeredEngine&) = delete;


    //- Destructor
    ~layeredEngine();


    // Member Functions

        //- Update the mesh for both mesh motion and topology change
        virtual bool update();

        //- Update corresponding to the given map
        virtual void updateMesh(const mapPolyMesh&);

        //- Update corresponding to the given distribution map
        virtual void distribute(const mapDistributePolyMesh&);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const layeredEngine&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvMeshMovers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
