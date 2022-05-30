/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022 OpenFOAM Foundation
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
    Foam::fvMeshTopoChangers::meshToMesh

Description
    fvMeshTopoChanger which maps to new mesh

SourceFiles
    fvMeshTopoChangersMeshToMesh.C

\*---------------------------------------------------------------------------*/

#ifndef fvMeshTopoChangersMeshToMesh_H
#define fvMeshTopoChangersMeshToMesh_H

#include "fvMeshTopoChanger.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fvMeshTopoChangers
{

/*---------------------------------------------------------------------------*\
                Class fvMeshTopoChangers::meshToMesh Declaration
\*---------------------------------------------------------------------------*/

class meshToMesh
:
    public fvMeshTopoChanger
{
    // Private Data

        //- Optional names of the patches which "cut" the internal mesh
        //  and hence cannot be mapped from existing patches
        wordList cuttingPatches_;

        //- List of mesh mapping times
        scalarList times_;

        //- Time delta used for time -> index
        scalar timeDelta_;

        //- Hash set of mesh mapping time indices
        labelHashSet timeIndices_;

        //- The time index used for updating
        label timeIndex_;


    // Private Member Functions

        //- Read the parameters from dictionary
        void readDict();

        //- Find the U field name corresponding to Uf
        word Uname(const surfaceVectorField& Uf) const;

        //- Interpolate U's to Uf's
        void interpolateUfs();


public:

    //- Runtime type information
    TypeName("meshToMesh");


    // Constructors

        //- Construct from fvMesh
        explicit meshToMesh(fvMesh& io);

        //- Disallow default bitwise copy construction
        meshToMesh(const meshToMesh&) = delete;


    //- Destructor
    virtual ~meshToMesh();


    // Member Functions

        const scalarList& times() const
        {
            return times_;
        }

        scalar timeDelta() const
        {
            return timeDelta_;
        }

        //- Update the mesh for both mesh motion and topology change
        virtual bool update();

        //- Update corresponding to the given map
        virtual void topoChange(const polyTopoChangeMap&);

        //- Update from another mesh using the given map
        virtual void mapMesh(const polyMeshMap&);

        //- Update corresponding to the given distribution map
        virtual void distribute(const polyDistributionMap&);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const meshToMesh&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvMeshTopoChangers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
