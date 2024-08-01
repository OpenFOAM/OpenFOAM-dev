/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2024 OpenFOAM Foundation
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
    fvMeshTopoChanger which maps the fields to a new mesh or sequence of meshes

    which can optionally be mapped to repeatedly for example in multi-cycle
    engine cases or cycled through for symmetric forward and reverse motion.

Usage
    \table
        Property  | Description                   | Required | Default value
        libs      | Libraries to load             | no       |
        times     | List of times for the meshes  | yes      |
        repeat    | Repetition period             | no       |
        cycle     | Cycle period                  | no       |
        begin     | Begin time for the meshes     | no       | Time::beginTime()
        timeDelta | Time tolerance used for time -> index | yes      |
    \endtable

    Examples of the mesh-to-mesh mapping for the multi-cycle
    tutorials/incompressibleFluid/movingCone case:
    \verbatim
    topoChanger
    {
        type    meshToMesh;

        libs    ("libmeshToMeshTopoChanger.so");

        times   (0.0015 0.003);

        cycle   #calc "1.0/300.0";
        begin   0;

        timeDelta 1e-6;
    }
    \endverbatim

SourceFiles
    meshToMesh_fvMeshTopoChanger.C

\*---------------------------------------------------------------------------*/

#ifndef meshToMesh_fvMeshTopoChanger_H
#define meshToMesh_fvMeshTopoChanger_H

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

        //- Mesh-to-mesh control dictionary
        dictionary dict_;

        //- List of mesh mapping times
        scalarList times_;

        //- Time delta used for time -> index
        scalar timeDelta_;

        //- Hash set of mesh mapping time indices
        HashSet<int64_t, Hash<int64_t>> timeIndices_;

        //- Optional begin time for the meshes
        scalar begin_;

        //- Optional repetition period
        scalar repeat_;

        //- Optional cycle period
        scalar cycle_;

        //- The time index used for updating
        label timeIndex_;


    // Private Member Functions

        //- Return true if the set of meshes are being traversed in the forward
        //  sequence or false if cycling is currently traversing the meshes in
        //  reverse order
        bool forward() const;

        //- Return the user time mapped to the mesh sequence
        //  handling the repeat or cycle option
        scalar meshTime() const;

        //- Interpolate U's to Uf's
        void interpolateUfs();


public:

    //- Runtime type information
    TypeName("meshToMesh");


    // Constructors

        //- Construct from fvMesh and dictionary
        meshToMesh(fvMesh&, const dictionary& dict);

        //- Disallow default bitwise copy construction
        meshToMesh(const meshToMesh&) = delete;


    //- Destructor
    virtual ~meshToMesh();


    // Member Functions

        scalar timeToNextMesh() const;

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
