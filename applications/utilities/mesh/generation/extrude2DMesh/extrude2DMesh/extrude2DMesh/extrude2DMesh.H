/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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
    Foam::extrude2DMesh

Description
    Given a 2D mesh insert all the topology changes to extrude. Does not work
    in parallel.

SourceFiles
    extrude2DMesh.C

\*---------------------------------------------------------------------------*/

#ifndef extrude2DMesh_H
#define extrude2DMesh_H

#include "typeInfo.H"
#include "label.H"
#include "scalar.H"
#include "labelList.H"
#include "dictionary.H"
#include "extrudeModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class polyMesh;
class polyTopoChange;
class mapPolyMesh;
class mapDistributePolyMesh;
class polyBoundaryMesh;

/*---------------------------------------------------------------------------*\
                           Class extrude2DMesh Declaration
\*---------------------------------------------------------------------------*/

class extrude2DMesh
{
    // Private Data

        //- Reference to 2D mesh
        polyMesh& mesh_;

        const dictionary dict_;

        // const dictionary patchDict_;

        const extrudeModel& model_;

        const word modelType_;

        const word patchType_;

        label frontPatchi_;
        label backPatchi_;

    // Private Member Functions

        //- Check the mesh is 2D
        void check2D() const;

        //- Find extrusion direction
        // void findExtrudeDirection();


public:

    //- Runtime type information
    ClassName("extrude2DMesh");


    // Constructors
        extrude2DMesh
        (
            polyMesh&,
            const dictionary& dict,
            const extrudeModel& model
        );

        //- Disallow default bitwise copy construction
        extrude2DMesh(const extrude2DMesh&) = delete;


    //- Destructor
    ~extrude2DMesh();


    // Member Functions

        //- Add front and back patches
        void addFrontBackPatches();

        //- Play commands into polyTopoChange to extrude mesh.
        void setRefinement(polyTopoChange&);

        //- Force recalculation of locally stored data on topological change
        void updateMesh(const mapPolyMesh&)
        {}

        //- Force recalculation of locally stored data for mesh distribution
        void distribute(const mapDistributePolyMesh&)
        {}

        label frontPatchi() const
        {
            return frontPatchi_;
        }

        label backPatchi() const
        {
            return backPatchi_;
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const extrude2DMesh&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
