/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2023 OpenFOAM Foundation
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
    Foam::cellEdgeAddressingData
    Foam::cellEdgeAddressing
    Foam::cellEdgeAddressingList

Description
    Engine for providing cell-local cell-edge to face-edge addressing

SourceFiles
    cellEdgeAddressing.C
    cellEdgeAddressingI.H

\*---------------------------------------------------------------------------*/

#ifndef cellEdgeAddressing_H
#define cellEdgeAddressing_H

#include "DemandDrivenMeshObject.H"
#include "polyMesh.H"
#include "CompactListList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class polyDistributionMap;
class polyTopoChangeMap;
class polyMeshMap;
class cellEdgeAddressingList;

/*---------------------------------------------------------------------------*\
                      Class cellEdgeAddressingData Declaration
\*---------------------------------------------------------------------------*/

class cellEdgeAddressingData
{
private:

    // Private Data

        //- Map from cell-face-index and face-edge-index to cell-edge-index
        CompactListList<label> cfiAndFeiToCei_;

        //- Map from cell-edge-index to cell-face-index and face-edge-index.
        //  Note that there are two sets of indices for each cell-edge,
        //  corresponding to the two faces that are connected to that edge.
        List<Pair<labelPair>> ceiToCfiAndFei_;

        //- For each cell-face, whether or not the cell owns it. The same as
        //  the concept of owner in the polyMesh, except this is calculated
        //  using the direction of the face-edges. I.e., if the first face is
        //  know to be "owned" by the face then any face that is edge connected
        //  to the first face is also owned by the cell if the connecting
        //  face-edge is numbered in an opposing direction.
        boolList cOwns_;


public:

    // Constructors

        //- Construct from a cell and its faces
        cellEdgeAddressingData
        (
            const cell& c,
            const faceList& fs,
            const bool cOwnsFirst
        );


    // Member Functions

        //- Map from cell-face-index and face-edge-index to cell-edge-index
        inline const CompactListList<label>& cfiAndFeiToCei() const;

        //- Map from cell-edge-index to cell-face-index and face-edge-index
        inline const List<Pair<labelPair>>& ceiToCfiAndFei() const;

        //- For each cell-face, whether or not the cell owns it
        inline const boolList& cOwns() const;
};


/*---------------------------------------------------------------------------*\
                    Class cellEdgeAddressing Declaration
\*---------------------------------------------------------------------------*/

class cellEdgeAddressing
{
private:

    // Private Data

        //- The cell edge addressing list
        const cellEdgeAddressingList& cAddrs_;

        //- The cell index
        const label celli_;


public:

    // Constructors

        //- Construct from the cell edge addressing list and a cell index
        inline cellEdgeAddressing
        (
            const cellEdgeAddressingList& cAddrs,
            const label celli
        );


    // Member Functions

        //- Map from cell-face-index and face-edge-index to cell-edge-index
        inline const CompactListList<label>& cfiAndFeiToCei() const;

        //- Map from cell-edge-index to cell-face-index and face-edge-index
        inline const List<Pair<labelPair>>& ceiToCfiAndFei() const;

        //- For each cell-face, whether or not the cell owns it
        inline const boolList& cOwns() const;
};


/*---------------------------------------------------------------------------*\
                    Class cellEdgeAddressingList Declaration
\*---------------------------------------------------------------------------*/

class cellEdgeAddressingList
:
    public DemandDrivenMeshObject
    <
        polyMesh,
        UpdateableMeshObject,
        cellEdgeAddressingList
    >
{
    // Private Data

        //- Edge addressing for each cell
        mutable PtrList<cellEdgeAddressingData> list_;


protected:

    friend class DemandDrivenMeshObject
    <
        polyMesh,
        UpdateableMeshObject,
        cellEdgeAddressingList
    >;

    // Protected Constructors

        //- Construct for a mesh
        cellEdgeAddressingList(const polyMesh& mesh);


public:

    // Member Functions

        //- Update following mesh motion
        virtual bool movePoints();

        //- Update following mesh distribution
        virtual void distribute(const polyDistributionMap& map);

        //- Update following topology change
        virtual void topoChange(const polyTopoChangeMap& map);

        //- Update following mapping
        virtual void mapMesh(const polyMeshMap& map);


    // Member Operators

        //- Get the addressing data for a given cell
        const cellEdgeAddressingData& data(const label celli) const;

        //- Get the addressing for a given cell
        inline const cellEdgeAddressing operator[](const label celli) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "cellEdgeAddressingI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
