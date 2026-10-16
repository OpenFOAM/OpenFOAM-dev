/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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
    Foam::fvCellZoneIndexBoundaryField

Description
    Mesh object which provides a boundary label field containing the indices of
    the zones in the adjacent cells. Relates to a supplied sub-set of the
    zones. The resulting indices are into that subset rather than into the list
    of all the cell zones in the mesh. If the zones overlap, the latter zone
    takes precedence.

    Useful when applying a zonal model to both the internal and the boundary of
    a geometric field. Handles couplings so that processor faces are given the
    zone of their neighbouring cell. This ensures that a zonal model is applied
    correctly on processor patches.

SourceFiles
    zonalBoundaryIndices.C

\*---------------------------------------------------------------------------*/

#ifndef fvCellZoneIndexBoundaryField_H
#define fvCellZoneIndexBoundaryField_H

#include "DemandDrivenMeshObject.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class fvCellZoneIndexBoundaryField Declaration
\*---------------------------------------------------------------------------*/

class fvCellZoneIndexBoundaryField
:
    public DemandDrivenMeshObject
    <
        fvMesh,
        MoveableMeshObject,
        fvCellZoneIndexBoundaryField
    >,
    public labelListList
{
protected:

    friend class DemandDrivenMeshObject
    <
        fvMesh,
        MoveableMeshObject,
        fvCellZoneIndexBoundaryField
    >;


    // Protected Constructors

        //- Construct from a mesh and a list of zones
        explicit fvCellZoneIndexBoundaryField
        (
            const word& name,
            const fvMesh& mesh,
            const labelList& zoneIndices
        );


public:

    //- Runtime type information
    TypeName("fvCellZoneIndexBoundaryField");


    // Constructors

        //- Disallow default bitwise copy construction
        fvCellZoneIndexBoundaryField
        (
            const fvCellZoneIndexBoundaryField&
        ) = delete;

        //- Construct and return for the given set of zones
        static const fvCellZoneIndexBoundaryField& New
        (
            const fvMesh& mesh,
            const labelList& zoneIndices
        );


    //- Destructor
    virtual ~fvCellZoneIndexBoundaryField();


    // Member Functions

        //- Update for mesh motion
        virtual bool movePoints();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
