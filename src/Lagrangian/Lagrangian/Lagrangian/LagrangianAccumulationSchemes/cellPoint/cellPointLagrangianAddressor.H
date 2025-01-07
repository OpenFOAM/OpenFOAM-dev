/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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
    Foam::cellPointLagrangianAddressor

Description
    Helper class for the cell-point Lagrangian accumulation scheme. This is a
    mesh object in order to provide a place to store re-usable accumulation
    weights and workspace.

SourceFiles
    cellPointLagrangianAddressor.C

\*---------------------------------------------------------------------------*/

#ifndef cellPointLagrangianAddressor_H
#define cellPointLagrangianAddressor_H

#include "CompactListList.H"
#include "DemandDrivenMeshObject.H"
#include "polyMesh.H"
#include "remote.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class cellPointLagrangianAddressor Declaration
\*---------------------------------------------------------------------------*/

class cellPointLagrangianAddressor
:
    public DemandDrivenMeshObject
    <
        polyMesh,
        MoveableMeshObject,
        cellPointLagrangianAddressor
    >
{
    // Private Data

        //- Map from mesh point to coupled point index
        labelList pointCoupledPoint_;

        //- Map from coupled point to mesh point index
        labelList coupledPointPoint_;

        //- For every coupled point, the list of points it is connected to
        CompactListList<remote> coupledPointCoupledPoint_;


protected:

    // Protected Constructors

        //- Allow the base mesh object to use the protected constructor
        friend class DemandDrivenMeshObject
        <
            polyMesh,
            MoveableMeshObject,
            cellPointLagrangianAddressor
        >;

        //- Construct from mesh
        explicit cellPointLagrangianAddressor(const polyMesh& mesh);


public:

    // Static Data Members

        //- Runtime type information
        TypeName("cellPointLagrangianAddressor");


    // Constructors

        //- Disallow default bitwise copy construction
        cellPointLagrangianAddressor(const cellPointLagrangianAddressor&) =
            delete;


    //- Destructor
    ~cellPointLagrangianAddressor();


    // Member Functions

        //- Update for mesh motion
        virtual bool movePoints();

        //- Synchronise the identified set of points
        void sync(labelList& pointIndex, DynamicList<label>& indexPoint) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const cellPointLagrangianAddressor&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
