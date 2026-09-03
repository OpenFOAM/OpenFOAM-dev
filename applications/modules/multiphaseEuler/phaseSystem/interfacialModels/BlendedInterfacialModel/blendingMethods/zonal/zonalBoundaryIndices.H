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
    Foam::blendingMethods::zonal::boundaryIndices

Description
    Mesh object for storing zone boundary indexing information

SourceFiles
    zonalBoundaryIndices.C

\*---------------------------------------------------------------------------*/

#ifndef zonalBoundaryIndices_H
#define zonalBoundaryIndices_H

#include "zonal.H"
#include "DemandDrivenMeshObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace blendingMethods
{

/*---------------------------------------------------------------------------*\
                    Class zonal::boundaryIndices Declaration
\*---------------------------------------------------------------------------*/

class zonal::boundaryIndices
:
    public DemandDrivenMeshObject<fvMesh, MoveableMeshObject, boundaryIndices>,
    public List<List<label>>
{
    // Private Data

        //- The indices of the zones
        const labelList zoneIndices_;


protected:

    friend class DemandDrivenMeshObject
    <
        fvMesh,
        MoveableMeshObject,
        boundaryIndices
    >;


    // Protected Constructors

        //- Construct from a mesh and a list of zones
        explicit boundaryIndices
        (
            const word& name,
            const fvMesh& mesh,
            const labelList& zoneIndices
        );


public:

    //- Runtime type information
    TypeName("zonal::boundaryIndices");


    // Constructors

        //- Disallow default bitwise copy construction
        boundaryIndices(const boundaryIndices&) = delete;

        //- Construct and return for the given set of zones
        static const boundaryIndices& New
        (
            const fvMesh& mesh,
            const labelList& zoneIndices
        );


    //- Destructor
    virtual ~boundaryIndices();


    // Member Functions

        //- Update for mesh motion
        virtual bool movePoints();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace blendingMethods
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
