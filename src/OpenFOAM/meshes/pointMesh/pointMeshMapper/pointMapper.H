/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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
    Foam::pointMapper

Description
    This object provides mapping and fill-in information for point data
    between the two meshes after the topological change.  It is
    constructed from polyTopoChangeMap.

SourceFiles
    pointMapper.C

\*---------------------------------------------------------------------------*/

#ifndef pointMapper_H
#define pointMapper_H

#include "morphFieldMapper.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class pointMesh;
class polyTopoChangeMap;
class polyMesh;

/*---------------------------------------------------------------------------*\
                         Class pointMapper Declaration
\*---------------------------------------------------------------------------*/

class pointMapper
:
    public morphFieldMapper
{
    // Private Data

        //- Reference to pointMesh
        const pointMesh& pMesh_;

        //- Reference to polyTopoChangeMap
        const polyTopoChangeMap& map_;

        //- Are there any inserted (unmapped) points
        bool insertedPoints_;

        //- Is the mapping direct
        bool direct_;


    // Demand-driven private data

        //- Direct addressing (only one for of addressing is used)
        mutable labelList* directAddrPtr_;

        //- Interpolated addressing (only one for of addressing is used)
        mutable labelListList* interpolationAddrPtr_;

        //- Interpolation weights
        mutable scalarListList* weightsPtr_;

        //- Inserted points
        mutable labelList* insertedPointLabelsPtr_;


    // Private Member Functions

        //- Calculate addressing for mapping with inserted points
        void calcAddressing() const;

        //- Clear out local storage
        void clearOut();


public:

    // Constructors

        //- Construct from polyTopoChangeMap
        pointMapper(const pointMesh&, const polyTopoChangeMap& map);

        //- Disallow default bitwise copy construction
        pointMapper(const pointMapper&) = delete;


    //- Destructor
    virtual ~pointMapper();


    // Member Functions

        //- Is the mapping direct
        virtual bool direct() const
        {
            return direct_;
        }

        //- Return direct addressing
        virtual const labelUList& directAddressing() const;

        //- Return interpolated addressing
        virtual const labelListList& addressing() const;

        //- Return interpolation weights
        virtual const scalarListList& weights() const;

        //- Return size before mapping
        virtual label sizeBeforeMapping() const;

        //- Are there any inserted points
        virtual bool insertedObjects() const
        {
            return insertedPoints_;
        }

        //- Return list of inserted points
        virtual const labelList& insertedObjectLabels() const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const pointMapper&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
