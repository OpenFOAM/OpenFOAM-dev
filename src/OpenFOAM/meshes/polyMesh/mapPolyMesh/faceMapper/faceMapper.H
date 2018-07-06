/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    Foam::faceMapper

Description
    This object provides mapping and fill-in information for face data
    between the two meshes after the topological change.  It is
    constructed from mapPolyMesh.

SourceFiles
    faceMapper.C

\*---------------------------------------------------------------------------*/

#ifndef faceMapper_H
#define faceMapper_H

#include "morphFieldMapper.H"
#include "HashSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class polyMesh;
class mapPolyMesh;

/*---------------------------------------------------------------------------*\
                           Class faceMapper Declaration
\*---------------------------------------------------------------------------*/

class faceMapper
:
    public morphFieldMapper
{
    // Private data

        //- Reference to polyMesh
        const polyMesh& mesh_;

        //- Reference to mapPolyMesh
        const mapPolyMesh& mpm_;

        //- Are there any inserted (unmapped) faces
        bool insertedFaces_;

        //- Is the mapping direct
        bool direct_;


    // Demand-driven private data

        //- Direct addressing (only one for of addressing is used)
        mutable labelList* directAddrPtr_;

        //- Interpolated addressing (only one for of addressing is used)
        mutable labelListList* interpolationAddrPtr_;

        //- Interpolation weights
        mutable scalarListList* weightsPtr_;

        //- Inserted faces
        mutable labelList* insertedFaceLabelsPtr_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        faceMapper(const faceMapper&);

        //- Disallow default bitwise assignment
        void operator=(const faceMapper&);


        //- Calculate addressing for mapping with inserted faces
        void calcAddressing() const;

        //- Clear out local storage
        void clearOut();


public:

    // Static data members

    // Constructors

        //- Construct from mapPolyMesh
        faceMapper(const mapPolyMesh& mpm);


    //- Destructor
    virtual ~faceMapper();


    // Member Functions

        //- Return size
        virtual label size() const;

        //- Return size of field before mapping
        virtual label sizeBeforeMapping() const;

        //- Return number of internal faces before mapping
        virtual label internalSizeBeforeMapping() const;

        //- Is the mapping direct
        virtual bool direct() const
        {
            return direct_;
        }

        virtual bool hasUnmapped() const
        {
            return insertedObjects();
        }

        //- Return direct addressing
        virtual const labelUList& directAddressing() const;

        //- Return interpolated addressing
        virtual const labelListList& addressing() const;

        //- Return interpolaion weights
        virtual const scalarListList& weights() const;

        //- Return flux flip map
        virtual const labelHashSet& flipFaceFlux() const;

        //- Return number of old internalFaces
        virtual label nOldInternalFaces() const;

        //- Return old patch starts
        virtual const labelList& oldPatchStarts() const;

        //- Return old patch sizes
        virtual const labelList& oldPatchSizes() const;

        //- Are there any inserted faces
        virtual bool insertedObjects() const
        {
            return insertedFaces_;
        }

        //- Return list of inserted faces
        virtual const labelList& insertedObjectLabels() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
