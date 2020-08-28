/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
    Foam::extrudePatchMesh

Description
    Mesh at a patch created on the fly. The following entry should be used
    on the field boundary dictionary:

    Example:
    \verbatim
        // New Shell mesh data

        extrudeModel    linearNormal;
        linearNormalCoeffs
        {
            thickness       40e-6;
        }
        nLayers         50;
        expansionRatio  1;
        columnCells      true;

        // Patch information
        bottomCoeffs
        {
            name        "bottom";
            type        mappedWall;
            sampleMode  nearestPatchFace;
            samplePatch fixedWalls;
            offsetMode  uniform;
            offset      (0 0 0);
        }

        topCoeffs
        {
            name        "top";
            type        patch;
        }

        sideCoeffs
        {
            name        "side";
            type        empty;
        }
    \endverbatim

\*---------------------------------------------------------------------------*/

#ifndef extrudePatchMesh_H
#define extrudePatchMesh_H

#include "extrudeModel.H"
#include "fvMesh.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


/*---------------------------------------------------------------------------*\
                      Class extrudePatchMesh Declaration
\*---------------------------------------------------------------------------*/

class extrudePatchMesh
:
    public fvMesh
{
    // Private Data

        //- Enumeration of patch IDs
        enum patchID
        {
            bottomPatchID,
            topPatchID,
            sidePatchID
        };

        //- Const reference to the patch from which this mesh is extruded
        const polyPatch& extrudedPatch_;

        //- Model dictionary
        dictionary dict_;


    // Private Member Functions

        //- Extrude mesh using polyPatches
        void extrudeMesh(const List<polyPatch*>& regionPatches);


public:

    //- Runtime type information
    TypeName("extrudePatchMesh");


    // Constructors

        //- Construct from mesh, patch and dictionary
        extrudePatchMesh
        (
            const fvMesh&,
            const fvPatch&,
            const dictionary&,
            const word
        );

        //- Construct from mesh, patch, dictionary and new mesh
        //  polyPatch information
        extrudePatchMesh
        (
            const fvMesh&,
            const fvPatch&,
            const dictionary&,
            const word,
            const List<polyPatch*>& polyPatches
        );


    //- Destructor
    virtual ~extrudePatchMesh();


    // Member Functions


        // Access functions

            //- Return region mesh
            const fvMesh& regionMesh() const
            {
                return *this;
            }

            //- Return bottom patch
            const polyPatch& bottomPatch() const
            {
                return this->boundaryMesh()[bottomPatchID];
            }

            //- Return top patch
            const polyPatch& topPatch() const
            {
                 return this->boundaryMesh()[topPatchID];
            }

            //- Return sides patch
            const polyPatch& sidesPatch() const
            {
                 return this->boundaryMesh()[sidePatchID];
            }

            //- Return extruded patch
            const polyPatch& extrudedPatch() const
            {
                return extrudedPatch_;
            }

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
