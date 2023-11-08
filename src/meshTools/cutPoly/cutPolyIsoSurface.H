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
    Foam::cutPolyIsoSurface

Description
    Iso-surface class based on the cutPoly cutting routines

SourceFiles
    cutPolyIsoSurface.C
    cutPolyIsoSurfaceTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef cutPolyIsoSurface_H
#define cutPolyIsoSurface_H

#include "polyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class cutPolyIsoSurface Declaration
\*---------------------------------------------------------------------------*/

class cutPolyIsoSurface
{
    // Private Data

        //- Points of the iso-surface
        pointField points_;

        //- The mesh edge associated with each point
        edgeList pointEdges_;

        //- The local coordinate of each point within its associated edge
        scalarList pointEdgeLambdas_;

        //- Faces of the iso-surface
        faceList faces_;

        //- The mesh cell index associated with each face
        labelList faceCells_;


public:

    //- Runtime type information
    TypeName("cutPolyIsoSurface");


    // Constructors

        //- Construct from a mesh, point values and an iso-value
        cutPolyIsoSurface
        (
            const polyMesh& mesh,
            const scalarField& pAlphas,
            const scalar isoAlpha,
            const word& cellZoneName = word::null
        );

        //- Construct by combining a list of iso-surfaces
        cutPolyIsoSurface(const PtrList<cutPolyIsoSurface>& isos);

        //- Disallow default bitwise copy construction
        cutPolyIsoSurface(const cutPolyIsoSurface&) = delete;


    //- Destructor
    virtual ~cutPolyIsoSurface();


    // Member Functions

        // Access

            //- Points of the iso-surface
            inline const pointField& points() const
            {
                return points_;
            }

            //- Faces of the iso-surface
            inline const faceList& faces() const
            {
                return faces_;
            }

            //- The mesh cell index associated with each face
            inline const labelList& faceCells() const
            {
                return faceCells_;
            }


        //- Return a field of values for the surface faces, sampled from their
        //  associated cells
        template<class Type>
        tmp<Field<Type>> sample(const Field<Type> cPsis) const;

        //- Return a field of values for the surface points, interpolated
        //  within their associated edges
        template<class Type>
        tmp<Field<Type>> interpolate(const Field<Type>& pPsis) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "cutPolyIsoSurfaceTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
