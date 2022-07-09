/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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
    Foam::closedTriSurfaceMesh

Description
    Like triSurface, a surface geometry formed of discrete facets, e.g. \
    triangles and/or quadrilaterals, defined in a file using formats such as
    Wavefront OBJ, or stereolithography STL. A closedTriSurfaceMesh is a
    geometry surface which is meant to be closed but contains some
    imperfections, e.g. small holes or multiple parts, which mean it is not
    strictly closed.

Usage
    \table
        Property     | Description                | Required
        file         | Name of the geometry file  | yes
        scale        | Scaling factor for surface | no
        minQuality   | Threshold triangle quality | no
    \endtable

    Note: when calculating surface normal vectors, triangles are ignored with
    quality < minQuality.

    Example specification in snappyHexMeshDict/geometry:
    \verbatim
        type       closedTriSurfaceMesh;
        file       "surfaceFile.obj";
    \endverbatim

SourceFiles
    closedTriSurfaceMesh.C

\*---------------------------------------------------------------------------*/

#ifndef closedTriSurfaceMesh_H
#define closedTriSurfaceMesh_H

#include "triSurfaceMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class closedTriSurfaceMesh Declaration
\*---------------------------------------------------------------------------*/

class closedTriSurfaceMesh
:
    public triSurfaceMesh
{
public:

    //- Runtime type information
    TypeName("closedTriSurfaceMesh");


    // Constructors

        //- Construct from triSurface
        closedTriSurfaceMesh(const IOobject&, const triSurface&);

        //- Construct read.
        closedTriSurfaceMesh(const IOobject& io);

        //- Construct from IO and dictionary (used by searchableSurface).
        //  Dictionary may contain a 'scale' entry (eg, 0.001: mm -> m)
        closedTriSurfaceMesh
        (
            const IOobject& io,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        closedTriSurfaceMesh(const closedTriSurfaceMesh&) = delete;


    // Destructor

        virtual ~closedTriSurfaceMesh();


    // Member Functions

        //- Whether supports volume type, forcing to true to force getVolumeType
        //  queries for this type
        virtual bool hasVolumeType() const
        {
            return true;
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const closedTriSurfaceMesh&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
