/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2025 OpenFOAM Foundation
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
    Foam::searchableSurfaces::closedTriSurface

Description
    Like triSurface, a surface geometry formed of discrete facets, e.g. \
    triangles and/or quadrilaterals, defined in a file using formats such as
    Wavefront OBJ, or stereolithography STL. A closedTriSurface is a
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
        type       closedTriSurface;
        file       "surfaceFile.obj";
    \endverbatim

SourceFiles
    closedTriSurface.C

\*---------------------------------------------------------------------------*/

#ifndef closedTriSurface_searchableSurface_H
#define closedTriSurface_searchableSurface_H

#include "triSurface_searchableSurface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace searchableSurfaces
{

/*---------------------------------------------------------------------------*\
                    Class closedTriSurface Declaration
\*---------------------------------------------------------------------------*/

class closedTriSurface
:
    public triSurface
{
public:

    //- Runtime type information
    TypeName("closedTriSurface");


    // Constructors

        //- Construct from triSurface
        closedTriSurface(const IOobject&, const triSurface&);

        //- Construct read.
        closedTriSurface(const IOobject& io);

        //- Construct from IO and dictionary (used by searchableSurface).
        //  Dictionary may contain a 'scale' entry (eg, 0.001: mm -> m)
        closedTriSurface
        (
            const IOobject& io,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        closedTriSurface(const closedTriSurface&) = delete;


    // Destructor

        virtual ~closedTriSurface();


    // Member Functions

        //- Whether supports volume type, forcing to true to force getVolumeType
        //  queries for this type
        virtual bool hasVolumeType() const
        {
            return true;
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const closedTriSurface&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End of namespace searchableSurfaces
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
