/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    Foam::rawTopoChangerFvMesh

Description
    topoChangerFvMesh without any added functionality.

    Note: run without FOAM_SETNAN. Temporary has unitialised patch
    data when faces get created out of nothing.

SourceFiles
    rawTopoChangerFvMesh.C

\*---------------------------------------------------------------------------*/

#ifndef rawTopoChangerFvMesh_H
#define rawTopoChangerFvMesh_H

#include "topoChangerFvMesh.H"
#include "PackedBoolList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

/*---------------------------------------------------------------------------*\
                      Class rawTopoChangerFvMesh Declaration
\*---------------------------------------------------------------------------*/

class rawTopoChangerFvMesh
:
    public topoChangerFvMesh
{
    // Private Member Functions

        //- Set unmapped values
        template<class Type, template<class> class PatchField, class GeoMesh>
        static void setUnmappedValues
        (
            GeometricField<Type, PatchField, GeoMesh>& fld,
            const PackedBoolList& mappedFace,
            const GeometricField<Type, PatchField, GeoMesh>& baseFld
        );

        template<class Type, template<class> class PatchField, class GeoMesh>
        void zeroUnmappedValues(const PackedBoolList&) const;

        //- Disallow default bitwise copy construct
        rawTopoChangerFvMesh(const rawTopoChangerFvMesh&);

        //- Disallow default bitwise assignment
        void operator=(const rawTopoChangerFvMesh&);

public:

    //- Runtime type information
    TypeName("rawTopoChangerFvMesh");


    // Constructors

        //- Construct from database
        explicit rawTopoChangerFvMesh(const IOobject& io);

    //- Destructor
    virtual ~rawTopoChangerFvMesh();


    // Member Functions

        //- Update the mesh for both mesh motion and topology change
        virtual bool update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "rawTopoChangerFvMeshTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
