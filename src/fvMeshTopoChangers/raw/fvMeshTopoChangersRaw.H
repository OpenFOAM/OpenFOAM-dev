/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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
    Foam::fvMeshTopoChangers::raw

Description
    fvMeshTopoChanger without any added functionality.

    Note: run without FOAM_SETNAN. Temporary has unitialised patch
    data when faces get created out of nothing.

SourceFiles
    fvMeshTopoChangersRaw.C

\*---------------------------------------------------------------------------*/

#ifndef fvMeshTopoChangersRaw_H
#define fvMeshTopoChangersRaw_H

#include "fvMeshTopoChanger.H"
#include "polyTopoChanger.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fvMeshTopoChangers
{

/*---------------------------------------------------------------------------*\
                   Class fvMeshTopoChangers::raw Declaration
\*---------------------------------------------------------------------------*/

class raw
:
    public fvMeshTopoChanger
{
    // Private Data

        polyTopoChanger topoChanger_;


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


public:

    //- Runtime type information
    TypeName("raw");


    // Constructors

        //- Construct from fvMesh
        explicit raw(fvMesh& io);

        //- Disallow default bitwise copy construction
        raw(const raw&) = delete;


    //- Destructor
    virtual ~raw();


    // Member Functions

        //- Update the mesh for both mesh motion and topology change
        virtual bool update();

        //- Update corresponding to the given map
        virtual void updateMesh(const mapPolyMesh&);

        //- Update corresponding to the given distribution map
        virtual void distribute(const mapDistributePolyMesh&);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const raw&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvMeshTopoChangers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fvMeshTopoChangersRawTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
