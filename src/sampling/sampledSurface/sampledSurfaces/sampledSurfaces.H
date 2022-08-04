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
    Foam::functionObjects::sampledSurfaces

Description
    Set of surfaces to sample

    Example of function object specification:
    \verbatim
    surfaces1
    {
        type            surfaces;
        libs            ("libsampling.so");

        writeControl    writeTime;

        fields          (p U);

        surfaceFormat   vtk;
        interpolationScheme cellPoint;

        surfaces
        (
            p1e5
            {
                type            isoSurface;
                isoField        p;
                isoValue        1e5;
                interpolate     yes;
            }
        );
    }
    \endverbatim

Usage
    \table
        Property     | Description             | Required    | Default value
        type         | type name: surfaces     | yes         |
        surfaceFormat | the format in which to write the surface file | yes   |
        interpolationScheme | the method by which values are interpolated \\
                              from the mesh to the surface | yes
        writeEmpty   | write out files for empty surfaces | no | no
        surfaces     | the list of surfaces    | yes         |
    \endtable

See also
    Foam::sampledSurfaces

SourceFiles
    sampledSurfaces.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_sampledSurfaces_H
#define functionObjects_sampledSurfaces_H

#include "fvMeshFunctionObject.H"
#include "sampledSurface.H"
#include "surfaceWriter.H"
#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "IOobjectList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class Time;
class fvMesh;
class dictionary;

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class sampledSurfaces Declaration
\*---------------------------------------------------------------------------*/

class sampledSurfaces
:
    public fvMeshFunctionObject,
    public PtrList<sampledSurface>
{
    // Private classes

        //- Class used for surface merging information
        class mergeInfo
        {
        public:
            pointField points;
            faceList   faces;
            labelList  pointsMap;

            //- Clear all storage
            void clear()
            {
                points.clear();
                faces.clear();
                pointsMap.clear();
            }
        };


    // Static Data Members

        //- Tolerance for merging points (fraction of mesh bounding box)
        static scalar mergeTol_;


    // Private Data

        //- Output path
        fileName outputPath_;


        // Read from dictionary

            //- Names of fields to sample
            wordList fields_;

            //- Interpolation scheme to use
            word interpolationScheme_;

            //- Should we create files for empty surfaces?
            Switch writeEmpty_;


        // Surfaces

            //- Information for merging surfaces
            List<mergeInfo> mergeList_;


        // Calculated

            //- Surface formatter
            autoPtr<surfaceWriter> formatter_;


    // Private Member Functions

        //- Does any of the surfaces need an update?
        bool needsUpdate() const;

        //- Mark the surfaces as needing an update.
        //  May also free up unneeded data.
        //  Return false if all surfaces were already marked as expired.
        bool expire();

        //- Update the surfaces as required and merge surface points (parallel).
        //  Return false if no surfaces required an update.
        bool update();

        //- Sample all fields of a type on a given surface
        template<class Type>
        PtrList<Field<Type>> sampleLocalType
        (
            const label surfi,
            const wordList& fieldNames,
            HashPtrTable<interpolation<Type>>& interpolations
        );

        //- Sample all fields of a type on a given surface and combine on the
        //  master for writing
        template<class Type>
        PtrList<Field<Type>> sampleType
        (
            const label surfi,
            const wordList& fieldNames,
            HashPtrTable<interpolation<Type>>& interpolations
        );


public:

    //- Runtime type information
    TypeName("surfaces");


    // Constructors

        //- Construct from Time and dictionary
        sampledSurfaces
        (
            const word& name,
            const Time& time,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        sampledSurfaces(const sampledSurfaces&) = delete;


    //- Destructor
    virtual ~sampledSurfaces();


    // Member Functions

        //- Read the sampledSurfaces dictionary
        virtual bool read(const dictionary&);

        //- Return the list of fields required
        virtual wordList fields() const;

        //- Execute, currently does nothing
        virtual bool execute();

        //- Sample and write
        virtual bool write();

        //- Update for mesh point-motion - expires the surfaces
        virtual void movePoints(const polyMesh&);

        //- Update topology using the given map - expires the surfaces
        virtual void topoChange(const polyTopoChangeMap&);

        //- Update from another mesh using the given map
        virtual void mapMesh(const polyMeshMap&);

        //- Redistribute or update using the given distribution map
        virtual void distribute(const polyDistributionMap&);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const sampledSurfaces&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "sampledSurfacesTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
