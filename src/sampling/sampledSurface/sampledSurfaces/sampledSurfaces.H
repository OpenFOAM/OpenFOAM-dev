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
        surfaces     | the list of surfaces    | yes         |
    \endtable

See also
    Foam::sampledSurfaces

SourceFiles
    sampledSurfaces.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_sampledSurfaces_H
#define functionObjects_sampledSurfaces_H

#include "functionObject.H"
#include "sampledSurface.H"
#include "surfaceWriter.H"
#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "wordReList.H"
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
    public functionObject,
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

        //- Output verbosity
        static bool verbose_;

        //- Tolerance for merging points (fraction of mesh bounding box)
        static scalar mergeTol_;


    // Private Data

        //- Const reference to fvMesh
        const fvMesh& mesh_;

        //- Load fields from files (not from objectRegistry)
        const bool loadFromFiles_;

        //- Output path
        fileName outputPath_;


        // Read from dictonary

            //- Names of fields to sample
            wordReList fieldSelection_;

            //- Interpolation scheme to use
            word interpolationScheme_;


        // surfaces

            //- Information for merging surfaces
            List<mergeInfo> mergeList_;


        // Calculated

            //- Surface formatter
            autoPtr<surfaceWriter> formatter_;


    // Private Member Functions


        //- Return number of fields
        label classifyFields();

        //- Write geometry only
        void writeGeometry() const;

        //- Write sampled fieldName on surface and on outputDir path
        template<class Type>
        void writeSurface
        (
            const Field<Type>& values,
            const label surfI,
            const word& fieldName,
            const fileName& outputDir
        );

        //- Sample and write a particular volume field
        template<class Type>
        void sampleAndWrite
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        );

        //- Sample and write a particular surface field
        template<class Type>
        void sampleAndWrite
        (
            const GeometricField<Type, fvsPatchField, surfaceMesh>&
        );

        //- Sample and write all sampled fields
        template<class Type> void sampleAndWrite(const IOobjectList& objects);


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

        //- Construct for given objectRegistry and dictionary
        //  allow the possibility to load fields from files
        sampledSurfaces
        (
            const word& name,
            const objectRegistry&,
            const dictionary&,
            const bool loadFromFiles = false
        );

        //- Disallow default bitwise copy construction
        sampledSurfaces(const sampledSurfaces&) = delete;


    //- Destructor
    virtual ~sampledSurfaces();


    // Member Functions

        //- Does any of the surfaces need an update?
        virtual bool needsUpdate() const;

        //- Mark the surfaces as needing an update.
        //  May also free up unneeded data.
        //  Return false if all surfaces were already marked as expired.
        virtual bool expire();

        //- Update the surfaces as required and merge surface points (parallel).
        //  Return false if no surfaces required an update.
        virtual bool update();

        //- Set verbosity level
        void verbose(const bool verbosity = true);

        //- Read the sampledSurfaces dictionary
        virtual bool read(const dictionary&);

        //- Execute, currently does nothing
        virtual bool execute();

        //- Sample and write
        virtual bool write();

        //- Update for changes of mesh - expires the surfaces
        virtual void updateMesh(const mapPolyMesh&);

        //- Update for mesh point-motion - expires the surfaces
        virtual void movePoints(const polyMesh&);

        //- Update for changes of mesh due to readUpdate - expires the surfaces
        virtual void readUpdate(const polyMesh::readUpdateState state);


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
