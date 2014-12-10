/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
    Foam::wallBoundedStreamLine

Group
    grpFieldFunctionObjects

Description
    This function object generates streamline data by sampling a set of
    user-specified fields along a particle track, transported by a
    user-specified velocity field, constrained to a patch.

    Example of function object specification:
    \verbatim
    wallBoundedStreamLine1
    {
        type            wallBoundedStreamLine;
        functionObjectLibs ("libfieldFunctionObjects.so");
        ...
        setFormat       vtk;
        UName           UNear;
        trackForward    yes;
        fields
        (
            UNear
            p
        );
        lifeTime        10000;
        trackLength     1e-3;
        nSubCycle       5;
        cloudName       particleTracks;
        seedSampleSet   patchSeed;
        patchSeedCoeffs
        {
            type        patchSeed;
            patches     (wall);
            axis        x;
            maxPoints   20000;
        }
    }
    \endverbatim

    \heading Function object usage
    \table
        Property     | Description             | Required    | Default value
        type         | type name: wallBoundedStreamLine| yes |
        setFormat    | output data type        | yes         |
        UName        | tracking velocity field name | yes    |
        fields       | fields to sample        | yes         |
        lifetime     | maximum number of particle tracking steps | yes |
        trackLength  | tracking segment length | no          |
        nSubCycle    | number of tracking steps per cell | no|
        cloudName    | cloud name to use       | yes         |
        seedSampleSet| seeding method (see below)| yes       |
    \endtable

    \linebreak
    Where \c seedSampleSet is typically one of
    \plaintable
        uniform | uniform particle seeding
        cloud   | cloud of points
        patchSeed | seeding via patch faces
        triSurfaceMeshPointSet | points according to a tri-surface mesh
    \endplaintable

Note
    When specifying the track resolution, the \c trackLength OR \c nSubCycle
    option should be used

SeeAlso
    Foam::functionObject
    Foam::OutputFilterFunctionObject
    Foam::sampledSet
    Foam::streamLine

SourceFiles
    wallBoundedStreamLine.C

\*---------------------------------------------------------------------------*/

#ifndef wallBoundedStreamLine_H
#define wallBoundedStreamLine_H

#include "volFieldsFwd.H"
#include "pointFieldFwd.H"
#include "Switch.H"
#include "DynamicList.H"
#include "scalarList.H"
#include "vectorList.H"
#include "polyMesh.H"
#include "writer.H"
#include "indirectPrimitivePatch.H"
#include "tetIndices.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class objectRegistry;
class dictionary;
class mapPolyMesh;
class meshSearch;
class sampledSet;

/*---------------------------------------------------------------------------*\
                    Class wallBoundedStreamLine Declaration
\*---------------------------------------------------------------------------*/

class wallBoundedStreamLine
{
    // Private data

        //- Input dictionary
        dictionary dict_;

        //- Name of this set of field averages.
        word name_;

        //- Database this class is registered to
        const objectRegistry& obr_;

        //- Load fields from files (not from objectRegistry)
        bool loadFromFiles_;

        //- On/off switch
        bool active_;


        //- List of fields to sample
        wordList fields_;

        //- Field to transport particle with
        word UName_;

        //- Interpolation scheme to use
        word interpolationScheme_;

        //- Whether to use +u or -u
        bool trackForward_;

        //- Maximum lifetime (= number of cells) of particle
        label lifeTime_;

        //- Track length
        scalar trackLength_;

        //- Optional specified name of particles
        word cloudName_;

        //- Type of seed
        word seedSet_;

        //- Names of scalar fields
        wordList scalarNames_;

        //- Names of vector fields
        wordList vectorNames_;


        // Demand driven

            //- Mesh searching enigne
            autoPtr<meshSearch> meshSearchPtr_;

            //- Seed set engine
            autoPtr<sampledSet> sampledSetPtr_;

            //- Axis of the sampled points to output
            word sampledSetAxis_;

            //- File output writer
            autoPtr<writer<scalar> > scalarFormatterPtr_;

            autoPtr<writer<vector> > vectorFormatterPtr_;


        // Generated data

            //- All tracks. Per particle the points it passed through
            DynamicList<List<point> > allTracks_;

            //- Per scalarField, per particle, the sampled value.
            List<DynamicList<scalarList> > allScalars_;

            //- Per scalarField, per particle, the sampled value.
            List<DynamicList<vectorList> > allVectors_;


        //- Construct patch out of all wall patch faces
        autoPtr<indirectPrimitivePatch> wallPatch() const;

        //- Find wall tet on cell
        tetIndices findNearestTet
        (
            const PackedBoolList& isWallPatch,
            const point& seedPt,
            const label cellI
        ) const;

        //- Do all seeding and tracking
        void track();

        //- Disallow default bitwise copy construct
        wallBoundedStreamLine(const wallBoundedStreamLine&);

        //- Disallow default bitwise assignment
        void operator=(const wallBoundedStreamLine&);


public:

    //- Runtime type information
    TypeName("wallBoundedStreamLine");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        wallBoundedStreamLine
        (
            const word& name,
            const objectRegistry&,
            const dictionary&,
            const bool loadFromFiles = false
        );


    //- Destructor
    virtual ~wallBoundedStreamLine();


    // Member Functions

        //- Return name of the set of field averages
        virtual const word& name() const
        {
            return name_;
        }

        //- Read the field average data
        virtual void read(const dictionary&);

        //- Execute the averaging
        virtual void execute();

        //- Execute the averaging at the final time-loop, currently does nothing
        virtual void end();

        //- Called when time was set at the end of the Time::operator++
        virtual void timeSet();

        //- Calculate the field average data and write
        virtual void write();

        //- Update for changes of mesh
        virtual void updateMesh(const mapPolyMesh&);

        //- Update for mesh point-motion
        virtual void movePoints(const polyMesh&);

        ////- Update for changes of mesh due to readUpdate
        //virtual void readUpdate(const polyMesh::readUpdateState state);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
