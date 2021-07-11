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
    Foam::functionObjects::streamlines

Description
    Generates streamline data by sampling a set of user-specified fields along
    a particle track, transported by a user-specified velocity field.

    Example of function object specification:
    \verbatim
    streamlines1
    {
        type            streamlines;
        libs            ("libfieldFunctionObjects.so");

        writeControl    writeTime;

        setFormat       vtk;
        U               U;
        direction       both;

        fields
        (
            U
            p
        );

        lifeTime        10000;
        trackLength     1e-3;
        nSubCycle       5;
        cloudName       particleTracks;

        seedSampleSet
        {
            type        uniform;
            axis        x;  // distance;
            start       (-0.0205 0.0001 0.00001);
            end         (-0.0205 0.0005 0.00001);
            nPoints     100;
        }
    }
    \endverbatim

Usage
    \table
        Property     | Description             | Required    | Default value
        type         | Type name: streamlines  | yes         |
        setFormat    | Output data type        | yes         |
        U            | Tracking velocity field name | no     | U
        direction    | Direction in which to track | yes     |
        fields       | Fields to sample        | yes         |
        writeTime    | Write the flow time along the streamlines | no | no
        lifetime     | Maximum number of particle tracking steps | yes |
        trackLength  | Tracking segment length | no          |
        nSubCycle    | Number of tracking steps per cell | no |
        cloudName    | Cloud name to use       | yes         |
        seedSampleSet| Seeding method (see below)| yes       |
    \endtable

    Where \c seedSampleSet \c type is typically one of
    \plaintable
        uniform | uniform particle seeding
        cloud   | cloud of points
        triSurfaceMeshPointSet | points according to a tri-surface mesh
    \endplaintable

Note
    When specifying the track resolution, the \c trackLength OR \c nSubCycle
    option should be used

See also
    Foam::functionObject
    Foam::functionObjects::timeControl
    Foam::sampledSet

SourceFiles
    streamlines.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_streamlines_H
#define functionObjects_streamlines_H

#include "fvMeshFunctionObject.H"
#include "volFieldsFwd.H"
#include "DynamicList.H"
#include "scalarList.H"
#include "vectorList.H"
#include "setWriter.H"
#include "indirectPrimitivePatch.H"
#include "NamedEnum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class meshSearch;
class sampledSet;

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                         Class streamlines Declaration
\*---------------------------------------------------------------------------*/

class streamlines
:
    public fvMeshFunctionObject
{
public:

    // Public data types

        //- Track direction enumerations
        enum class trackDirection
        {
            forward,
            backward,
            both
        };

        //- Track direction enumeration names
        static const NamedEnum<trackDirection, 3> trackDirectionNames_;


private:

    // Private Data

        //- Input dictionary
        dictionary dict_;

        //- List of fields to sample
        wordList fields_;

        //- Field to transport particle with
        word UName_;

        //- Interpolation scheme to use
        word interpolationScheme_;

        //- The direction in which to track
        trackDirection trackDirection_;

        //- Maximum lifetime (= number of cells) of particle
        label lifeTime_;

        //- Number of subcycling steps
        label nSubCycle_;

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

        //- Write the streamlines ages
        Switch writeAge_;


        // Demand driven

            //- Mesh searching engine
            autoPtr<meshSearch> meshSearchPtr_;

            //- Seed set engine
            autoPtr<sampledSet> sampledSetPtr_;

            //- Axis of the sampled points to output
            word sampledSetAxis_;

            //- File writer for scalar data
            autoPtr<setWriter<scalar>> scalarFormatterPtr_;

            //- File writer for vector data
            autoPtr<setWriter<vector>> vectorFormatterPtr_;


        // Generated data

            //- All tracks. Per particle the points it passed through
            DynamicList<List<point>> allTracks_;

            //- All ages. Per particle the age when it passed through the points
            DynamicList<List<scalar>> allAges_;

            //- Per scalarField, per particle, the sampled value.
            List<DynamicList<scalarList>> allScalars_;

            //- Per scalarField, per particle, the sampled value.
            List<DynamicList<vectorList>> allVectors_;


        //- Construct patch out of all wall patch faces
        autoPtr<indirectPrimitivePatch> wallPatch() const;

        //- Do all seeding and tracking
        void track();


public:

    //- Runtime type information
    TypeName("streamlines");


    // Constructors

        //- Construct from Time and dictionary
        streamlines
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        streamlines(const streamlines&) = delete;


    //- Destructor
    virtual ~streamlines();


    // Member Functions

        //- Read the field average data
        virtual bool read(const dictionary&);

        //- Do nothing
        virtual bool execute();

        //- Calculate and write the streamlines
        virtual bool write();

        //- Update for changes of mesh
        virtual void updateMesh(const mapPolyMesh&);

        //- Update for mesh point-motion
        virtual void movePoints(const polyMesh&);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const streamlines&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
