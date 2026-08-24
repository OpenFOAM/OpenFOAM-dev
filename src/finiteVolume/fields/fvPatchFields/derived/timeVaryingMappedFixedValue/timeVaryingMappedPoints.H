/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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
    Foam::timeVaryingMappedPoints

Description
    Provides the sample points, times and field files from which time-varying
    mapped data is interpolated, and constructs the corresponding
    point-to-point mapper.

    The points are read from a static points file in the data directory, or,
    if the sampled surface changes in time, from per-time points files in the
    time directories which take precedence over the static file. For a given
    sample time the most recent points file at or before that time is used.

    The mapper is only reconstructed when the selected points file changes.

See also
    Foam::timeVaryingMappedFvPatchField

SourceFiles
    timeVaryingMappedPoints.C

\*---------------------------------------------------------------------------*/

#ifndef timeVaryingMappedPoints_H
#define timeVaryingMappedPoints_H

#include "instantList.H"
#include "pointFieldFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class Time;
class pointToPointPlanarInterpolation;

/*---------------------------------------------------------------------------*\
                   Class timeVaryingMappedPoints Declaration
\*---------------------------------------------------------------------------*/

class timeVaryingMappedPoints
{
    // Private Data

        //- Optional name of the data directory,
        //  defaults to "constant/boundaryData/<patch name>"
        fileName dataDir_;

        //- Optional name of the points file, defaults to "points"
        fileName pointsName_;

        //- Optional name of the sample sub-directory from which the fields
        //  are read, defaults to ""
        fileName sampleName_;

        //- Fraction of perturbation (fraction of bounding box) to add
        scalar perturb_;

        //- Interpolation scheme to use
        word mapMethod_;

        //- List of the data directory time directories
        instantList sampleTimes_;

        //- 2D interpolation (for 'planarInterpolation' mapMethod)
        autoPtr<pointToPointPlanarInterpolation> mapperPtr_;

        //- Path of the points file from which the mapper was constructed
        fileName mapperPointsFile_;

        //- Highest sample time index already scanned for a points file
        label pointsSampleIndex_;

        //- Points file found by the most recent scan
        fileName pointsScanFile_;


    // Private Member Functions

        //- Find the most recent points file at or before the given sample
        //  time index, falling back to the static points file
        fileName findPointsFile(const label sampleIndex);


public:

    // Constructors

        //- Construct from dictionary, time and patch name
        timeVaryingMappedPoints
        (
            const dictionary&,
            const Time&,
            const word& patchName
        );

        //- Copy constructor. Copies the settings and sample times but not
        //  the mapper and points file scan caches
        timeVaryingMappedPoints(const timeVaryingMappedPoints&);


    // Member Functions

        //- Return the name of the data directory
        const fileName& dataDir() const
        {
            return dataDir_;
        }

        //- Return the sample times, reading the time directories on first
        //  use and after clear()
        const instantList& sampleTimes(const Time&);

        //- Find the field file for the given sample time index, checking
        //  the type sub-directory (e.g. "vectorField") as a fallback
        fileName findFieldFile
        (
            const label sampleIndex,
            const word& fieldTableName,
            const string& fieldTypeDirName
        ) const;

        //- Return the mapper for the given sample time, updating it if the
        //  points file differs from that from which it was constructed
        const pointToPointPlanarInterpolation& mapper
        (
            const label sampleIndex,
            const pointField& destPoints
        );

        //- Clear the sample times, the mapper and the points file scan cache
        void clear();

        //- Write the non-default dictionary entries
        void write(Ostream&, const Time&, const word& patchName) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const timeVaryingMappedPoints&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
