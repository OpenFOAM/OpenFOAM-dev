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

\*---------------------------------------------------------------------------*/

#include "timeVaryingMappedPoints.H"
#include "timeVaryingMappedFvPatchField.H"
#include "pointToPointPlanarInterpolation.H"
#include "Time.H"
#include "IFstream.H"
#include "OSspecific.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::fileName Foam::timeVaryingMappedPoints::findPointsFile
(
    const label sampleIndex
)
{
    const label startIndex =
        sampleIndex < pointsSampleIndex_ ? 0 : pointsSampleIndex_ + 1;
    for (label i = sampleIndex; i >= startIndex; i--)
    {
        const fileName timePointsFile
        (
            dataDir_/sampleTimes_[i].name()/sampleName_/pointsName_.name()
        );

        if (isFile(timePointsFile))
        {
            pointsSampleIndex_ = sampleIndex;
            pointsScanFile_ = timePointsFile;

            return pointsScanFile_;
        }
    }

    if (sampleIndex >= pointsSampleIndex_ && !pointsScanFile_.empty())
    {
        pointsSampleIndex_ = sampleIndex;

        return pointsScanFile_;
    }

    // Static points file name used for all sample times
    const fileName pointsFile(dataDir_/pointsName_);
    if (isFile(pointsFile))
    {
        pointsSampleIndex_ = sampleIndex;
        pointsScanFile_ = pointsFile;

        return pointsScanFile_;
    }
    else
    {
        FatalErrorInFunction
            << "Cannot find points file "
            << dataDir_/sampleTimes_[sampleIndex].name()
              /sampleName_/pointsName_.name()
            << " " << pointsFile
            << exit(FatalError);

        return fileName::null;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::timeVaryingMappedPoints::timeVaryingMappedPoints
(
    const dictionary& dict,
    const Time& time,
    const word& patchName
)
:
    dataDir_
    (
        dict.lookupOrDefault
        (
            "dataDir",
            time.constant()/"boundaryData"/patchName
        )
    ),
    pointsName_(dict.lookupOrDefault<fileName>("points", "points")),
    sampleName_(dict.lookupOrDefault("sample", word::null)),
    perturb_(dict.lookupOrDefault("perturb", 1e-5)),
    mapMethod_
    (
        dict.lookupOrDefault<word>
        (
            "mapMethod",
            "planarInterpolation"
        )
    ),
    mapperPtr_(nullptr),
    pointsSampleIndex_(-1)
{
    dataDir_.expand();
    pointsName_.expand();
    sampleName_.expand();

    if
    (
        mapMethod_ != "planarInterpolation"
     && mapMethod_ != "nearest"
    )
    {
        FatalIOErrorInFunction(dict)
            << "mapMethod should be one of 'planarInterpolation'"
            << ", 'nearest'" << exit(FatalIOError);
    }
}


Foam::timeVaryingMappedPoints::timeVaryingMappedPoints
(
    const timeVaryingMappedPoints& mp
)
:
    dataDir_(mp.dataDir_),
    pointsName_(mp.pointsName_),
    sampleName_(mp.sampleName_),
    perturb_(mp.perturb_),
    mapMethod_(mp.mapMethod_),
    sampleTimes_(mp.sampleTimes_),
    mapperPtr_(nullptr),
    pointsSampleIndex_(-1)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::instantList& Foam::timeVaryingMappedPoints::sampleTimes
(
    const Time& time
)
{
    if (sampleTimes_.empty())
    {
        sampleTimes_ = time.findTimes(dataDir_);

        if (timeVaryingMappedFvPatchFieldName::debug)
        {
            Info<< "timeVaryingMappedFvPatchField : In directory "
                << dataDir_ << " found times "
                << pointToPointPlanarInterpolation::timeNames(sampleTimes_)
                << endl;
        }
    }

    return sampleTimes_;
}


Foam::fileName Foam::timeVaryingMappedPoints::findFieldFile
(
    const label sampleIndex,
    const word& fieldTableName,
    const string& fieldTypeDirName
) const
{
    const word& timeName = sampleTimes_[sampleIndex].name();

    const fileName fieldFileName
    (
        dataDir_/timeName/sampleName_/fieldTableName
    );

    const fileName typeFieldFileName
    (
        dataDir_/timeName/sampleName_/fieldTypeDirName/fieldTableName
    );

    if (exists(fieldFileName))
    {
        return fieldFileName;
    }
    else if (exists(typeFieldFileName))
    {
        return typeFieldFileName;
    }
    else
    {
        FatalErrorInFunction
            << "Cannot find field file "
            << fieldFileName << " " << typeFieldFileName
            << exit(FatalError);

        return fileName::null;
    }
}


const Foam::pointToPointPlanarInterpolation&
Foam::timeVaryingMappedPoints::mapper
(
    const label sampleIndex,
    const pointField& destPoints
)
{
    const fileName pointsFile(findPointsFile(sampleIndex));

    // If the points file has changed update the mapper
    if (mapperPtr_.empty() || mapperPointsFile_ != pointsFile)
    {
        // Reread the sample points
        const pointField samplePoints((IFstream(pointsFile)()));

        if (timeVaryingMappedFvPatchFieldName::debug)
        {
            Info<< "timeVaryingMappedFvPatchField :"
                << " Read " << samplePoints.size() << " sample points from "
                << pointsFile << endl;
        }

        // tbd: run-time selection
        const bool nearestOnly
        (
           !mapMethod_.empty()
         && mapMethod_ != "planarInterpolation"
        );

        // Allocate the interpolator
        mapperPtr_.reset
        (
            new pointToPointPlanarInterpolation
            (
                samplePoints,
                destPoints,
                perturb_,
                nearestOnly
            )
        );

        mapperPointsFile_ = pointsFile;
    }

    return mapperPtr_();
}


void Foam::timeVaryingMappedPoints::clear()
{
    sampleTimes_.clear();
    mapperPtr_.clear();
    mapperPointsFile_.clear();
    pointsSampleIndex_ = -1;
    pointsScanFile_.clear();
}


void Foam::timeVaryingMappedPoints::write
(
    Ostream& os,
    const Time& time,
    const word& patchName
) const
{
    writeEntryIfDifferent
    (
        os,
        "dataDir",
        time.constant()/"boundaryData"/patchName,
        dataDir_
    );

    writeEntryIfDifferent(os, "points", fileName("points"), pointsName_);
    writeEntryIfDifferent(os, "sample", fileName::null, sampleName_);
    writeEntryIfDifferent(os, "perturb", scalar(1e-5), perturb_);

    writeEntryIfDifferent
    (
        os,
        "mapMethod",
        word("planarInterpolation"),
        mapMethod_
    );
}


// ************************************************************************* //
