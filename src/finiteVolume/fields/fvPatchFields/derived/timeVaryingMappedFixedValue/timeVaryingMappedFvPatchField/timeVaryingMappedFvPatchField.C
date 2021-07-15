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

\*---------------------------------------------------------------------------*/

#include "timeVaryingMappedFvPatchField.H"
#include "Time.H"
#include "AverageField.H"
#include "IFstream.H"
#include "OSspecific.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
Foam::fileName
Foam::timeVaryingMappedFvPatchField<Type>::findFieldFile
(
    const word& timeName
) const
{
    const fileName fieldFileName
    (
        dataDir_/timeName/sampleName_/fieldTableName_
    );

    const fileName typeFieldFileName
    (
        dataDir_/timeName/sampleName_
       /pTraits<Type>::typeName + Field<Type>::typeName
       /fieldTableName_
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


template<class Type>
void Foam::timeVaryingMappedFvPatchField<Type>::checkTable()
{
    // Initialise
    if (mapperPtr_.empty())
    {
        // Reread values and interpolate
        const fileName samplePointsFile(dataDir_/pointsName_);

        pointField samplePoints((IFstream(samplePointsFile)()));

        if (debug)
        {
            Info<< "timeVaryingMappedFvPatchField :"
                << " Read " << samplePoints.size() << " sample points from "
                << samplePointsFile << endl;
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
                patch_.patch().faceCentres(),
                perturb_,
                nearestOnly
            )
        );

        // Read the times for which data is available
        sampleTimes_ = Time::findTimes(dataDir_);

        if (debug)
        {
            Info<< "timeVaryingMappedFvPatchField : In directory "
                << dataDir_ << " found times "
                << pointToPointPlanarInterpolation::timeNames(sampleTimes_)
                << endl;
        }
    }


    // Find current time in sampleTimes
    label lo = -1;
    label hi = -1;

    bool foundTime = mapperPtr_().findTime
    (
        sampleTimes_,
        startSampleTime_,
        time().value(),
        lo,
        hi
    );

    if (!foundTime)
    {
        FatalErrorInFunction
            << "Cannot find starting sampling values for current time "
            << time().value() << nl
            << "Have sampling values for times "
            << pointToPointPlanarInterpolation::timeNames(sampleTimes_) << nl
            << "In directory " <<  dataDir_ << " of field " << fieldTableName_
            << exit(FatalError);
    }


    // Update sampled data fields.

    if (lo != startSampleTime_)
    {
        startSampleTime_ = lo;

        if (startSampleTime_ == endSampleTime_)
        {
            // No need to reread since are end values
            if (debug)
            {
                Pout<< "checkTable : Setting startValues to (already read) "
                    << dataDir_/sampleTimes_[startSampleTime_].name()
                    << endl;
            }
            startSampledValues_ = endSampledValues_;
            startAverage_ = endAverage_;
        }
        else
        {
            if (debug)
            {
                Pout<< "checkTable : Reading startValues from "
                    << dataDir_/sampleTimes_[lo].name()
                    << endl;
            }

            // Reread values and interpolate
            const fileName valsFile
            (
                findFieldFile(sampleTimes_[startSampleTime_].name())
            );

            Field<Type> vals;

            if (setAverage_)
            {
                AverageField<Type> avals((IFstream(valsFile)()));
                vals = avals;
                startAverage_ = avals.average();
            }
            else
            {
                (IFstream(valsFile)()) >> vals;
            }

            if (vals.size() != mapperPtr_().sourceSize())
            {
                FatalErrorInFunction
                    << "Number of values (" << vals.size()
                    << ") differs from the number of points ("
                    <<  mapperPtr_().sourceSize()
                    << ") in file " << valsFile << exit(FatalError);
            }

            startSampledValues_ = mapperPtr_().interpolate(vals);
        }
    }

    if (hi != endSampleTime_)
    {
        endSampleTime_ = hi;

        if (endSampleTime_ == -1)
        {
            // endTime no longer valid. Might as well clear endValues.
            if (debug)
            {
                Pout<< "checkTable : Clearing endValues" << endl;
            }
            endSampledValues_.clear();
        }
        else
        {
            if (debug)
            {
                Pout<< "checkTable : Reading endValues from "
                    << dataDir_/sampleTimes_[endSampleTime_].name()
                    << endl;
            }

            // Reread values and interpolate
            const fileName valsFile
            (
                findFieldFile(sampleTimes_[endSampleTime_].name())
            );

            Field<Type> vals;

            if (setAverage_)
            {
                AverageField<Type> avals((IFstream(valsFile)()));
                vals = avals;
                endAverage_ = avals.average();
            }
            else
            {
                (IFstream(valsFile)()) >> vals;
            }

            if (vals.size() != mapperPtr_().sourceSize())
            {
                FatalErrorInFunction
                    << "Number of values (" << vals.size()
                    << ") differs from the number of points ("
                    <<  mapperPtr_().sourceSize()
                    << ") in file " << valsFile << exit(FatalError);
            }

            endSampledValues_ = mapperPtr_().interpolate(vals);
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::timeVaryingMappedFvPatchField<Type>::timeVaryingMappedFvPatchField
(
    const fvPatch& p,
    const word& fieldName
)
:
    patch_(p),
    fieldTableName_(word::null),
    dataDir_(time().constant()/"boundaryData"/p.name()),
    pointsName_("points"),
    sampleName_(word::null),
    setAverage_(false),
    perturb_(0),
    mapperPtr_(nullptr),
    sampleTimes_(0),
    startSampleTime_(-1),
    startSampledValues_(0),
    startAverage_(Zero),
    endSampleTime_(-1),
    endSampledValues_(0),
    endAverage_(Zero),
    offset_()
{}


template<class Type>
Foam::timeVaryingMappedFvPatchField<Type>::timeVaryingMappedFvPatchField
(
    const fvPatch& p,
    const dictionary& dict,
    const word& fieldName
)
:
    patch_(p),
    fieldTableName_(dict.lookupOrDefault("fieldTable", fieldName)),
    dataDir_
    (
        dict.lookupOrDefault
        (
            "dataDir",
            time().constant()/"boundaryData"/p.name()
        )
    ),
    pointsName_(dict.lookupOrDefault<fileName>("points", "points")),
    sampleName_(dict.lookupOrDefault("sample", word::null)),
    setAverage_(dict.lookupOrDefault("setAverage", false)),
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
    sampleTimes_(0),
    startSampleTime_(-1),
    startSampledValues_(0),
    startAverage_(Zero),
    endSampleTime_(-1),
    endSampledValues_(0),
    endAverage_(Zero),
    offset_()
{
    dataDir_.expand();
    pointsName_.expand();
    sampleName_.expand();

    if (dict.found("offset"))
    {
        offset_ = Function1<Type>::New("offset", dict);
    }

    if
    (
        mapMethod_ != "planarInterpolation"
     && mapMethod_ != "nearest"
    )
    {
        FatalIOErrorInFunction
        (
            dict
        )   << "mapMethod should be one of 'planarInterpolation'"
            << ", 'nearest'" << exit(FatalIOError);
    }
}


template<class Type>
Foam::timeVaryingMappedFvPatchField<Type>::
timeVaryingMappedFvPatchField
(
    const timeVaryingMappedFvPatchField<Type>& ptf
)
:
    patch_(ptf.patch_),
    fieldTableName_(ptf.fieldTableName_),
    dataDir_(ptf.dataDir_),
    pointsName_(ptf.pointsName_),
    sampleName_(ptf.sampleName_),
    setAverage_(ptf.setAverage_),
    perturb_(ptf.perturb_),
    mapMethod_(ptf.mapMethod_),
    mapperPtr_(nullptr),
    sampleTimes_(ptf.sampleTimes_),
    startSampleTime_(ptf.startSampleTime_),
    startSampledValues_(ptf.startSampledValues_),
    startAverage_(ptf.startAverage_),
    endSampleTime_(ptf.endSampleTime_),
    endSampledValues_(ptf.endSampledValues_),
    endAverage_(ptf.endAverage_),
    offset_(ptf.offset_, false)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::timeVaryingMappedFvPatchField<Type>::autoMap
(
    const fvPatchFieldMapper& m
)
{
    if (startSampledValues_.size())
    {
        m(startSampledValues_, startSampledValues_);
        m(endSampledValues_, endSampledValues_);
    }
    // Clear interpolator
    mapperPtr_.clear();
    startSampleTime_ = -1;
    endSampleTime_ = -1;
}


template<class Type>
void Foam::timeVaryingMappedFvPatchField<Type>::rmap
(
    const timeVaryingMappedFvPatchField<Type>& tiptf,
    const labelList& addr
)
{
    startSampledValues_.rmap(tiptf.startSampledValues_, addr);
    endSampledValues_.rmap(tiptf.endSampledValues_, addr);

    // Clear interpolator
    mapperPtr_.clear();
    startSampleTime_ = -1;
    endSampleTime_ = -1;
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::timeVaryingMappedFvPatchField<Type>::map()
{
    checkTable();

    // Interpolate between the sampled data

    tmp<Field<Type>> tfld(new Field<Type>(patch_.size()));
    Field<Type>& fld = tfld.ref();

    Type wantedAverage;

    if (endSampleTime_ == -1)
    {
        // Only start value
        if (debug)
        {
            Pout<< "updateCoeffs : Sampled, non-interpolated values"
                << " from start time:"
                << sampleTimes_[startSampleTime_].name() << nl;
        }

        fld = startSampledValues_;
        wantedAverage = startAverage_;
    }
    else
    {
        const scalar start = sampleTimes_[startSampleTime_].value();
        const scalar end = sampleTimes_[endSampleTime_].value();

        const scalar s = (time().value() - start)/(end - start);

        if (debug)
        {
            Pout<< "updateCoeffs : Sampled, interpolated values"
                << " between start time:"
                << sampleTimes_[startSampleTime_].name()
                << " and end time:" << sampleTimes_[endSampleTime_].name()
                << " with weight:" << s << endl;
        }

        fld = (1 - s)*startSampledValues_ + s*endSampledValues_;
        wantedAverage = (1 - s)*startAverage_ + s*endAverage_;
    }

    // Enforce average. Either by scaling (if scaling factor > 0.5) or by
    // offsetting.
    if (setAverage_)
    {
        Type averagePsi = gSum(patch_.magSf()*fld)/gSum(patch_.magSf());

        if (debug)
        {
            Pout<< "updateCoeffs :"
                << " actual average:" << averagePsi
                << " wanted average:" << wantedAverage
                << endl;
        }

        if (mag(averagePsi) < vSmall)
        {
            // Field too small to scale. Offset instead.
            const Type offset = wantedAverage - averagePsi;
            if (debug)
            {
                Pout<< "updateCoeffs :"
                    << " offsetting with:" << offset << endl;
            }
            fld += offset;
        }
        else
        {
            const scalar scale = mag(wantedAverage)/mag(averagePsi);

            if (debug)
            {
                Pout<< "updateCoeffs :"
                    << " scaling with:" << scale << endl;
            }
            fld *= scale;
        }
    }

    // Apply offset to mapped values
    if (offset_.valid())
    {
        const scalar t = time().timeOutputValue();
        fld += offset_->value(t);
    }

    if (debug)
    {
        Pout<< "updateCoeffs : set fixedValue to min:" << gMin(fld)
            << " max:" << gMax(fld)
            << " avg:" << gAverage(fld) << endl;
    }

    return tfld;
}


template<class Type>
void Foam::timeVaryingMappedFvPatchField<Type>::write
(
    Ostream& os
) const
{
    writeEntryIfDifferent
    (
        os,
        "dataDir",
        time().constant()/"boundaryData"/patch_.name(),
        dataDir_
    );

    writeEntryIfDifferent(os, "points", fileName("points"), pointsName_);
    writeEntryIfDifferent(os, "sample", fileName::null, sampleName_);
    writeEntryIfDifferent(os, "setAverage", Switch(false), setAverage_);
    writeEntryIfDifferent(os, "perturb", scalar(1e-5), perturb_);

    writeEntry(os, "fieldTable", fieldTableName_);

    writeEntryIfDifferent
    (
        os,
        "mapMethod",
        word("planarInterpolation"),
        mapMethod_
    );

    if (offset_.valid())
    {
        writeEntry(os, offset_());
    }
}


// ************************************************************************* //
