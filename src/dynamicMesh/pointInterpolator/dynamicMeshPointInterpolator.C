/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2019 OpenFOAM Foundation
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

#include "dynamicMeshPointInterpolator.H"
#include "pointFields.H"
#include "interpolationWeights.H"
#include "uniformInterpolate.H"
#include "ReadFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dynamicMeshPointInterpolator::dynamicMeshPointInterpolator
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    mesh_(mesh),
    fieldName_(dict.lookup("field")),
    interpolationScheme_(dict.lookup("interpolationScheme"))
{
    const pointMesh& pMesh = pointMesh::New(mesh_);

    // Read time values
    instantList allTimes = Time::findTimes(pMesh().time().path());

    // Only keep those that contain the field
    DynamicList<word> names(allTimes.size());
    DynamicList<scalar> values(allTimes.size());

    forAll(allTimes, i)
    {
        IOobject io
        (
            fieldName_,
            allTimes[i].name(),
            pMesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        );
        if (io.typeHeaderOk<pointVectorField>(false))
        {
            names.append(allTimes[i].name());
            values.append(allTimes[i].value());
        }
    }
    timeNames_.transfer(names);
    timeVals_.transfer(values);

    Info<< mesh_.type() << " : found " << fieldName_ << " for times "
        << timeNames_ << endl;

    if (timeNames_.size() < 1)
    {
        FatalErrorInFunction
            << "Did not find any times with " << fieldName_
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dynamicMeshPointInterpolator::~dynamicMeshPointInterpolator()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::pointVectorField>
Foam::dynamicMeshPointInterpolator::curPointField() const
{
    if (!interpolatorPtr_.valid())
    {
        interpolatorPtr_ = interpolationWeights::New
        (
            interpolationScheme_,
            timeVals_
        );
    }

    const pointMesh& pMesh = pointMesh::New(mesh_);

    const Time& t = pMesh().time();

    // Update indices of times and weights
    bool timesChanged
    (
        interpolatorPtr_->valueWeights
        (
            t.timeOutputValue(),
            currentIndices_,
            currentWeights_
        )
    );

    const wordList currentTimeNames
    (
        UIndirectList<word>(timeNames_, currentIndices_)
    );


    // Load if necessary fields for this interpolation
    if (timesChanged)
    {
        objectRegistry& fieldsCache = const_cast<objectRegistry&>
        (
            pMesh.thisDb().subRegistry("fieldsCache", true)
        );

        // Save old times so we now which ones have been loaded and need
        // 'correctBoundaryConditions'. Bit messy.
        HashSet<word> oldTimes(fieldsCache.toc());

        ReadFields<pointVectorField>
        (
            fieldName_,
            pMesh,
            currentTimeNames
        );

        forAllConstIter(objectRegistry, fieldsCache, fieldsCacheIter)
        {
            if (!oldTimes.found(fieldsCacheIter.key()))
            {
                // Newly loaded fields. Make sure the internal
                // values are consistent with the boundary conditions.
                // This is quite often not the case since these
                // fields typically are constructed 'by hand'

                const objectRegistry& timeCache = dynamic_cast
                <
                    const objectRegistry&
                >(*fieldsCacheIter());

                timeCache.lookupObjectRef<pointVectorField>(fieldName_)
                    .correctBoundaryConditions();
            }
        }
    }


    // Interpolate the point field
    return
    (
        uniformInterpolate<pointVectorField>
        (
            IOobject
            (
                word("uniformInterpolate(") + fieldName_ + ')',
                pMesh.time().timeName(),
                pMesh.thisDb(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            fieldName_,
            currentTimeNames,
            currentWeights_
        )
    );
}


void Foam::dynamicMeshPointInterpolator::write(Ostream& os) const
{
    writeEntry(os, "field", fieldName_);
    writeEntry(os, "interpolationScheme", interpolationScheme_);
}


// ************************************************************************* //
