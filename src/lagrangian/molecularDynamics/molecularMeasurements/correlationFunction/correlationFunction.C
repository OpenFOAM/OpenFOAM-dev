/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "correlationFunction.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class Type>
const char* const
    Foam::correlationFunction<Type>::typeName("correlationFunction");


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::correlationFunction<Type>::setTimesAndSizes
(
    const label tZeroBufferSize
)
{
    sampleSteps_  = ceil(sampleInterval_/mesh_.time().deltaTValue());

    sampleInterval_ = sampleSteps_*mesh_.time().deltaTValue();

    label bufferLength(ceil(duration_/sampleInterval_));

    duration_ = bufferLength*sampleInterval_;

    label bufferingInterval(ceil(averagingInterval_/sampleInterval_));

    averagingInterval_ = bufferingInterval*sampleInterval_;

    label nBuffers(ceil(duration_/averagingInterval_));

    this->setSizes
    (
        nBuffers,
        bufferLength,
        bufferingInterval
    );

    tZeroBuffers_ =
        Field<Field<Type>>
        (
            nBuffers,
            Field<Type>
            (
                tZeroBufferSize,
                Zero
            )
        );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::correlationFunction<Type>::correlationFunction
(
    const polyMesh& mesh,
    const dictionary& cfDict,
    const label tZeroBufferSize
)
:
    bufferedAccumulator<scalar>(),
    mesh_(mesh)
{
    duration_ = readScalar(cfDict.lookup("duration"));

    sampleInterval_ = readScalar(cfDict.lookup("sampleInterval"));

    averagingInterval_ = readScalar(cfDict.lookup("averagingInterval"));

    setTimesAndSizes(tZeroBufferSize);
}


template<class Type>
Foam::correlationFunction<Type>::correlationFunction
(
    const polyMesh& mesh,
    const label tZeroBufferSize,
    const scalar duration,
    const scalar sampleInterval,
    const scalar averagingInterval
)
:
    bufferedAccumulator<scalar>(),
    mesh_(mesh),
    duration_(duration),
    sampleInterval_(sampleInterval),
    averagingInterval_(averagingInterval)
{
    setTimesAndSizes(tZeroBufferSize);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::correlationFunction<Type>::~correlationFunction()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::correlationFunction<Type>::calculateCorrelationFunction
(
    const Field<Type>& currentValues
)
{
    if (measurandFieldSize() != currentValues.size())
    {
        FatalErrorInFunction
            << "Trying to supply a Field of length"
            << currentValues.size()
            << " to calculate the correlation function. "
            << "Expecting a Field of length "
            << measurandFieldSize() << nl
            << abort(FatalError);
    }

    List<scalar> cFSums(nBuffers(),0.0);

    forAll(tZeroBuffers_, tZB)
    {
        scalar& cFSum = cFSums[tZB];

        const Field<Type>& tZeroBuffer = tZeroBuffers_[tZB];

        forAll(currentValues, cV)
        {
            const Type& tZeroBufferValue = tZeroBuffer[cV];

            const Type& currentValue = currentValues[cV];

            forAll(currentValue, component)
            {
                cFSum +=
                (
                    tZeroBufferValue[component]*currentValue[component]
                );
            }
        }

        cFSum /= (measurandFieldSize()*currentValues[0].size());
    }

    label bufferToRefill = addToBuffers(cFSums);

    if (bufferToRefill != -1)
    {
        tZeroBuffers_[bufferToRefill] = currentValues;
    }
}


template<class Type>
void Foam::correlationFunction<Type>::calculateCorrelationFunction
(
    const Type& currentValue
)
{
    if (measurandFieldSize() != 1)
    {
        FatalErrorInFunction
            << "Trying to supply a single value to calculate the correlation "
            << "function.  Expecting a Field of length "
            << measurandFieldSize()
            << abort(FatalError);
    }

    calculateCorrelationFunction(Field<Type>(1, currentValue));
}


template<class Type>
Foam::scalar Foam::correlationFunction<Type>::integral() const
{
    Field<scalar> averageCF(averaged());

    scalar cFIntegral = 0.0;

    for (label v = 0; v < averageCF.size() - 1; v++)
    {
        cFIntegral +=
            0.5
           *sampleInterval_
           *(averageCF[v+1] + averageCF[v]);
    }

    return cFIntegral;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    #include "correlationFunctionIO.C"

// ************************************************************************* //
