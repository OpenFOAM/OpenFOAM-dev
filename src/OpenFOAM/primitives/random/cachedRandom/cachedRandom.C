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

\*---------------------------------------------------------------------------*/

#include "cachedRandom.H"
#include "OSspecific.H"
#include "PstreamReduceOps.H"

// * * * * * * * * * * * * * private Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::cachedRandom::scalar01()
{
    if (sampleI_ < 0)
    {
        return osRandomDouble();
    }

    if (sampleI_ == samples_.size() - 1)
    {
        scalar s = samples_[sampleI_];
        sampleI_ = 0;
        return s;
    }
    else
    {
        scalar s = samples_[sampleI_];
        sampleI_++;
        return s;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cachedRandom::cachedRandom(const label seed, const label count)
:
    seed_(1),
    samples_(0),
    sampleI_(-1)
{
    if (seed > 1)
    {
        seed_ = seed;
    }

    // Samples will be cached if count > 0
    if (count > 0)
    {
        samples_.setSize(count);
        sampleI_ = 0;
    }

    // Initialise samples
    osRandomSeed(seed_);
    forAll(samples_, i)
    {
        samples_[i] = osRandomDouble();
    }
}


Foam::cachedRandom::cachedRandom(const cachedRandom& cr, const bool reset)
:
    seed_(cr.seed_),
    samples_(cr.samples_),
    sampleI_(cr.sampleI_)
{
    if (sampleI_ == -1)
    {
        WarningIn
        (
            "Foam::cachedRandom::cachedRandom(const cachedRandom& cr)"
        )   << "Copy constructor called, but samples not being cached. "
            << "This may lead to non-repeatable behaviour" << endl;

        osRandomSeed(seed_);
    }

    if (reset && samples_.size())
    {
        sampleI_ = 0;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::cachedRandom::~cachedRandom()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<>
Foam::label Foam::cachedRandom::sample01()
{
    return round(scalar01());
}


template<>
Foam::scalar Foam::cachedRandom::sample01()
{
    return scalar01();
}


template<>
Foam::label Foam::cachedRandom::position(const label& start, const label& end)
{
    return start + round(scalar01()*(end - start));
}


template<>
Foam::scalar Foam::cachedRandom::position
(
    const scalar& start,
    const scalar& end
)
{
    return start + scalar01()*(end - start);
}


template<>
Foam::label Foam::cachedRandom::globalSample01()
{
    scalar value = -GREAT;

    if (Pstream::master())
    {
        value = scalar01();
    }

    reduce(value, maxOp<scalar>());

    return round(value);
}


template<>
Foam::scalar Foam::cachedRandom::globalSample01()
{
    scalar value = -GREAT;

    if (Pstream::master())
    {
        value = scalar01();
    }

    reduce(value, maxOp<scalar>());

    return value;
}


template<>
Foam::label Foam::cachedRandom::globalPosition
(
    const label& start,
    const label& end
)
{
    label value = labelMin;

    if (Pstream::master())
    {
        value = round(scalar01()*(end - start));
    }

    reduce(value, maxOp<label>());

    return start + value;
}


template<>
Foam::scalar Foam::cachedRandom::globalPosition
(
    const scalar& start,
    const scalar& end
)
{
    scalar value = -GREAT;

    if (Pstream::master())
    {
        value = scalar01()*(end - start);
    }

    reduce(value, maxOp<scalar>());

    return start + value;
}


void Foam::cachedRandom::operator=(const cachedRandom& cr)
{
    seed_ = cr.seed_;
    samples_ = cr.samples_;
    sampleI_ = cr.sampleI_;
}


// ************************************************************************* //
