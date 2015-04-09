/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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

#include "Random.H"
#include "OSspecific.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#if INT_MAX    != 2147483647
#    error "INT_MAX    != 2147483647"
#    error "The random number generator may not work!"
#endif


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::Random::Random(const label seed)
{
    if (seed > 1)
    {
        Seed = seed;
    }
    else
    {
        Seed = 1;
    }

    osRandomSeed(Seed);
}


int Foam::Random::bit()
{
    if (osRandomInteger() > INT_MAX/2)
    {
        return 1;
    }
    else
    {
        return 0;
    }
}


Foam::scalar Foam::Random::scalar01()
{
    return osRandomDouble();
}


Foam::vector Foam::Random::vector01()
{
    vector rndVec;
    for (direction cmpt=0; cmpt<vector::nComponents; cmpt++)
    {
        rndVec.component(cmpt) = scalar01();
    }

    return rndVec;
}


Foam::sphericalTensor Foam::Random::sphericalTensor01()
{
    sphericalTensor rndTen;
    rndTen.ii() = scalar01();

    return rndTen;
}


Foam::symmTensor Foam::Random::symmTensor01()
{
    symmTensor rndTen;
    for (direction cmpt=0; cmpt<symmTensor::nComponents; cmpt++)
    {
        rndTen.component(cmpt) = scalar01();
    }

    return rndTen;
}


Foam::tensor Foam::Random::tensor01()
{
    tensor rndTen;
    for (direction cmpt=0; cmpt<tensor::nComponents; cmpt++)
    {
        rndTen.component(cmpt) = scalar01();
    }

    return rndTen;
}


Foam::label Foam::Random::integer(const label lower, const label upper)
{
    return lower + (osRandomInteger() % (upper+1-lower));
}


Foam::vector Foam::Random::position(const vector& start, const vector& end)
{
    vector rndVec(start);

    for (direction cmpt=0; cmpt<vector::nComponents; cmpt++)
    {
        rndVec.component(cmpt) +=
            scalar01()*(end.component(cmpt) - start.component(cmpt));
    }

    return rndVec;
}


void Foam::Random::randomise(scalar& s)
{
     s = scalar01();
}


void Foam::Random::randomise(vector& v)
{
    v = vector01();
}


void Foam::Random::randomise(sphericalTensor& st)
{
    st = sphericalTensor01();
}


void Foam::Random::randomise(symmTensor& st)
{
    st = symmTensor01();
}


void Foam::Random::randomise(tensor& t)
{
    t = tensor01();
}


Foam::scalar Foam::Random::GaussNormal()
{
    static int iset = 0;
    static scalar gset;
    scalar fac, rsq, v1, v2;

    if (iset == 0)
    {
        do
        {
            v1 = 2.0*scalar01() - 1.0;
            v2 = 2.0*scalar01() - 1.0;
            rsq = sqr(v1) + sqr(v2);
        } while (rsq >= 1.0 || rsq == 0.0);

        fac = sqrt(-2.0*log(rsq)/rsq);
        gset = v1*fac;
        iset = 1;

        return v2*fac;
    }
    else
    {
        iset = 0;

        return gset;
    }
}


// ************************************************************************* //
