/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

Application
    testDataEntry

Description
    Tests DataEntry

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "DataEntry.H"
#include "IOdictionary.H"
#include "linearInterpolationWeights.H"
#include "splineInterpolationWeights.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"

{
    scalarField samples(4);
    samples[0] = 0;
    samples[1] = 1;
    samples[2] = 2;
    samples[3] = 3;
    scalarField values(4);
    values = 1.0;
    //values[0] = 0.0;
    //values[1] = 1.0;

    //linearInterpolationWeights interpolator
    splineInterpolationWeights interpolator
    (
        samples,
        interpolationWeights::WARN
    );
    labelList indices;
    scalarField weights;

    interpolator.integrationWeights(1.1, 1.2, indices, weights);
    Pout<< "indices:" << indices << endl;
    Pout<< "weights:" << weights << endl;

    scalar baseSum = interpolator.weightedSum
    (
        weights,
        UIndirectList<scalar>(values, indices)
    );
    Pout<< "baseSum=" << baseSum << nl << nl << endl;


//    interpolator.integrationWeights(-0.01, 0, indices, weights);
//    scalar partialSum = interpolator.weightedSum
//    (
//        weights,
//        UIndirectList<scalar>(values, indices)
//    );
//    Pout<< "partialSum=" << partialSum << nl << nl << endl;
//
//
//    interpolator.integrationWeights(-0.01, 1, indices, weights);
//    //Pout<< "samples:" << samples << endl;
//    //Pout<< "indices:" << indices << endl;
//    //Pout<< "weights:" << weights << endl;
//    scalar sum = interpolator.weightedSum
//    (
//        weights,
//        UIndirectList<scalar>(values, indices)
//    );
//    Pout<< "integrand=" << sum << nl << nl << endl;


    return 1;
}

    IOdictionary dataEntryProperties
    (
        IOobject
        (
            "dataEntryProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    autoPtr<DataEntry<scalar> > dataEntry
    (
        DataEntry<scalar>::New
        (
            "dataEntry",
            dataEntryProperties
        )
    );

    scalar x0 = readScalar(dataEntryProperties.lookup("x0"));
    scalar x1 = readScalar(dataEntryProperties.lookup("x1"));

    Info<< "Data entry type: " << dataEntry().type() << nl << endl;

    Info<< "Inputs" << nl
        << "    x0 = " << x0 << nl
        << "    x1 = " << x1 << nl
        << endl;

    Info<< "Interpolation" << nl
        << "    f(x0) = " << dataEntry().value(x0) << nl
        << "    f(x1) = " << dataEntry().value(x1) << nl
        << endl;

    Info<< "Integration" << nl
        << "    int(f(x)) lim(x0->x1) = " << dataEntry().integrate(x0, x1) << nl
        << endl;

    return 0;
}


// ************************************************************************* //
