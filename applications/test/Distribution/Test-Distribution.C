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

Application
    Test-Distribution

Description
    Test the Distribution class

    Plot normal distribution test in gnuplot using:

    \verbatim
    normalDistribution(mean, sigma, x) = \
        sqrt(1.0/(2.0*pi*sigma**2))*exp(-(x - mean)**2.0/(2.0*sigma**2))

    plot normalDistribution(8.5, 2.5, x), "Distribution_scalar_test_x" w p
    \endverbatim

\*---------------------------------------------------------------------------*/

#include "vector.H"
#include "labelVector.H"
#include "tensor.H"
#include "Distribution.H"
#include "Random.H"
#include "dimensionedTypes.H"
#include "argList.H"
#include "PstreamReduceOps.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace Foam;

int main(int argc, char *argv[])
{
    #include "setRootCase.H"

    Random R(918273);

    {
        // scalar
        label randomDistributionTestSize = 50000000;

        Distribution<scalar> dS(scalar(5e-2));

        Info<< nl << "Distribution<scalar>" << nl
            << "Sampling "
            << randomDistributionTestSize
            << " times from a standard normal distribution."
            << endl;

        for (label i = 0; i < randomDistributionTestSize; i++)
        {
            dS.add(2.5*R.scalarNormal() + 8.5);
        }

        Info<< "Mean " << dS.mean() << nl
            << "Median " << dS.median()
            << endl;

        dS.write("Distribution_scalar_test_1");

        Distribution<scalar> dS2(scalar(1e-2));

        Info<< nl << "Distribution<scalar>" << nl
            << "Sampling "
            << randomDistributionTestSize
            << " times from a standard normal distribution."
            << endl;

        for (label i = 0; i < randomDistributionTestSize; i++)
        {
            dS2.add(1.5*R.scalarNormal() -6.0);
        }

        Info<< "Mean " << dS2.mean() << nl
            << "Median " << dS2.median()
            << endl;

        dS2.write("Distribution_scalar_test_2");

        Info<< nl << "Adding previous two Distribution<scalar>" << endl;

        dS = dS + dS2;

        dS.write("Distribution_scalar_test_1+2");
    }

    if (Pstream::parRun())
    {
        // scalar in parallel
        label randomDistributionTestSize = 100000000;

        Distribution<scalar> dS(scalar(1e-1));

        Pout<< "Distribution<scalar>" << nl
            << "Sampling "
            << randomDistributionTestSize
            << " times from uniform distribution."
            << endl;

        for (label i = 0; i < randomDistributionTestSize; i++)
        {
            dS.add(R.scalar01() + 10*Pstream::myProcNo());
        }

        Pout<< "Mean " << dS.mean() << nl
            << "Median " << dS.median()
            << endl;

        reduce(dS, sumOp<Distribution<scalar>>());

        if (Pstream::master())
        {
            Info<< "Reducing parallel Distribution<scalar>" << nl
                << "Mean " << dS.mean() << nl
                << "Median " << dS.median()
                << endl;

            dS.write("Distribution_scalar_test_parallel_reduced");
        }
    }

    {
        // vector
        Distribution<vector> dV(vector(0.1, 0.05, 0.15));

        label randomDistributionTestSize = 1000000;

        Info<< nl << "Distribution<vector>" << nl
            << "Sampling "
            << randomDistributionTestSize
            << " times from uniform and a standard normal distribution."
            << endl;

        for (label i = 0; i < randomDistributionTestSize; i++)
        {
            dV.add(R.sample01<vector>());

            // Adding separate standard normal components with component
            // weights

            dV.add
            (
                vector
                (
                    R.scalarNormal()*3.0 + 1.5,
                    R.scalarNormal()*0.25 + 4.0,
                    R.scalarNormal()*3.0 - 1.5
                ),
                vector(1.0, 2.0, 5.0)
            );
        }

        Info<< "Mean " << dV.mean() << nl
            << "Median " << dV.median()
            << endl;

        dV.write("Distribution_vector_test");
    }

    // {
    //     // labelVector
    //     Distribution<labelVector> dLV(labelVector::one*10);

    //     label randomDistributionTestSize = 2000000;

    //     Info<< nl << "Distribution<labelVector>" << nl
    //         << "Sampling "
    //         << randomDistributionTestSize
    //         << " times from uniform distribution."
    //         << endl;

    //     for (label i = 0; i < randomDistributionTestSize; i++)
    //     {
    //         dLV.add
    //         (
    //             labelVector
    //             (
    //                 R.sampleAB<label>(-1000, 1001),
    //                 R.sampleAB<label>(-5000, 5001),
    //                 R.sampleAB<label>(-2000, 7001)
    //             )
    //         );
    //     }

    //     Info<< "Mean " << dLV.mean() << nl
    //         << "Median " << dLV.median()
    //         << endl;

    //     dLV.write("Distribution_labelVector_test");
    // }

    {
        // tensor
        Distribution<tensor> dT(tensor::one*1e-2);

        label randomDistributionTestSize = 2000000;

        Info<< nl << "Distribution<tensor>" << nl
            << "Sampling "
            << randomDistributionTestSize
            << " times from uniform distribution."
            << endl;

        for (label i = 0; i < randomDistributionTestSize; i++)
        {
            dT.add(R.sample01<tensor>());
        }

        Info<< "Mean " << dT.mean() << nl
            << "Median " << dT.median()
            << endl;

        dT.write("Distribution_tensor_test");
    }

    {
        // symmTensor
        Distribution<symmTensor> dSyT(symmTensor::one*1e-2);

        label randomDistributionTestSize = 2000000;

        Info<< nl << "Distribution<symmTensor>" << nl
            << "Sampling "
            << randomDistributionTestSize
            << " times from uniform distribution."
            << endl;

        for (label i = 0; i < randomDistributionTestSize; i++)
        {
            dSyT.add(R.sample01<symmTensor>());
        }

        Info<< "Mean " << dSyT.mean() << nl
            << "Median " << dSyT.median()
            << endl;

        dSyT.write("Distribution_symmTensor_test");
    }

    {
        // sphericalTensor
        Distribution<sphericalTensor> dSpT(sphericalTensor::one*1e-2);

        label randomDistributionTestSize = 50000000;

        Info<< nl << "Distribution<sphericalTensor>" << nl
            << "Sampling "
            << randomDistributionTestSize
            << " times from uniform distribution."
            << endl;

        for (label i = 0; i < randomDistributionTestSize; i++)
        {
            dSpT.add(R.sample01<sphericalTensor>());
        }

        Info<< "Mean " << dSpT.mean() << nl
            << "Median " << dSpT.median()
            << endl;

        dSpT.write("Distribution_sphericalTensor_test");
    }

    Info<< nl << "End" << nl << endl;

    return 0;
}


// ************************************************************************* //
