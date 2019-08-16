/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019 OpenFOAM Foundation
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
    Test-Distribution2

Description
    Test the general distributionModel.

\*---------------------------------------------------------------------------*/

#include "Distribution.H"
#include "Random.H"
#include "dimensionedTypes.H"
#include "argList.H"
#include "distributionModel.H"
#include "IFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace Foam;

int main(int argc, char *argv[])
{
    #include "setRootCase.H"

    Random R(918273);

    dictionary dict(IFstream("testDict")());

    Info << nl << "Testing general distribution" << endl;
    Info << nl << "Continuous probability density function:" << endl;

    autoPtr<distributionModel> dist1
    (
        distributionModel::New
        (
            dict.subDict("densityFunction"),
            R
        )
    );

    label randomDistributionTestSize = 50000000;
    Distribution<scalar> dS(scalar(1e-6));

    Info<< nl
        << "Sampling " << randomDistributionTestSize << " times." << endl;

    for (label i = 0; i < randomDistributionTestSize; i++)
    {
        dS.add(dist1->sample());
    }

    Info<< "Produced mean " << dS.mean() << endl;
    dS.write("densityTest");
    dS.clear();

    Info << nl << "Discrete probability density function:" << endl;
    dist1.clear();

    dist1 =
        distributionModel::New
        (
            dict.subDict("cumulativeFunction"),
            R
        );

    Info<< nl
        << "Sampling " << randomDistributionTestSize << " times." << endl;

    for (label i = 0; i < randomDistributionTestSize; i++)
    {
        dS.add(dist1->sample());
    }

    Info<< "Produced mean " << dS.mean() << endl;
    dS.write("cumulativeTest");

    Info<< nl << "End" << nl << endl;

    return 0;
}


// ************************************************************************* //
