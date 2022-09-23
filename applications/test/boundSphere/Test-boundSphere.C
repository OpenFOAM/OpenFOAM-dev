/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022 OpenFOAM Foundation
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

#include "argList.H"
#include "boundSphere.H"
#include "clock.H"
#include "cpuTime.H"
#include "OBJstream.H"
#include "triFace.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::validArgs.append("number of points");
    argList::addOption
    (
        "nTests",
        "value",
        "number of tests - default is 1"
    );
    argList::addOption
    (
        "seed",
        "value",
        "random number generator seed - default is clock time"
    );

    Foam::argList args(argc, argv);
    Info<< nl;

    const label nPs = args.argRead<label>(1);
    const label nTests = args.optionLookupOrDefault<label>("nTests", 1);

    // Initialise the random number generator
    label seed;
    if (args.optionFound("seed"))
    {
        seed = args.optionRead<label>("seed");
    }
    else
    {
        seed = clock::getTime();
        Info<< "Seeding random number generator with value " << seed
            << nl << endl;
    }

    // Do tests
    List<point> ps;
    Tuple2<point, scalar> sphere;
    cpuTime time;
    for (label testi = 0; testi < nTests; ++ testi)
    {
        Random rndGen(seed + testi);

        // Get random points
        ps.resize(nPs);
        forAll(ps, i)
        {
            ps[i] = point::uniform(1);

            while (magSqr(ps[i]) > 1)
            {
                ps[i] = 2*(rndGen.sample01<point>() - point::uniform(0.5));
            }
        }

        // Determine the bounds
        sphere = boundSphere(ps, rndGen);

        // Check
        forAll(ps, i)
        {
            if
            (
                magSqr(ps[i] - sphere.first())
              > (1 + rootSmall)*sqr(sphere.second())
            )
            {
                FatalErrorInFunction
                    << "Sphere does not bound point " << i << '/' << nPs << " ("
                    << magSqr(ps[i] - sphere.first())/sqr(sphere.second())
                    << ")" << exit(FatalError);
            }
        }

        if (nTests != 1)
        {
            Info<< "\rSuccessfully computed " << testi + 1 << " bound spheres";
        }
    }

    // Report
    if (nTests == 1)
    {
        Info<< "Bounding sphere with centre = " << sphere.first()
            << ", radius = " << sphere.second() << ", calculated in "
            << time.cpuTimeIncrement() << " s" << nl << endl;

        // Write the points
        OBJstream osPs(args[0] + "_points.obj");
        Info<< "Writing points to " << osPs.name() << endl;
        forAll(ps, i)
        {
            osPs.write(ps[i]);
        }

        // Write the sphere
        OBJstream osSphere(args[0] + "_sphere.obj");
        Info<< "Writing sphere to " << osSphere.name() << endl;

        using namespace constant::mathematical;

        static const label nRings = 20;
        static const label nRing1Points = 4;

        DynamicList<point> spherePs;
        DynamicList<face> sphereTris;
        auto appendP = [&](const scalar x, const scalar y, const scalar z)
        {
            spherePs.append(sphere.first() + sphere.second()*point(x, y, z));
            spherePs.append(sphere.first() + sphere.second()*point(x, y, -z));
        };
        auto appendTri = [&](const label a, const label b, const label c)
        {
            sphereTris.append(triFace(2*a, 2*b, 2*c));
            sphereTris.append(triFace(2*a + 1, 2*b + 1, 2*c + 1));
        };

        appendP(0, 0, 1);

        label nPoints0 = 0, nPoints = 1;

        for (label ringi = 1; ringi < nRings; ++ ringi)
        {
            const scalar rXY = Foam::sin(pi/2*ringi/(nRings - 1));

            const label ringPN0 = max((ringi - 1)*nRing1Points, 1);
            const label ringPN = ringi*nRing1Points;

            label ringPi0 = 0, ringPi = 0;

            for (label i = 0; i < nRing1Points; ++ i)
            {
                const scalar theta0 = 2*pi*i/nRing1Points;
                const scalar theta1 = 2*pi*(i + 1)/nRing1Points;

                for (label j = 0; j < ringi; ++ j)
                {
                    const scalar f = scalar(j)/ringi;
                    const scalar theta = (1 - f)*theta0 + f*theta1;

                    const scalar x = rXY*Foam::cos(theta);
                    const scalar y = rXY*Foam::sin(theta);
                    const scalar z = Foam::sqrt(max(1 - sqr(x) - sqr(y), 0));

                    appendP(x, y, z);

                    appendTri
                    (
                        nPoints0 + (ringPi0 % ringPN0),
                        nPoints + ringPi,
                        nPoints + ((ringPi + 1) % ringPN)
                    );

                    if (ringi > 1 && j < ringi - 1)
                    {
                        appendTri
                        (
                            nPoints0 + (ringPi0 % ringPN0),
                            nPoints + ((ringPi + 1) % ringPN),
                            nPoints0 + ((ringPi0 + 1) % ringPN0)
                        );
                    }

                    if (j < ringi - 1) ringPi0 ++;
                    ringPi ++;
                }
            }

            nPoints0 = nPoints;
            nPoints += ringPi;
        }

        osSphere.write(faceList(sphereTris), pointField(spherePs), false);
    }
    else
    {
        Info<< nl << nl << "Bounding spheres calculated in "
            << time.cpuTimeIncrement() << " s" << endl;
    }

    Info<< nl << "End" << nl << endl;

    return 0;
}


// ************************************************************************* //

