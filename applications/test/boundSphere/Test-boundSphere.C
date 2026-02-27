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

#include "argList.H"
#include "boundSphere.H"
#include "clock.H"
#include "OBJstream.H"
#include "triFace.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void writePoints(const word& name, const List<point>& ps)
{
    OBJstream osPs(name + "_points.obj");
    Sout<< "Writing points to " << osPs.name() << nl;
    forAll(ps, i)
    {
        osPs.write(ps[i]);
    }
}

void writeSphere(const word& name, const boundSphere& sphere)
{
    OBJstream osSphere(name + "_sphere.obj");
    Info<< "Writing sphere to " << osSphere.name() << nl << endl;

    using namespace constant::mathematical;

    static const label nRings = 20;
    static const label nRing1Points = 4;

    DynamicList<point> spherePs;
    DynamicList<face> sphereTris;
    auto appendP = [&](const scalar x, const scalar y, const scalar z)
    {
        spherePs.append(sphere.c() + sphere.r()*point(x, y, z));
        spherePs.append(sphere.c() + sphere.r()*point(x, y, -z));
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


void checkSphere
(
    const word& name,
    const List<point>& ps,
    const boundSphere& sphere,
    const FixedList<label, 4>& boundaryPis
)
{
    forAll(ps, pi)
    {
        if (findIndex(boundaryPis, pi) != -1) continue;

        if (!sphere.contains(ps[pi]))
        {
            const scalar rSqr = sphere.rSqr();
            const scalar prSqr = magSqr(ps[pi] - sphere.c());

            FatalErrorInFunction
                << name << " sphere does not bound point "
                << pi << '/' << ps.size() << " (" << prSqr/rSqr << ")"
                << exit(FatalError);
        }
    }
}


List<point> randomPoints(const label nPs, randomGenerator& rndGen)
{
    List<point> ps(nPs);

    forAll(ps, i)
    {
        ps[i] = point::uniform(1);

        while (magSqr(ps[i]) > 1)
        {
            ps[i] = 2*(rndGen.sample01<point>() - point::uniform(0.5));
        }
    }

    return ps;
}


void doLocalTest
(
    const word& name,
    const label nPs,
    randomGenerator& rndGen,
    const bool bruteForce
)
{
    const List<point> ps(randomPoints(nPs, rndGen));

    // Efficiently calculate the bound sphere
    FixedList<label, 4> welzlBoundaryPis;
    const boundSphere welzlSphere =
        boundSphere::local(ps, rndGen, welzlBoundaryPis);

    // Check that the sphere contains all the points
    checkSphere("Welzl", ps, welzlSphere, welzlBoundaryPis);

    if (bruteForce)
    {
        // Inefficiently calculate the bounding sphere
        FixedList<label, 4> bruteForceBoundaryPis;
        const boundSphere bruteForceSphere =
            boundSphere::bruteForce(ps, bruteForceBoundaryPis);

        // Check that the sphere contains all the points
        checkSphere("Brute-force", ps, bruteForceSphere, bruteForceBoundaryPis);

        // Check that the boundary points are the same
        if (bruteForceBoundaryPis != welzlBoundaryPis)
        {
            FatalError
                << "Welzl's algorithm boundary point set "
                << welzlBoundaryPis << " differs from the brute-force "
                << "boundary point set " << bruteForceBoundaryPis
                << exit(FatalError);
        }
    }

    if (name != word::null)
    {
        writePoints(name, ps);
        writeSphere(name, welzlSphere);
    }
}


void doGlobalTest
(
    const word& name,
    const label nPs,
    randomGenerator& rndGen,
    const bool bruteForce
)
{
    const List<point> ps(randomPoints(nPs, rndGen));

    // Efficiently calculate the bound sphere
    FixedList<RemoteData<point>, 4> welzlBoundaryPs;
    const boundSphere welzlSphere =
        boundSphere::global(ps, rndGen, welzlBoundaryPs, true);

    // Check that the sphere contains all the points
    {
        FixedList<label, 4> welzlBoundaryPis({-1, -1, -1, -1});
        for (label i = 0, j = 0; i < 4 && welzlBoundaryPs[i].proci != -1; ++ i)
        {
            if (welzlBoundaryPs[i].proci == Pstream::myProcNo())
            {
                welzlBoundaryPis[j ++] = welzlBoundaryPs[i].elementi;
            }
        }

        checkSphere("Welzl", ps, welzlSphere, welzlBoundaryPis);
    }

    FixedList<remote, 4> welzlBoundaryPis;
    forAll(welzlBoundaryPis, i) welzlBoundaryPis[i] = welzlBoundaryPs[i];

    if (bruteForce)
    {
        // Gather all the points
        List<List<point>> procPs(Pstream::nProcs());
        procPs[Pstream::myProcNo()] = ps;
        Pstream::gatherList(procPs);

        if (Pstream::master())
        {
            // Flatten the gathered points into a single list and maintain a
            // set of offsets into where each processor's "block" of points
            // start, so that we can work the local index of each point
            labelList procPsOffsets(Pstream::nProcs() + 1, 0);
            forAll(procPs, proci)
            {
                procPsOffsets[proci + 1] =
                    procPsOffsets[proci] + procPs[proci].size();
            }
            List<point> allPs(procPsOffsets.last());
            forAll(procPs, proci)
            {
                SubList<point>
                (
                    allPs,
                    procPsOffsets[proci + 1] - procPsOffsets[proci],
                    procPsOffsets[proci]
                ) = procPs[proci];
            }

            // Inefficiently calculate the bounding sphere
            FixedList<label, 4> bruteForceBoundaryAllPis;
            const boundSphere bruteForceSphere =
                boundSphere::bruteForce(allPs, bruteForceBoundaryAllPis);

            // Check that the sphere contains all the points
            checkSphere
            (
                "Brute-force",
                allPs,
                bruteForceSphere,
                bruteForceBoundaryAllPis
            );

            // Convert the boundary all-point indices into remote-point indices
            FixedList<remote, 4> bruteForceBoundaryPis;
            for
            (
                label i = 0, proci = 0;
                i < 4 && bruteForceBoundaryAllPis[i] != -1;
                ++ i
            )
            {
                while (bruteForceBoundaryAllPis[i] >= procPsOffsets[proci + 1])
                    proci ++;

                bruteForceBoundaryPis[i].proci = proci;
                bruteForceBoundaryPis[i].elementi =
                    bruteForceBoundaryAllPis[i] - procPsOffsets[proci];
            }

            // Check that the boundary points are the same
            if (bruteForceBoundaryPis != welzlBoundaryPis)
            {
                FatalError
                    << "Welzl's algorithm boundary point set "
                    << welzlBoundaryPis << " differs from the brute-force "
                    << "boundary point set " << bruteForceBoundaryPis
                    << exit(FatalError);
            }
        }
    }

    if (name != word::null)
    {
        writePoints(name + "_proc" + Foam::name(Pstream::myProcNo()), ps);
        if (Pstream::master())
        {
            writeSphere(name, welzlSphere);
        }
    }
}


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
    argList::addBoolOption
    (
        "bruteForce",
        "also check against a bound sphere obtained via a brute-force search"
        " - default true if the number of points is less than 64"
    );

    Foam::argList args(argc, argv);

    const label nPs = args.argRead<label>(1);
    const label nTests = args.optionLookupOrDefault<label>("nTests", 1);
    const bool bruteForce = args.optionFound("bruteForce") || nPs < 64;

    Info<< endl;

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
    randomGenerator rndGen(seed);

    // Run tests
    for (label testi = 0; testi < nTests; ++ testi)
    {
        if (Pstream::parRun())
        {
            doGlobalTest
            (
                nTests == 1 ? argv[0] : word::null,
                nPs/Pstream::nProcs(),
                rndGen,
                bruteForce
            );
        }
        else
        {
            doLocalTest
            (
                nTests == 1 ? argv[0] : word::null,
                nPs,
                rndGen,
                bruteForce
            );
        }

        if (nTests != 1)
        {
            Info<< "\rSuccessfully computed " << testi + 1
                << " bound spheres" << flush;
        }
    }

    if (nTests != 1)
    {
        Info<< nl << endl;
    }

    Info<< "End" << nl << endl;

    return 0;
}


// ************************************************************************* //

