/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023-2024 OpenFOAM Foundation
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
    Test-distribution

Description
    Test distributions

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "clock.H"
#include "distribution.H"
#include "IFstream.H"
#include "OFstream.H"
#include "rawSetWriter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace Foam;

int main(int argc, char *argv[])
{
    argList::validArgs.append("dictionary");
    argList::validArgs.append("sampleQ");
    argList args(argc, argv);

    // Create the input dictionary
    const word dictName(args.argRead<word>(1));
    Info<< "Reading " << dictName << nl << endl;
    dictionary dict = IFstream(dictName)();

    // Get the number of size exponents to consider
    const label nQ = args.argRead<label>(2);

    // Discretisation of sample space. Make arguments?
    static const label nAnalytic = 1000;
    static const label nSampled = 100;
    static const label nSamples = 10000000;

    // Create a zero-size-exponent, zero-sampling-size-exponent distribution
    // from which to get X-coordinates
    dict.set("Q", 0);
    autoPtr<distribution> distribution00Ptr
    (
        distribution::New(unitAny, dict, 0, clock::getTime())
    );

    // Get the X-coordinates for the plots
    const scalarField xAnalytic(distribution00Ptr->x(nAnalytic));

    const scalar xSampled0 = distribution00Ptr->min();
    const scalar xSampled1 = distribution00Ptr->max();
    const scalar dxSampled = max(xSampled1 - xSampled0, rootVSmall);
    scalarField xSampled(nSampled);
    forAll(xSampled, i)
    {
        const scalar f = scalar(i)/(nSampled - 1);
        xSampled[i] = (1 - f)*xSampled0 + f*xSampled1;
    }

    // Declare the Y-names and Y-coordinates for the plots
    wordList yNames;
    #define DeclareY(Type, nullArg) \
        PtrList<Field<Type>> Type##YAnalytic, Type##YSampled;
    FOR_ALL_FIELD_TYPES(DeclareY);
    #undef DeclareY

    // Create plot file and output initial plot commands
    OFstream plot(dictName + ".gnuplot");
    plot<< "set terminal eps size " << 4*(nQ + 1) << "," << 3*(nQ + 1) << nl
        << "set output '" << dictName << ".eps'" << nl
        << "set multiplot layout " << nQ + 1 << "," << nQ + 1 << nl;

    // Loop distribution size exponents
    for (label Q = 0; Q <= nQ; ++ Q)
    {
        // Create a Q-size-exponent, zero-sampling-size-exponent distribution
        dict.set("Q", Q);
        autoPtr<distribution> distributionQ0Ptr
        (
            Q == 0
          ? distribution00Ptr->clone(0)
          : distribution::New(unitAny, dict, 0, clock::getTime(), false, false)
        );

        // Resize
        const label Q0i = yNames.size();
        yNames.append("Q=" + name(Q));
        #define AppendY(Type, nullArg) \
            Type##YAnalytic.resize(Type##YAnalytic.size() + 1); \
            Type##YSampled.resize(Type##YSampled.size() + 1);
        FOR_ALL_FIELD_TYPES(AppendY);
        #undef AppendY

        // Compute the analytic PDF
        scalarYAnalytic.set
        (
            Q0i,
            distributionQ0Ptr->PDF
            (
                distributionQ0Ptr->x(nAnalytic)
            ).ptr()
        );

        // Compute the sampled PDF
        scalarYSampled.set(Q0i, new scalarField(nSampled, 0));

        scalarField samples(distributionQ0Ptr->sample(nSamples));

        forAll(samples, samplei)
        {
            const scalar x = samples[samplei];
            const scalar f = (x - xSampled0)/dxSampled;
            const label i = min(floor(f*(nSampled - 1)), nSampled - 2);
            const scalar g = f*(nSampled - 1) - scalar(i);

            scalarYSampled[Q0i][i] += (1 - g);
            scalarYSampled[Q0i][i + 1] += g;
        }

        scalarYSampled[Q0i] /=
            sum(scalarYSampled[Q0i])*dxSampled/(nSampled - 1);
        scalarYSampled[Q0i].first() *= 2;
        scalarYSampled[Q0i].last() *= 2;

        // Loop sampling size exponents
        for (label sampleQ = 0; sampleQ <= nQ; ++ sampleQ)
        {
            // Create a Q-size-exponent, Q-sampling-size-exponent distribution
            autoPtr<distribution> distributionQSampleQPtr
            (
                distributionQ0Ptr->clone(sampleQ)
            );
            {
                OStringStream oss;
                distributionQSampleQPtr->write(oss, unitAny);
                writeEntry(oss, "sampleQ", sampleQ);
                writeEntry(oss, "min", distributionQSampleQPtr->min());
                writeEntry(oss, "mean", distributionQSampleQPtr->mean());
                writeEntry(oss, "max", distributionQSampleQPtr->max());
                Info<< dictionary(IStringStream(oss.str())());
            }

            // Resize
            const label QSampleQi = yNames.size();
            const label QSampleQiStar = yNames.size() + 1;;
            yNames.append("Q=" + name(Q) + ", sampleQ=" + name(sampleQ));
            yNames.append(yNames.last() + ", *x^{" + name(-sampleQ) + "}");
            #define AppendY(Type, nullArg) \
                Type##YAnalytic.resize(Type##YAnalytic.size() + 2); \
                Type##YSampled.resize(Type##YSampled.size() + 2);
            FOR_ALL_FIELD_TYPES(AppendY);

            // Compute the analytic PDF
            scalarYAnalytic.set
            (
                QSampleQi,
                distributionQSampleQPtr->PDF
                (
                    distributionQSampleQPtr->x(nAnalytic)
                ).ptr()
            );
            scalarYAnalytic.set
            (
                QSampleQiStar,
                scalarYAnalytic[Q0i]
            );

            // Compute the sampled PDF
            scalarYSampled.set(QSampleQi, new scalarField(nSampled, 0));
            scalarYSampled.set(QSampleQiStar, new scalarField(nSampled, 0));

            scalarField samples(distributionQSampleQPtr->sample(nSamples));

            forAll(samples, samplei)
            {
                const scalar x = samples[samplei];
                const scalar xPowNegSampleQ = Foam::pow(x, - sampleQ);
                const scalar f = (x - xSampled0)/dxSampled;
                const label i = min(floor(f*(nSampled - 1)), nSampled - 2);
                const scalar g = f*(nSampled - 1) - scalar(i);

                scalarYSampled[QSampleQi][i] += (1 - g);
                scalarYSampled[QSampleQi][i + 1] += g;

                scalarYSampled[QSampleQiStar][i] += (1 - g)*xPowNegSampleQ;
                scalarYSampled[QSampleQiStar][i + 1] += g*xPowNegSampleQ;
            }

            scalarYSampled[QSampleQi] /=
                sum(scalarYSampled[QSampleQi])*dxSampled/(nSampled - 1);
            scalarYSampled[QSampleQi].first() *= 2;
            scalarYSampled[QSampleQi].last() *= 2;

            scalarYSampled[QSampleQiStar] /=
                sum(scalarYSampled[QSampleQiStar])*dxSampled/(nSampled - 1);
            scalarYSampled[QSampleQiStar].first() *= 2;
            scalarYSampled[QSampleQiStar].last() *= 2;

            plot<< "plot \\" << nl
                << "    '" << dictName << ".analytic.xy' us 1:"
                << Q0i + 2 << " w l t 'Analytic "
                << yNames[Q0i] << "', \\" << nl
                << "    '" << dictName << ".sampled.xy' us 1:"
                << Q0i + 2 << " w l t 'Sampled "
                << yNames[Q0i] << "', \\" << nl
                << "    '" << dictName << ".analytic.xy' us 1:"
                << QSampleQi + 2 << " w l t 'Analytic "
                << yNames[QSampleQi] << "', \\" << nl
                << "    '" << dictName << ".sampled.xy' us 1:"
                << QSampleQi + 2 << " w l t 'Sampled "
                << yNames[QSampleQi] << "', \\" << nl
                << "    '" << dictName << ".sampled.xy' us 1:"
                << QSampleQiStar + 2 << " w l t 'Sampled "
                << yNames[QSampleQiStar] << "'" << nl;
        }
    }

    plot<< "unset multiplot" << nl;

    forAll(yNames, i)
    {
        yNames[i].replaceAll(" ", "");
    }

    IOstream::defaultPrecision(15);

    rawSetWriter(IOstream::ASCII, IOstream::UNCOMPRESSED).write
    (
        ".",
        dictName + ".analytic",
        coordSet(true, "x", xAnalytic),
        yNames
        #define TypeYAnalyticParameter(Type, nullArg) , Type##YAnalytic
        FOR_ALL_FIELD_TYPES(TypeYAnalyticParameter)
        #undef TypeYAnalyticParameter
    );

    rawSetWriter(IOstream::ASCII, IOstream::UNCOMPRESSED).write
    (
        ".",
        dictName + ".sampled",
        coordSet(true, "x", xSampled),
        yNames
        #define TypeYSampledParameter(Type, nullArg) , Type##YSampled
        FOR_ALL_FIELD_TYPES(TypeYSampledParameter)
        #undef TypeYSampledParameter
    );

    Info<< nl << "End" << nl << endl;

    return 0;
}


// ************************************************************************* //
