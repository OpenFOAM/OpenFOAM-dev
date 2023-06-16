/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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
    noise

Description
    Utility to perform noise analysis of pressure data using the noiseFFT
    library.

    Control settings are read from the $FOAM_CASE/system/noiseDict dictionary,
    or user-specified dictionary using the -dict option.  Pressure data is
    read using a Table Function1:

Usage
    \verbatim
    pRef        101325;
    N           65536;
    nw          100;
    f1          25;
    fU          10000;
    graphFormat raw;

    pressureData
    {
        file                "pressureData";
        nHeaderLine         1;          // number of header lines
        refColumn           0;          // reference column index
        componentColumns    (1);        // component column indices
        separator           " ";        // optional (defaults to ",")
        mergeSeparators     no;         // merge multiple separators
        outOfBounds         clamp;      // optional out-of-bounds handling
        interpolationScheme linear;     // optional interpolation scheme
    }
    \endverbatim

    where
    \table
        Property    | Description                   | Required  | Default value
        pRef        | Reference pressure            | no        | 0
        N           | Number of samples in sampling window | no | 65536
        nw          | Number of sampling windows    | no        | 100
        fl          | Lower frequency band          | no        | 25
        fU          | Upper frequency band          | no        | 10000
        graphFormat | Output graph format          | no        | raw
    \endtable

    Current graph outputs include:
    - FFT of the pressure data
    - narrow-band PFL (pressure-fluctuation level) spectrum
    - one-third-octave-band PFL spectrum
    - one-third-octave-band pressure spectrum

See also
    Table.H
    noiseFFT.H

\*---------------------------------------------------------------------------*/

#include "noiseFFT.H"
#include "argList.H"
#include "Time.H"
#include "Table.H"
#include "systemDict.H"
#include "setWriter.H"
#include "writeFile.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::scalar checkUniformTimeStep(const scalarField& t)
{
    // check that a uniform time step has been applied
    scalar deltaT = -1.0;
    if (t.size() > 1)
    {
        for (label i = 1; i < t.size(); i++)
        {
            const scalar dT = t[i] - t[i-1];
            if (deltaT < 0)
            {
                deltaT = dT;
            }

            if (mag(deltaT - dT) > rootSmall)
            {
                FatalErrorInFunction
                    << "Unable to process data with a variable time step"
                    << exit(FatalError);
            }
        }
    }
    else
    {
        FatalErrorInFunction
            << "Unable to create FFT with a single value"
            << exit(FatalError);
    }

    return deltaT;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noParallel();
    #include "addDictOption.H"
    #include "setRootCase.H"
    #include "createTime.H"

    const word dictName("noiseDict");

    IOdictionary dict(systemDict(dictName, args, runTime));

    // Reference pressure
    const scalar pRef = dict.lookupOrDefault("pRef", 0.0);

    // Number of samples in sampling window
    const label N = dict.lookupOrDefault("N", 65536);

    // Number of sampling windows
    const label nw = dict.lookupOrDefault("nw", 100);

    // Lower frequency of frequency band
    const scalar f1 = dict.lookupOrDefault("f1", 25.0);

    // Upper frequency of frequency band
    const scalar fU = dict.lookupOrDefault("fU", 10000.0);

    // Graph format
    const word graphFormat
    (
        dict.lookupOrDefault<word>("graphFormat", "raw")
    );

    Info<< "Reading data file" << endl;
    Function1s::Table<scalar> pData
    (
        "pressure",
        dict.subDict("pressureData")
    );

    // time history data
    const scalarField t(pData.x());

    // pressure data
    const scalarField p(pData.y());

    if (t.size() < N)
    {
        FatalErrorInFunction
            << "Block size N = " << N
            << " is larger than number of data = " << t.size()
            << exit(FatalError);
    }

    Info<< "    read " << t.size() << " values" << nl << endl;


    Info<< "Creating noise FFT" << endl;
    noiseFFT nfft(checkUniformTimeStep(t), p);

    nfft -= pRef;

    const fileName pDateFileName(dict.subDict("pressureData").lookup("file"));
    const fileName baseFileName(pDateFileName.lessExt());
    const fileName outputPath
    (
        runTime.path()
       /functionObjects::writeFile::outputPrefix
       /baseFileName
    );

    autoPtr<setWriter> writer(setWriter::New(graphFormat));

    const Pair<scalarField> Pf(nfft.RMSmeanPf(N, min(nfft.size()/N, nw)));
    Info<< "    Creating graph for P(f)" << endl;
    writer->write
    (
        outputPath,
        "Pf",
        coordSet(true, "f [Hz]", Pf.first()),
        "P(f) [Pa]",
        Pf.second()
    );

    const Pair<scalarField> Lf(nfft.Lf(Pf));
    Info<< "    Creating graph for L(f)" << endl;
    writer->write
    (
        outputPath,
        "Lf",
        coordSet(true, "f [Hz]", Lf.first()),
        "L(f) [dB]",
        Lf.second()
    );

    const Pair<scalarField> Ldelta(nfft.Ldelta(Lf, f1, fU));
    Info<< "    Creating graph for Ldelta" << endl;
    writer->write
    (
        outputPath,
        "Ldelta",
        coordSet(true, "fm [Hz]", Ldelta.first()),
        "Ldelta(f) [dB]",
        Ldelta.second()
    );

    const Pair<scalarField> Pdelta(nfft.Pdelta(Pf, f1, fU));
    Info<< "    Creating graph for Pdelta" << endl;
    writer->write
    (
        outputPath,
        "Pdelta",
        coordSet(true, "fm [Hz]", Pdelta.first()),
        "P(f) [dB]",
        Pdelta.second()
    );

    Info<< nl << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
