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

Application
    noise

Description
    Utility to perform noise analysis of pressure data using the noiseFFT
    library.

    Control settings are read from the $FOAM_CASE/system/noiseDict dictionary,
    or user-specified dictionary using the -dict option.  Pressure data is
    read using a CSV reader:

    \heading Usage

    \verbatim
    pRef        101325;
    N           65536;
    nw          100;
    f1          25;
    fU          10000;
    graphFormat raw;

    csvFileData
    {
        fileName        "pressureData"
        nHeaderLine     1;
        refColumn       0;
        componentColumns (1);
        separator       " ";
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

SeeAlso
    CSV.H
    noiseFFT.H

\*---------------------------------------------------------------------------*/


#include "noiseFFT.H"
#include "argList.H"
#include "Time.H"
#include "functionObjectFile.H"
#include "CSV.H"

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
            scalar dT = t[i] - t[i-1];
            if (deltaT < 0)
            {
                deltaT = dT;
            }

            if (mag(deltaT - dT) > SMALL)
            {
                FatalErrorIn("checkUniformTimeStep(const scalarField&)")
                    << "Unable to process data with a variable time step"
                    << exit(FatalError);
            }
        }
    }
    else
    {
        FatalErrorIn("checkUniformTimeStep(const scalarField&)")
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
    #include "createFields.H"

    Info<< "Reading data file" << endl;
    CSV<scalar> pData("pressure", dict, "Data");

    // time history data
    const scalarField t(pData.x());

    // pressure data
    const scalarField p(pData.y());

    if (t.size() < N)
    {
        FatalErrorIn(args.executable())
            << "Block size N = " << N
            << " is larger than number of data = " << t.size()
            << exit(FatalError);
    }

    Info<< "    read " << t.size() << " values" << nl << endl;


    Info<< "Creating noise FFT" << endl;
    noiseFFT nfft(checkUniformTimeStep(t), p);

    nfft -= pRef;

    fileName baseFileName(pData.fName().lessExt());

    graph Pf(nfft.RMSmeanPf(N, min(nfft.size()/N, nw)));
    Info<< "    Creating graph for " << Pf.title() << endl;
    Pf.write(baseFileName + graph::wordify(Pf.title()), graphFormat);

    graph Lf(nfft.Lf(Pf));
    Info<< "    Creating graph for " << Lf.title() << endl;
    Lf.write(baseFileName + graph::wordify(Lf.title()), graphFormat);

    graph Ldelta(nfft.Ldelta(Lf, f1, fU));
    Info<< "    Creating graph for " << Ldelta.title() << endl;
    Ldelta.write(baseFileName + graph::wordify(Ldelta.title()), graphFormat);

    graph Pdelta(nfft.Pdelta(Pf, f1, fU));
    Info<< "    Creating graph for " << Pdelta.title() << endl;
    Pdelta.write(baseFileName + graph::wordify(Pdelta.title()), graphFormat);

    Info<< nl << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
