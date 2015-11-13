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

Application
    createTurbulenceFields

Description
    Creates a full set of turbulence fields.

    - Currently does not output nut and nuTilda

Source files:
    createFields.H

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    timeSelector::addOptions();

    argList::addOption
    (
        "fields",
        "wordReList",
        "specify which turbulence fields (k, epsilon, omega, R) to write"
        " - eg '(k omega)' or '(R)' or '(.*)'."
    );

    #include "setRootCase.H"
    #include "createTime.H"

    instantList timeDirs = timeSelector::select0(runTime, args);

    const bool selectedFields = args.optionFound("fields");
    wordReList fieldPatterns;
    if (selectedFields)
    {
        fieldPatterns = wordReList(args.optionLookup("fields")());
    }

    #include "createMesh.H"

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);

        Info<< "Time = " << runTime.timeName() << endl;

        #include "createFields.H"

        if (findStrings(fieldPatterns, "k"))
        {
            if (!IOobject("k", runTime.timeName(), mesh).headerOk())
            {
                Info<< "    Writing turbulence field k" << endl;
                const volScalarField k(RASModel->k());
                k.write();
            }
            else
            {
                Info<< "    Turbulence k field already exists" << endl;
            }
        }

        if (findStrings(fieldPatterns, "epsilon"))
        {
            if (!IOobject("epsilon", runTime.timeName(), mesh).headerOk())
            {
                Info<< "    Writing turbulence field epsilon" << endl;
                const volScalarField epsilon(RASModel->epsilon());
                epsilon.write();
            }
            else
            {
                Info<< "    Turbulence epsilon field already exists" << endl;
            }
        }

        if (findStrings(fieldPatterns, "R"))
        {
            if (!IOobject("R", runTime.timeName(), mesh).headerOk())
            {
                Info<< "    Writing turbulence field R" << endl;
                const volSymmTensorField R(RASModel->R());
                R.write();
            }
            else
            {
                Info<< "    Turbulence R field already exists" << endl;
            }
        }

        if (findStrings(fieldPatterns, "omega"))
        {
            if (!IOobject("omega", runTime.timeName(), mesh).headerOk())
            {
                const scalar Cmu = 0.09;

                // Assume k and epsilon are available
                const volScalarField k(RASModel->k());
                const volScalarField epsilon(RASModel->epsilon());

                volScalarField omega
                (
                    IOobject
                    (
                        "omega",
                        runTime.timeName(),
                        mesh
                    ),
                    epsilon/(Cmu*k),
                    epsilon.boundaryField().types()
                );

                Info<< "    Writing turbulence field omega" << endl;
                omega.write();
            }
            else
            {
                Info<< "    Turbulence omega field already exists" << endl;
            }
        }
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
