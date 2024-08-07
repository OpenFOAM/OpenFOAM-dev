/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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

Info << nl << "Creating autocorrelation functions." << endl;

IOdictionary mdTransportPropertiesDict
(
    IOobject
    (
        "mdTransportPropertiesDict",
        mesh.time().system(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE,
        false
    )
);

const dictionary& autocorrelationFunctionDict
(
    mdTransportPropertiesDict.subDict("autocorrelationFunctions")
);

//- Velocity autocorrelation function

Info << tab << "velocity" << endl;

const dictionary& velocityACFDict
(
    autocorrelationFunctionDict.subDict("velocity")
);

correlationFunction<vector> vacf
(
    mesh,
    velocityACFDict,
    molecules.size()
);

bool writeVacf(Switch(velocityACFDict.lookup("writeFile")));

//- Pressure autocorrelation function

Info << tab << "pressure" << endl;

const dictionary& pressureACFDict
(
    autocorrelationFunctionDict.subDict("pressure")
);

correlationFunction<vector> pacf
(
    mesh,
    pressureACFDict,
    1
);

bool writePacf(Switch(pressureACFDict.lookup("writeFile")));

//- Heat flux autocorrelation function

Info << tab << "heat flux" << endl;

const dictionary& heatFluxACFDict
(
    autocorrelationFunctionDict.subDict("heatFlux")
);

correlationFunction<vector> hfacf
(
    mesh,
    heatFluxACFDict,
    1
);

bool writeHFacf(Switch(heatFluxACFDict.lookup("writeFile")));
