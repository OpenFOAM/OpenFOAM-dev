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
    equilibriumCO

Description
    Calculates the equilibrium level of carbon monoxide.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "dictionary.H"
#include "IFstream.H"
#include "OSspecific.H"
#include "IOmanip.H"

#include "specie.H"
#include "perfectGas.H"
#include "thermo.H"
#include "janafThermo.H"
#include "absoluteEnthalpy.H"

#include "SLPtrList.H"

using namespace Foam;

typedef species::thermo<janafThermo<perfectGas<specie> >, absoluteEnthalpy>
    thermo;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"

    Info<< nl << "Reading thermodynamic data IOdictionary" << endl;

    IOdictionary thermoData
    (
        IOobject
        (
            "thermoData",
            runTime.constant(),
            runTime,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE,
            false
        )
    );



    scalar P = 1e5;
    scalar T = 3000.0;

    SLPtrList<thermo> EQreactions;

    EQreactions.append
    (
        new thermo
        (
            thermo(thermoData.subDict("CO2"))
         ==
            thermo(thermoData.subDict("CO"))
          + 0.5*thermo(thermoData.subDict("O2"))
        )
    );

    EQreactions.append
    (
        new thermo
        (
            thermo(thermoData.subDict("O2"))
         ==
            2.0*thermo(thermoData.subDict("O"))
        )
    );

    EQreactions.append
    (
        new thermo
        (
            thermo(thermoData.subDict("H2O"))
         ==
            thermo(thermoData.subDict("H2"))
          + 0.5*thermo(thermoData.subDict("O2"))
        )
    );

    EQreactions.append
    (
        new thermo
        (
            thermo(thermoData.subDict("H2O"))
         ==
            thermo(thermoData.subDict("H"))
          + thermo(thermoData.subDict("OH"))
        )
    );


    forAllConstIter(SLPtrList<thermo>, EQreactions, iter)
    {
        Info<< "Kc(EQreactions) = " << iter().Kc(P, T) << endl;
    }

    Info<< nl << "end" << endl;

    return 0;
}


// ************************************************************************* //
