/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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

Description

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "IFstream.H"
#include "OFstream.H"
#include "liquidProperties.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::validArgs.append("liquidName");
    argList::validArgs.append("pMin");
    argList::validArgs.append("pMax");
    argList::validArgs.append("nP");
    argList::validArgs.append("Tmin");
    argList::validArgs.append("Tmax");
    argList::validArgs.append("nT");
    argList args(argc, argv);
    const word liquidName(args[1]);

    const scalar pMin = args.argRead<scalar>(2);
    const scalar pMax = args.argRead<scalar>(3);
    const scalar nP = args.argRead<label>(4);
    const scalar Tmin = args.argRead<scalar>(5);
    const scalar Tmax = args.argRead<scalar>(6);
    const scalar nT = args.argRead<label>(7);

    autoPtr<liquidProperties> liquidPtr = liquidProperties::New(liquidName);

    OFstream plotFile(liquidName + ".dat");

    plotFile << "# p T rho Cp Hs Ha pv hl Cpg mu mug kappa kappag sigma" << nl;

    for (label pi = 0; pi < nP; ++ pi)
    {
        const scalar p = pMin + (pMax - pMin)*pi/(nP - 1);

        for (label Ti = 0; Ti < nT; ++ Ti)
        {
            const scalar T = Tmin + (Tmax - Tmin)*Ti/(nT - 1);

            plotFile
                << p << ' '
                << T << ' '
                << liquidPtr->rho(p, T) << ' '
                << liquidPtr->Cp(p, T) << ' '
                << liquidPtr->Hs(p, T) << ' '
                << liquidPtr->Ha(p, T) << ' '
                << liquidPtr->pv(p, T) << ' '
                << liquidPtr->hl(p, T) << ' '
                << liquidPtr->Cpg(p, T) << ' '
                << liquidPtr->mu(p, T) << ' '
                << liquidPtr->mug(p, T) << ' '
                << liquidPtr->kappa(p, T) << ' '
                << liquidPtr->kappag(p, T) << ' '
                << liquidPtr->sigma(p, T)
                << nl;
        }

        plotFile << nl;
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
