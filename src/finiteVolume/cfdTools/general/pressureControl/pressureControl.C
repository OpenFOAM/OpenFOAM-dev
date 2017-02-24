/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenFOAM Foundation
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

#include "pressureControl.H"
#include "findRefCell.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pressureControl::pressureControl
(
    const volScalarField& p,
    const volScalarField& rho,
    const dictionary& dict
)
:
    refCell_(0),
    refValue_(0),
    pMax_("pMax", dimPressure, 0),
    pMin_("pMin", dimPressure, GREAT)
{
    if (setRefCell(p, dict, refCell_, refValue_))
    {
        pMax_.value() = refValue_;
        pMin_.value() = refValue_;
    }

    const volScalarField::Boundary& pbf = p.boundaryField();
    const volScalarField::Boundary& rhobf = rho.boundaryField();

    scalar rhoRefMax = -GREAT;
    scalar rhoRefMin = GREAT;
    bool rhoLimits = false;

    forAll(pbf, patchi)
    {
        if (pbf[patchi].fixesValue())
        {
            rhoLimits = true;

            pMax_.value() = max(pMax_.value(), max(pbf[patchi]));
            pMin_.value() = min(pMin_.value(), min(pbf[patchi]));

            rhoRefMax = max(rhoRefMax, max(rhobf[patchi]));
            rhoRefMin = min(rhoRefMin, min(rhobf[patchi]));
        }
    }

    if (dict.found("pMax"))
    {
        pMax_.value() = readScalar(dict.lookup("pMax"));
    }
    else if (dict.found("pMaxFactor"))
    {
        const scalar pMaxFactor(readScalar(dict.lookup("pMaxFactor")));
        pMax_ *= pMaxFactor;
    }
    else if (dict.found("rhoMax"))
    {
        // For backward-compatibility infer the pMax from rhoMax

        IOWarningInFunction(dict)
            << "'rhoMax' specified rather than 'pMax' or 'pMaxFactor'" << nl
            << "    This is supported for backward-compatibility but "
               "'pMax' or 'pMaxFactor' are more reliable." << endl;

        if (!rhoLimits)
        {
            FatalIOErrorInFunction(dict)
                << "'rhoMax' specified rather than 'pMaxFactor'" << nl
                << "    but the corresponding reference density cannot"
                   " be evaluated from the boundary conditions." << nl
                << "Please specify 'pMaxFactor' rather than 'rhoMax'"
                << exit(FatalError);
        }

        dimensionedScalar rhoMax("rhoMax", dimDensity, dict);

        pMax_ *= max(rhoMax.value()/rhoRefMax, 1);
    }

    if (dict.found("pMin"))
    {
        pMin_.value() = readScalar(dict.lookup("pMin"));
    }
    else if (dict.found("pMinFactor"))
    {
        const scalar pMinFactor(readScalar(dict.lookup("pMinFactor")));
        pMin_ *= pMinFactor;
    }
    else if (dict.found("rhoMin"))
    {
        // For backward-compatibility infer the pMin from rhoMin

        IOWarningInFunction(dict)
            << "'rhoMin' specified rather than 'pMin' or 'pMinFactor'" << nl
            << "    This is supported for backward-compatibility but"
               "'pMin' or 'pMinFactor' are more reliable." << endl;

        if (!rhoLimits)
        {
            FatalIOErrorInFunction(dict)
                << "'rhoMin' specified rather than 'pMinFactor'" << nl
                << "    but the corresponding reference density cannot"
                   " be evaluated from the boundary conditions." << nl
                << "Please specify 'pMinFactor' rather than 'rhoMin'"
                << exit(FatalError);
        }

        dimensionedScalar rhoMin("rhoMin", dimDensity, dict);

        pMin_ *= min(rhoMin.value()/rhoRefMin, 1);
    }

    Info<< "pressureControl" << nl
        << "    pMax/pMin " << pMax_.value() << " " << pMin_.value()
        << nl << endl;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::pressureControl::limit(volScalarField& p) const
{
    Info<< "pressureControl: p max/min "
        << max(p).value() << " "
        << min(p).value() << endl;

    p = max(p, pMin_);
    p = min(p, pMax_);
}


// ************************************************************************* //
