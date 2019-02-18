/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2019 OpenFOAM Foundation
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
    const dictionary& dict,
    const bool pRefRequired
)
:
    refCell_(-1),
    refValue_(0),
    pMax_("pMax", dimPressure, great),
    pMin_("pMin", dimPressure, 0),
    limitMaxP_(false),
    limitMinP_(false)
{
    bool pLimits = false;
    scalar pMax = -great;
    scalar pMin = great;

    // Set the reference cell and value for closed domain simulations
    if (pRefRequired && setRefCell(p, dict, refCell_, refValue_))
    {
        pLimits = true;

        pMax = refValue_;
        pMin = refValue_;
    }

    if (dict.found("pMax") && dict.found("pMin"))
    {
        pMax_.value() = readScalar(dict.lookup("pMax"));
        limitMaxP_ = true;
        pMin_.value() = readScalar(dict.lookup("pMin"));
        limitMinP_ = true;
    }
    else
    {
        const volScalarField::Boundary& pbf = p.boundaryField();
        const volScalarField::Boundary& rhobf = rho.boundaryField();

        scalar rhoRefMax = -great;
        scalar rhoRefMin = great;
        bool rhoLimits = false;

        forAll(pbf, patchi)
        {
            if (pbf[patchi].fixesValue())
            {
                pLimits = true;
                rhoLimits = true;

                pMax = max(pMax, max(pbf[patchi]));
                pMin = min(pMin, min(pbf[patchi]));

                rhoRefMax = max(rhoRefMax, max(rhobf[patchi]));
                rhoRefMin = min(rhoRefMin, min(rhobf[patchi]));
            }
        }

        reduce(rhoLimits, andOp<bool>());
        if (rhoLimits)
        {
            reduce(pMax, maxOp<scalar>());
            reduce(pMin, minOp<scalar>());

            reduce(rhoRefMax, maxOp<scalar>());
            reduce(rhoRefMin, minOp<scalar>());
        }

        if (dict.found("pMax"))
        {
            pMax_.value() = readScalar(dict.lookup("pMax"));
            limitMaxP_ = true;
        }
        else if (dict.found("pMaxFactor"))
        {
            if (!pLimits)
            {
                FatalIOErrorInFunction(dict)
                    << "'pMaxFactor' specified rather than 'pMax'" << nl
                    << "    but the corresponding reference pressure cannot"
                       " be evaluated from the boundary conditions." << nl
                    << "    Please specify 'pMax' rather than 'pMaxFactor'"
                    << exit(FatalIOError);
            }

            const scalar pMaxFactor(readScalar(dict.lookup("pMaxFactor")));
            pMax_.value() = pMaxFactor*pMax;
            limitMaxP_ = true;
        }
        else if (dict.found("rhoMax"))
        {
            // For backward-compatibility infer the pMax from rhoMax

            IOWarningInFunction(dict)
                 << "'rhoMax' specified rather than 'pMax' or 'pMaxFactor'"
                 << nl
                 << "    This is supported for backward-compatibility but "
                    "'pMax' or 'pMaxFactor' are more reliable." << endl;

            if (!pLimits)
            {
                FatalIOErrorInFunction(dict)
                    << "'rhoMax' specified rather than 'pMax'" << nl
                    << "    but the corresponding reference pressure cannot"
                       " be evaluated from the boundary conditions." << nl
                    << "    Please specify 'pMax' rather than 'rhoMax'"
                    << exit(FatalIOError);
            }

            if (!rhoLimits)
            {
                FatalIOErrorInFunction(dict)
                    << "'rhoMax' specified rather than 'pMaxFactor'" << nl
                    << "    but the corresponding reference density cannot"
                       " be evaluated from the boundary conditions." << nl
                    << "    Please specify 'pMaxFactor' rather than 'rhoMax'"
                    << exit(FatalIOError);
            }

            dimensionedScalar rhoMax("rhoMax", dimDensity, dict);

            pMax_.value() = max(rhoMax.value()/rhoRefMax, 1)*pMax;
            limitMaxP_ = true;
        }

        if (dict.found("pMin"))
        {
            pMin_.value() = readScalar(dict.lookup("pMin"));
            limitMinP_ = true;
        }
        else if (dict.found("pMinFactor"))
        {
            if (!pLimits)
            {
                FatalIOErrorInFunction(dict)
                    << "'pMinFactor' specified rather than 'pMin'" << nl
                    << "    but the corresponding reference pressure cannot"
                       " be evaluated from the boundary conditions." << nl
                    << "    Please specify 'pMin' rather than 'pMinFactor'"
                    << exit(FatalIOError);
            }

            const scalar pMinFactor(readScalar(dict.lookup("pMinFactor")));
            pMin_.value() = pMinFactor*pMin;
            limitMinP_ = true;
        }
        else if (dict.found("rhoMin"))
        {
            // For backward-compatibility infer the pMin from rhoMin

            IOWarningInFunction(dict)
                << "'rhoMin' specified rather than 'pMin' or 'pMinFactor'" << nl
                << "    This is supported for backward-compatibility but"
                   "'pMin' or 'pMinFactor' are more reliable." << endl;

            if (!pLimits)
            {
                FatalIOErrorInFunction(dict)
                    << "'rhoMin' specified rather than 'pMin'" << nl
                    << "    but the corresponding reference pressure cannot"
                       " be evaluated from the boundary conditions." << nl
                    << "    Please specify 'pMin' rather than 'rhoMin'"
                    << exit(FatalIOError);
            }

            if (!rhoLimits)
            {
                FatalIOErrorInFunction(dict)
                    << "'rhoMin' specified rather than 'pMinFactor'" << nl
                    << "    but the corresponding reference density cannot"
                       " be evaluated from the boundary conditions." << nl
                    << "    Please specify 'pMinFactor' rather than 'rhoMin'"
                    << exit(FatalIOError);
            }

            dimensionedScalar rhoMin("rhoMin", dimDensity, dict);

            pMin_.value() = min(rhoMin.value()/rhoRefMin, 1)*pMin;
            limitMinP_ = true;
        }
    }

    if (limitMaxP_ || limitMinP_)
    {
        Info<< "pressureControl" << nl;

        if (limitMaxP_)
        {
            Info<< "    pMax " << pMax_.value() << nl;
        }

        if (limitMinP_)
        {
            Info<< "    pMin " << pMin_.value() << nl;
        }

        Info << endl;
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::pressureControl::limit(volScalarField& p) const
{
    if (limitMaxP_ || limitMinP_)
    {
        if (limitMaxP_)
        {
            const scalar pMax = max(p).value();

            if (pMax > pMax_.value())
            {
                Info<< "pressureControl: p max " << pMax << endl;
                p = min(p, pMax_);
            }
        }

        if (limitMinP_)
        {
            const scalar pMin = min(p).value();

            if (pMin < pMin_.value())
            {
                Info<< "pressureControl: p min " << pMin << endl;
                p = max(p, pMin_);
            }
        }

        p.correctBoundaryConditions();

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
