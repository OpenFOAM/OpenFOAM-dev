/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021 OpenFOAM Foundation
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

#include "limitPressure.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(limitPressure, 0);
    addToRunTimeSelectionTable
    (
        fvConstraint,
        limitPressure,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::limitPressure::readCoeffs()
{
    const dictionary& dict(coeffs());

    pName_ = dict.lookupOrDefault<word>("p", "p");

    if (dict.found("min") && dict.found("max"))
    {
        pMin_.value() = dict.lookup<scalar>("min");
        limitMinP_ = true;

        pMax_.value() = dict.lookup<scalar>("max");
        limitMaxP_ = true;
    }
    else
    {
        const volScalarField& p = mesh().lookupObject<volScalarField>(pName_);
        const volScalarField::Boundary& pbf = p.boundaryField();

        bool pLimits = false;
        scalar pMin = vGreat;
        scalar pMax = -vGreat;

        forAll(pbf, patchi)
        {
            if
            (
                pbf[patchi].fixesValue()
             || isA<calculatedFvPatchField<scalar>>(pbf[patchi])
            )
            {
                pLimits = true;

                pMin = min(pMin, min(pbf[patchi]));
                pMax = max(pMax, max(pbf[patchi]));
            }
        }

        reduce(pLimits, andOp<bool>());
        if (pLimits)
        {
            reduce(pMax, maxOp<scalar>());
            reduce(pMin, minOp<scalar>());
        }

        if (dict.found("min"))
        {
            pMin_.value() = dict.lookup<scalar>("min");
            limitMinP_ = true;
        }
        else if (dict.found("minFactor"))
        {
            if (!pLimits)
            {
                FatalIOErrorInFunction(dict)
                    << "'minFactor' specified rather than 'min'" << nl
                    << "    but the corresponding reference pressure cannot"
                       " be evaluated from the boundary conditions." << nl
                    << "    Please specify 'min' rather than 'minFactor'"
                    << exit(FatalIOError);
            }

            const scalar pMinFactor(dict.lookup<scalar>("minFactor"));
            pMin_.value() = pMinFactor*pMin;
            limitMinP_ = true;
        }

        if (dict.found("max"))
        {
            pMax_.value() = dict.lookup<scalar>("max");
            limitMaxP_ = true;
        }
        else if (dict.found("maxFactor"))
        {
            if (!pLimits)
            {
                FatalIOErrorInFunction(dict)
                    << "'maxFactor' specified rather than 'max'" << nl
                    << "    but the corresponding reference pressure cannot"
                       " be evaluated from the boundary conditions." << nl
                    << "    Please specify 'max' rather than 'maxFactor'"
                    << exit(FatalIOError);
            }

            const scalar pMaxFactor(dict.lookup<scalar>("maxFactor"));
            pMax_.value() = pMaxFactor*pMax;
            limitMaxP_ = true;
        }
    }

    if (limitMinP_ || limitMaxP_)
    {
        if (limitMinP_)
        {
            Info<< "    min " << pMin_.value() << nl;
        }

        if (limitMaxP_)
        {
            Info<< "    max " << pMax_.value() << nl;
        }

        Info << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::limitPressure::limitPressure
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    fvConstraint(name, modelType, dict, mesh),
    pName_(word::null),
    pMin_("pMin", dimPressure, 0),
    pMax_("pMax", dimPressure, great),
    limitMinP_(false),
    limitMaxP_(false)
{
    readCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::wordList Foam::fv::limitPressure::constrainedFields() const
{
    return wordList(1, pName_);
}


bool Foam::fv::limitPressure::constrain(volScalarField& p) const
{
    bool constrained = false;

    if (limitMinP_ || limitMaxP_)
    {
        if (limitMinP_)
        {
            const scalar pMin = min(p).value();

            if (pMin < pMin_.value())
            {
                Info<< "limitPressure: min " << pMin << endl;
                p = max(p, pMin_);
                constrained = true;
            }
        }

        if (limitMaxP_)
        {
            const scalar pMax = max(p).value();

            if (pMax > pMax_.value())
            {
                Info<< "limitPressure: max " << pMax << endl;
                p = min(p, pMax_);
                constrained = true;
            }
        }

        p.correctBoundaryConditions();
    }

    return constrained;
}


void Foam::fv::limitPressure::updateMesh(const mapPolyMesh&)
{}


void Foam::fv::limitPressure::distribute(const mapDistributePolyMesh&)
{}


bool Foam::fv::limitPressure::movePoints()
{
    return true;
}



bool Foam::fv::limitPressure::read(const dictionary& dict)
{
    if (fvConstraint::read(dict))
    {
        readCoeffs();
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
