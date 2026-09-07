/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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

#include "limitGroupSurfaceAreaVolumeRatio.H"
#include "fractal.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(limitGroupSurfaceAreaVolumeRatio, 0);
    addToRunTimeSelectionTable
    (
        fvConstraint,
        limitGroupSurfaceAreaVolumeRatio,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::limitGroupSurfaceAreaVolumeRatio::readCoeffs
(
    const dictionary& dict
)
{
    if (!dict.found("minDiameter"))
    {
        Info<< indent << "minDiameter not specified" << nl
            << indent << "limiting to the smallest diameter in the population"
            << " (= " << popBal_.dSph(0).value() << ")" << endl;
    }

    minDiameter_.readOrDefault(dict, popBal_.dSph(0).value());
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::limitGroupSurfaceAreaVolumeRatio::limitGroupSurfaceAreaVolumeRatio
(
    const word& name,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    fvConstraint(name, modelType, mesh, dict),
    popBal_
    (
        mesh().lookupObject<populationBalanceModel>
        (
            coeffs(dict).lookup<word>("populationBalance")
        )
    ),
    minDiameter_("minDiameter", dimLength, 0),
    constrainedFields_()
{
    readCoeffs(coeffs(dict));

    const populationBalance::shapeModel& shape = popBal_.shape();

    if (isA<populationBalance::shapeModels::fractal>(shape))
    {
        const populationBalance::shapeModels::fractal& fractal =
            refCast<const populationBalance::shapeModels::fractal>(shape);

        constrainedFields_.setSize(fractal.kappas().size());

        forAll(fractal.kappas(), i)
        {
            constrainedFields_[i] = fractal.kappas()[i].name();
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::wordList
Foam::fv::limitGroupSurfaceAreaVolumeRatio::constrainedFields() const
{
    return constrainedFields_;
}


bool Foam::fv::limitGroupSurfaceAreaVolumeRatio::constrain
(
    VolField<scalar>& kappa
) const
{
    const populationBalance::shapeModels::fractal& fractal =
        refCast<const populationBalance::shapeModels::fractal>(popBal_.shape());

    forAll(fractal.kappas(), i)
    {
        if (&kappa == &fractal.kappas()[i])
        {
            kappa = min(max(kappa, 6/popBal_.dSph(i)), 6/minDiameter_);
        }
    }

    return true;
}


bool Foam::fv::limitGroupSurfaceAreaVolumeRatio::movePoints()
{
    return true;
}


void Foam::fv::limitGroupSurfaceAreaVolumeRatio::topoChange
(
    const polyTopoChangeMap&
)
{}


void Foam::fv::limitGroupSurfaceAreaVolumeRatio::mapMesh
(
    const polyMeshMap&
)
{}


void Foam::fv::limitGroupSurfaceAreaVolumeRatio::distribute
(
    const polyDistributionMap&
)
{}


bool Foam::fv::limitGroupSurfaceAreaVolumeRatio::read(const dictionary& dict)
{
    if (fvConstraint::read(dict))
    {
        readCoeffs(coeffs(dict));
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
