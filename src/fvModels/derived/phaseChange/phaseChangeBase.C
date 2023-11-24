/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2023 OpenFOAM Foundation
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

#include "phaseChangeBase.H"
#include "fluidThermo.H"
#include "multicomponentThermo.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(phaseChangeBase, 0);
}
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

const Foam::volScalarField& Foam::fv::phaseChangeBase::p() const
{
    for (label i = 0; i < 2; ++ i)
    {
        if (isA<fluidThermo>(thermos()[i]))
        {
            return refCast<const fluidThermo>(thermos()[i]).p();
        }
    }

    return mesh().lookupObject<volScalarField>("p");
}


Foam::tmp<Foam::volScalarField>
Foam::fv::phaseChangeBase::vifToVf
(
    const tmp<volScalarField::Internal>& tvif
)
{
    tmp<volScalarField> tvf =
        volScalarField::New
        (
            tvif().name(),
            tvif().mesh(),
            tvif().dimensions(),
            extrapolatedCalculatedFvPatchField<scalar>::typeName
        );

    tvf->ref() = tvif();
    tvf->correctBoundaryConditions();

    tvif.clear();

    return tvf;
}


Foam::tmp<Foam::volScalarField::Internal>
Foam::fv::phaseChangeBase::vfToVif
(
    const tmp<volScalarField>& tvf
)
{
    tmp<volScalarField::Internal> tvif(tvf.ptr());

    tvf.clear();

    return tvif;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::phaseChangeBase::phaseChangeBase
(
    const word& name,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict,
    const Pair<bool> specieThermosRequired
)
:
    massTransferBase(name, modelType, mesh, dict),
    thermos_(mesh, phaseNames()),
    specieThermos_(thermos_),
    heNames_(thermos_.first().he().name(), thermos_.second().he().name())
{
    forAll(specieThermos_.valid(), i)
    {
        if (!specieThermos_.valid()[i] && specieThermosRequired[i])
        {
            FatalErrorInFunction
                << "Model " << name << " of type " << modelType
                << " requires a multicomponent thermo for phase "
                << phaseNames()[i] << exit(FatalError);
        }
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField::Internal>
Foam::fv::phaseChangeBase::Tchange() const
{
    const volScalarField::Internal mDot(this->mDot());

    return pos0(mDot)*thermos().first().T() + neg(mDot)*thermos().second().T();
}


Foam::tmp<Foam::volScalarField::Internal>
Foam::fv::phaseChangeBase::Lfraction() const
{
    const volScalarField& kappa1 = thermos().first().kappa();
    const volScalarField& kappa2 = thermos().second().kappa();

    return vfToVif(kappa2/(kappa1 + kappa2));
}


// ************************************************************************* //
