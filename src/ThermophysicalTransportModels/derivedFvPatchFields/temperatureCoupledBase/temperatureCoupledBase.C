/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

#include "temperatureCoupledBase.H"
#include "fluidThermo.H"
#include "solidThermo.H"
#include "thermophysicalTransportModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::temperatureCoupledBase::temperatureCoupledBase
(
    const fvPatch& patch
)
:
    patch_(patch),
    alphaAniName_(word::null)
{}


Foam::temperatureCoupledBase::temperatureCoupledBase
(
    const fvPatch& patch,
    const dictionary& dict
)
:
    patch_(patch),
    alphaAniName_(dict.lookupOrDefault<word>("alphaAni", word::null))
{}


Foam::temperatureCoupledBase::temperatureCoupledBase
(
    const fvPatch& patch,
    const temperatureCoupledBase& base
)
:
    patch_(patch),
    alphaAniName_(base.alphaAniName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField> Foam::temperatureCoupledBase::kappa
(
    const fvPatchScalarField& Tp
) const
{
    const fvMesh& mesh = patch_.boundaryMesh().mesh();
    const label patchi = patch_.index();

    const word& phase(Tp.internalField().group());

    const word fluidThermoName
    (
        IOobject::groupName(basicThermo::dictName, phase)
    );

    if (mesh.foundObject<fluidThermo>(fluidThermoName))
    {
        static word ttmName
        (
            IOobject::groupName
            (
                thermophysicalTransportModel::typeName,
                phase
            )
        );

        if (mesh.foundObject<thermophysicalTransportModel>(ttmName))
        {
            const thermophysicalTransportModel& ttm =
                mesh.lookupObject<thermophysicalTransportModel>(ttmName);

            return ttm.kappaEff(patchi);
        }
        else
        {
            const fluidThermo& thermo =
                mesh.lookupObject<fluidThermo>(fluidThermoName);

            return thermo.kappa(patchi);
        }
    }
    else if (mesh.foundObject<solidThermo>(basicThermo::dictName))
    {
        const solidThermo& thermo =
            mesh.lookupObject<solidThermo>(basicThermo::dictName);

        if (alphaAniName_ != word::null)
        {
            const symmTensorField& alphaAni =
                patch_.lookupPatchField<volSymmTensorField, scalar>
                (
                    alphaAniName_
                );

            const symmTensorField kappa(alphaAni*thermo.Cp(Tp, patchi));
            const vectorField n(patch_.nf());

            return n & kappa & n;
        }
        else
        {
            return thermo.kappa(patchi);
        }
    }
    else
    {
        FatalErrorInFunction
            << "Cannot find a fluidThermo or solidThermo instance"
            << exit(FatalError);

        return scalarField::null();
    }
}


void Foam::temperatureCoupledBase::write(Ostream& os) const
{
    if (alphaAniName_ != word::null)
    {
        writeEntry(os, "alphaAni", alphaAniName_);
    }
}


// ************************************************************************* //
