 /*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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
#include "volFields.H"
#include "fluidThermo.H"
#include "solidThermo.H"
#include "turbulenceModel.H"

// * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * * //

namespace Foam
{
    template<>
    const char* Foam::NamedEnum
    <
        Foam::temperatureCoupledBase::KMethodType,
        4
    >::names[] =
    {
        "fluidThermo",
        "solidThermo",
        "directionalSolidThermo",
        "lookup"
    };
}


const Foam::NamedEnum<Foam::temperatureCoupledBase::KMethodType, 4>
    Foam::temperatureCoupledBase::KMethodTypeNames_;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::temperatureCoupledBase::temperatureCoupledBase
(
    const fvPatch& patch,
    const word& calculationType,
    const word& kappaName,
    const word& alphaAniName
)
:
    patch_(patch),
    method_(KMethodTypeNames_[calculationType]),
    kappaName_(kappaName),
    alphaAniName_(alphaAniName)
{}


Foam::temperatureCoupledBase::temperatureCoupledBase
(
    const fvPatch& patch,
    const dictionary& dict
)
:
    patch_(patch),
    method_(KMethodTypeNames_.read(dict.lookup("kappa"))),
    kappaName_(dict.lookup("kappaName")),
    alphaAniName_(dict.lookupOrDefault<word>("alphaAniName","Anialpha"))
{}


Foam::temperatureCoupledBase::temperatureCoupledBase
(
    const fvPatch& patch,
    const temperatureCoupledBase& base
)
:
    patch_(patch),
    method_(base.method_),
    kappaName_(base.kappaName_),
    alphaAniName_(base.alphaAniName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField> Foam::temperatureCoupledBase::kappa
(
    const scalarField& Tp
) const
{
    const fvMesh& mesh = patch_.boundaryMesh().mesh();
    const label patchI = patch_.index();

    switch (method_)
    {
        case mtFluidThermo:
        {
            typedef compressible::turbulenceModel turbulenceModel;

            if (mesh.foundObject<turbulenceModel>("turbulenceModel"))
            {
                const turbulenceModel& turbModel =
                    mesh.lookupObject<turbulenceModel>("turbulenceModel");

                return turbModel.kappaEff(patchI);
            }
            else if (mesh.foundObject<fluidThermo>("thermophysicalProperties"))
            {
                const fluidThermo& thermo =
                    mesh.lookupObject<fluidThermo>("thermophysicalProperties");

                return thermo.kappa(patchI);
            }
            else
            {
                FatalErrorIn
                (
                    "temperatureCoupledBase::kappa(const scalarField&) const"
                )
                    << "Kappa defined to employ " << KMethodTypeNames_[method_]
                    << " method, but thermo package not available"
                    << exit(FatalError);
            }

            break;
        }

        case mtSolidThermo:
        {
            const solidThermo& thermo =
                mesh.lookupObject<solidThermo>("thermophysicalProperties");

            return thermo.kappa(patchI);
            break;
        }

        case mtDirectionalSolidThermo:
        {
            const solidThermo& thermo =
                mesh.lookupObject<solidThermo>("thermophysicalProperties");

            const symmTensorField& alphaAni =
                patch_.lookupPatchField<volSymmTensorField, scalar>
                (
                    alphaAniName_
                );

            const scalarField& pp = thermo.p().boundaryField()[patchI];

            const symmTensorField kappa(alphaAni*thermo.Cp(pp, Tp, patchI));

            const vectorField n(patch_.nf());

            return n & kappa & n;
        }

        case mtLookup:
        {
            if (mesh.foundObject<volScalarField>(kappaName_))
            {
                return patch_.lookupPatchField<volScalarField, scalar>
                (
                    kappaName_
                );
            }
            else if (mesh.foundObject<volSymmTensorField>(kappaName_))
            {
                const symmTensorField& KWall =
                    patch_.lookupPatchField<volSymmTensorField, scalar>
                    (
                        kappaName_
                    );

                const vectorField n(patch_.nf());

                return n & KWall & n;
            }
            else
            {
                FatalErrorIn
                (
                    "temperatureCoupledBase::kappa(const scalarField&) const"
                )
                    << "Did not find field " << kappaName_
                    << " on mesh " << mesh.name() << " patch " << patch_.name()
                    << nl
                    << "Please set 'kappa' to one of "
                    << KMethodTypeNames_.toc()
                    << " and 'kappaName' to the name of the volScalar"
                    << " or volSymmTensor field (if kappa=lookup)"
                    << exit(FatalError);
            }

            break;
        }

        default:
        {
            FatalErrorIn
            (
                "temperatureCoupledBase::kappa(const scalarField&) const"
            )
                << "Unimplemented method " << KMethodTypeNames_[method_] << nl
                << "Please set 'kappa' to one of " << KMethodTypeNames_.toc()
                << " and 'kappaName' to the name of the volScalar"
                << " or volSymmTensor field (if kappa=lookup)"
                << exit(FatalError);
        }
    }

    return scalarField(0);
}


void Foam::temperatureCoupledBase::write(Ostream& os) const
{
    os.writeKeyword("kappa") << KMethodTypeNames_[method_]
        << token::END_STATEMENT << nl;
    os.writeKeyword("kappaName") << kappaName_ << token::END_STATEMENT << nl;
}


// ************************************************************************* //
