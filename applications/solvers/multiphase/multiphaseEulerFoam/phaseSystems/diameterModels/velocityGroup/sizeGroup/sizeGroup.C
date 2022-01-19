/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2022 OpenFOAM Foundation
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

#include "sizeGroup.H"
#include "mixedFvPatchField.H"
#include "shapeModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

using Foam::constant::mathematical::pi;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diameterModels::sizeGroup::sizeGroup
(
    const word& name,
    const dictionary& dict,
    const phaseModel& phase,
    const velocityGroup& velocityGroup,
    const fvMesh& mesh
)
:
    volScalarField
    (
        IOobject
        (
            IOobject::groupName
            (
                name,
                velocityGroup.phase().name()
            ),
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(name, dimless, dict.lookup<scalar>("value")),
        velocityGroup.f().boundaryField().types()
    ),
    dict_(dict),
    phase_(phase),
    velocityGroup_(velocityGroup),
    dSph_("dSph", dimLength, dict),
    x_("x", pi/6*pow3(dSph_)),
    value_(dict.lookup<scalar>("value"))
{
    // Adjust refValue at mixedFvPatchField boundaries
    forAll(this->boundaryField(), patchi)
    {
        typedef mixedFvPatchField<scalar> mixedFvPatchScalarField;

        if
        (
            isA<mixedFvPatchScalarField>(this->boundaryFieldRef()[patchi])
        )
        {
            mixedFvPatchScalarField& f =
                refCast<mixedFvPatchScalarField>
                (
                    this->boundaryFieldRef()[patchi]
                );

            f.refValue() = value_;
        }
    }

    shapeModel_ = shapeModel::New(velocityGroup_.diameterProperties(), *this);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::diameterModels::sizeGroup::~sizeGroup()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::autoPtr<Foam::diameterModels::sizeGroup>
Foam::diameterModels::sizeGroup::clone() const
{
    notImplemented("sizeGroup::clone() const");
    return autoPtr<sizeGroup>(nullptr);
}


const Foam::autoPtr<Foam::label>& Foam::diameterModels::sizeGroup::i() const
{
    if (!i_.valid())
    {
        const populationBalanceModel& popBal =
            this->mesh().lookupObject<populationBalanceModel>
            (
                velocityGroup_.popBalName()
            );

        forAll(popBal.sizeGroups(), j)
        {
            if (&popBal.sizeGroups()[j] == &*this)
            {
                i_.set(new label(j));
            }
        }
    }

    return i_;
}


const Foam::tmp<Foam::volScalarField>
Foam::diameterModels::sizeGroup::a() const
{
    return shapeModel_->a();
}


const Foam::tmp<Foam::volScalarField>
Foam::diameterModels::sizeGroup::d() const
{
    return shapeModel_->d();
}


void Foam::diameterModels::sizeGroup::correct()
{
    shapeModel_->correct();
}


// ************************************************************************* //
