/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2023 OpenFOAM Foundation
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

using Foam::constant::mathematical::pi;

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::IOobject Foam::diameterModels::sizeGroup::fieldIo
(
    const word& name,
    const label i,
    const velocityGroup& group,
    const IOobject::readOption r,
    const bool registerObject
)
{
    return
        IOobject
        (
            IOobject::groupName
            (
                name + (i == -1 ? "Default" : Foam::name(i)),
                group.phase().name()
            ),
            group.phase().mesh().time().name(),
            group.phase().mesh(),
            r,
            IOobject::AUTO_WRITE,
            registerObject
        );
}


Foam::tmp<Foam::volScalarField> Foam::diameterModels::sizeGroup::field
(
    const word& name,
    const label i,
    const velocityGroup& group
)
{
    typeIOobject<volScalarField> io
    (
        fieldIo(name, i, group, IOobject::MUST_READ, false)
    );

    return
        tmp<volScalarField>
        (
            new volScalarField
            (
                io.headerOk()
              ? io
              : fieldIo(name, -1, group, IOobject::MUST_READ, false),
                group.phase().mesh()
            )
        );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diameterModels::sizeGroup::sizeGroup
(
    const label i,
    const dictionary& dict,
    const velocityGroup& group
)
:
    volScalarField(fieldIo("f", i, group), field("f", i, group)),
    i_(i),
    group_(group),
    dSph_("dSph", dimLength, dict),
    x_("x", pi/6*pow3(dSph_)),
    shapeModel_(shapeModel::New(group_.diameterProperties(), *this, dict))
{
    // Check and filter for old syntax (remove in due course)
    if (dict.found("value"))
    {
        FatalErrorInFunction
            << "A 'value' entry should not be specified for size-group #"
            << i << " of population balance "
            << group.popBalName()
            << ". Instead, the value should be initialised within the field, "
            << this->name() << " (or the default field, "
            << IOobject::groupName("fDefault", group.phase().name())
            << ", as appropriate)."
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::diameterModels::sizeGroup::~sizeGroup()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::autoPtr<Foam::diameterModels::sizeGroup>
Foam::diameterModels::sizeGroup::clone() const
{
    NotImplemented;
    return autoPtr<sizeGroup>(nullptr);
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
