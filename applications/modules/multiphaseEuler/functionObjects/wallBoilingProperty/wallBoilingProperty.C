/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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

#include "wallBoiling.H"
#include "wallBoilingProperty.H"
#include "wallBoilingPhaseChangeRateFvPatchScalarField.H"
#include "fvModels.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(wallBoilingProperty, 0);
    addToRunTimeSelectionTable(functionObject, wallBoilingProperty, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

const Foam::fv::wallBoiling&
Foam::functionObjects::wallBoilingProperty::model() const
{
    const Foam::fvModels& fvModels = Foam::fvModels::New(this->mesh());

    if (modelName_ != word::null)
    {
        return refCast<const fv::wallBoiling>(fvModels[modelName_]);
    }

    label wallBoilingFvModeli = -1;

    forAll(fvModels, fvModeli)
    {
        if (!isA<fv::wallBoiling>(fvModels[fvModeli])) continue;

        if (wallBoilingFvModeli != -1)
        {
            WarningInFunction
                << "Multiple wall boiling fvModels found for " << typeName
                << " function " << name() << endl;

            return NullObjectRef<fv::wallBoiling>();
        }

        wallBoilingFvModeli = fvModeli;
    }

    if (wallBoilingFvModeli == -1)
    {
        WarningInFunction
            << "Wall boiling fvModel not found for " << typeName
            << " function " << name() << endl;

        return NullObjectRef<fv::wallBoiling>();
    }

    return refCast<const fv::wallBoiling>(fvModels[wallBoilingFvModeli]);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::wallBoilingProperty::wallBoilingProperty
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::wallBoilingProperty::~wallBoilingProperty()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::wallBoilingProperty::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    modelName_ = dict.lookupOrDefault<word>("model", word::null);

    fieldName_ = dict.lookup<word>("field");

    return true;
}


bool Foam::functionObjects::wallBoilingProperty::execute()
{
    const fv::wallBoiling& model = this->model();

    if (isNull(model)) return false;

    tmp<volScalarField> tpropertyVf
    (
        new volScalarField
        (
            IOobject
            (
                model.name() + ":" + fieldName_,
                mesh().time().name(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh(),
            dimensionedScalar
            (
                wallBoilingPhaseChangeRateFvPatchScalarField::propertyDimensions
                (
                    fieldName_
                ),
                0
            )
        )
    );
    volScalarField& propertyVf = tpropertyVf.ref();

    forAll(propertyVf.boundaryField(), patchi)
    {
        if (!model.isBoiling(patchi)) continue;

        const wallBoilingPhaseChangeRateFvPatchScalarField& mDot =
            model.mDotPf(patchi);

        propertyVf.boundaryFieldRef()[patchi] = mDot.property(fieldName_);
    }

    store(tpropertyVf);

    return true;
}


bool Foam::functionObjects::wallBoilingProperty::write()
{
    const fv::wallBoiling& model = this->model();

    if (isNull(model)) return false;

    return writeObject(model.name() + ":" + fieldName_);
}


// ************************************************************************* //
