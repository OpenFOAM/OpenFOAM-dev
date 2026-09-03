/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2026 OpenFOAM Foundation
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

#include "blendingMethod.H"
#include "phaseSystem.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(blendingMethod, 0);
    defineRunTimeSelectionTable(blendingMethod, dictionary);
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::blendingParameter Foam::blendingMethod::readParameter
(
    const word& name,
    const dictionary& dict,
    const Pair<scalar>& bounds,
    const bool allowNone,
    const scalar noneValue
)
{
    if (allowNone)
    {
        ITstream& is = dict.lookup(name);

        const token t(is);

        if (t.isWord() && t.wordToken() == "none")
        {
            return {false, noneValue};
        }

        if (!t.isNumber())
        {
            FatalIOErrorInFunction(is)
                << "wrong token type - expected Scalar or the word 'none', "
                << "found " << t.info() << exit(FatalIOError);
        }
    }

    const scalar value = dict.lookup<scalar>(name);

    forAll(bounds, i)
    {
        const label s = i == 0 ? -1 : +1;

        if (s*value > s*bounds[i])
        {
            FatalErrorInFunction
                << "Blending parameter " << name << " is "
                << (i == 0 ? "less" : "greater") << " than "
                << bounds[i] << exit(FatalError);
        }
    }

    return {true, value};
}


Foam::blendingParameter Foam::blendingMethod::readParameter
(
    const word& name,
    const dictionary& dict,
    const Pair<scalar>& bounds
)
{
    return readParameter(name, dict, bounds, false, NaN);
}


Foam::blendingParameter Foam::blendingMethod::readParameter
(
    const word& name,
    const dictionary& dict,
    const Pair<scalar>& bounds,
    const scalar noneValue
)
{
    return readParameter(name, dict, bounds, true, noneValue);
}


Foam::Pair<Foam::blendingParameter> Foam::blendingMethod::readParameters
(
    const word& name,
    const dictionary& dict,
    const phaseInterface& interface,
    const Pair<scalar>& bounds
)
{
    const word name1 = IOobject::groupName(name, interface.phase1().name());
    const word name2 = IOobject::groupName(name, interface.phase2().name());

    return
        Pair<blendingParameter>
        (
            readParameter(name1, dict, bounds),
            readParameter(name2, dict, bounds)
        );
}


Foam::Pair<Foam::blendingParameter> Foam::blendingMethod::readParameters
(
    const word& name,
    const dictionary& dict,
    const phaseInterface& interface,
    const Pair<scalar>& bounds,
    const scalar noneValue
)
{
    const word name1 = IOobject::groupName(name, interface.phase1().name());
    const word name2 = IOobject::groupName(name, interface.phase2().name());

    return
        Pair<blendingParameter>
        (
            readParameter(name1, dict, bounds, noneValue),
            readParameter(name2, dict, bounds, noneValue)
        );
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::blendingMethod::constant
(
    const UPtrList<const volScalarField>& alphas,
    const scalar k
) const
{
    return
        volScalarField::New
        (
            name(k),
            alphas.first().mesh(),
            dimensionedScalar(dimless, k)
        );
}


Foam::tmp<Foam::volScalarField> Foam::blendingMethod::alpha
(
    const UPtrList<const volScalarField>& alphas,
    const label index,
    const bool protect
) const
{
    switch (index)
    {
        case 0:
        case 1:
        {
            const label phasei = interface_[index].index();
            return
                protect
              ? eval(max(alphas[phasei], interface_[index].residualAlpha()))
              : tmp<volScalarField>(alphas[phasei]);
        }
        case -1:
        {
            const label phasei0 = interface_[0].index();
            const label phasei1 = interface_[1].index();
            return
                protect
              ? eval
                (
                    max(alphas[phasei0], interface_[0].residualAlpha())
                  + max(alphas[phasei1], interface_[1].residualAlpha())
                )
              : eval(alphas[phasei0] + alphas[phasei1]);
        }
    }

    FatalErrorInFunction
        << "Index should be 0, 1, or -1"
        << exit(FatalError);

    return tmp<volScalarField>();
}


Foam::tmp<Foam::volScalarField> Foam::blendingMethod::parameter
(
    const UPtrList<const volScalarField>& alphas,
    const label index,
    const Pair<blendingParameter>& parameters
) const
{
    switch (index)
    {
        case 0:
        case 1:
        {
            return constant(alphas, parameters[index].value);
        }
        case -1:
        {
            const label phasei0 = interface_[0].index();
            const label phasei1 = interface_[1].index();
            return
                (
                    max(alphas[phasei0], interface_[0].residualAlpha())
                   *parameters[0].value
                  + max(alphas[phasei1], interface_[1].residualAlpha())
                   *parameters[1].value
                )
               /(
                   max(alphas[phasei0], interface_[0].residualAlpha())
                 + max(alphas[phasei1], interface_[1].residualAlpha())
                );

        }
    }

    FatalErrorInFunction
        << "Index should be 0, 1, or -1"
        << exit(FatalError);

    return tmp<volScalarField>();
}


Foam::tmp<Foam::volScalarField> Foam::blendingMethod::x
(
    const UPtrList<const volScalarField>& alphas,
    const label index
) const
{
    return
        index == -1
      ? alpha(alphas, -1, false)
      : alpha(alphas, index, true)/alpha(alphas, -1, true);
}


Foam::tmp<Foam::volScalarField> Foam::blendingMethod::fContinuousFiltered
(
    const UPtrList<const volScalarField>& alphas,
    const label index
) const
{
    if (!canBeContinuous(index)) return constant(alphas, 0);

    if (!canBeContinuous(1 - index)) return constant(alphas, 1);

    return fContinuous(alphas, index);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::blendingMethod::blendingMethod(const phaseInterface& interface)
:
    interface_(interface)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::blendingMethod::~blendingMethod()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::blendingMethod::f1DispersedIn2
(
    const UPtrList<const volScalarField>& alphas
) const
{
    return fContinuousFiltered(alphas, 1);
}


Foam::tmp<Foam::volScalarField> Foam::blendingMethod::f2DispersedIn1
(
    const UPtrList<const volScalarField>& alphas
) const
{
    return fContinuousFiltered(alphas, 0);
}


Foam::tmp<Foam::volScalarField> Foam::blendingMethod::f12DisplacedBy3
(
    const UPtrList<const volScalarField>& alphas,
    const label displacingPhasei
) const
{
    return fDisplaced(alphas, displacingPhasei);
}


// ************************************************************************* //
