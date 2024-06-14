/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2024 OpenFOAM Foundation
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
    const bool allowNone
)
{
    if (dict.found(name) || allowNone)
    {
        const token t(dict.lookup(name));

        if (allowNone && t.isWord() && t.wordToken() == "none")
        {
            return {false, NaN};
        }

        if (t.isNumber())
        {
            forAll(bounds, i)
            {
                const label s = i == 0 ? -1 : +1;

                if (s*t.number() > s*bounds[i])
                {
                    FatalErrorInFunction
                        << "Blending parameter " << name << " is "
                        << (i == 0 ? "less" : "greater") << " than "
                        << bounds[i] << exit(FatalError);
                }
            }

            return {true, t.number()};
        }

        FatalIOErrorInFunction(dict)
            << "wrong token type - expected Scalar or the word 'none', found "
            << t.info() << exit(FatalIOError);
    }

    return {false, NaN};
}


Foam::Pair<Foam::blendingParameter> Foam::blendingMethod::readParameters
(
    const word& name,
    const dictionary& dict,
    const phaseInterface& interface,
    const Pair<scalar>& bounds,
    const bool allowNone
)
{
    const word name1 = IOobject::groupName(name, interface.phase1().name());
    const word name2 = IOobject::groupName(name, interface.phase2().name());
    return
        Pair<blendingParameter>
        (
            readParameter(name1, dict, bounds, allowNone),
            readParameter(name2, dict, bounds, allowNone)
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
    const label set,
    const bool protect
) const
{
    tmp<volScalarField> talpha = constant(alphas, 0);

    forAllConstIter(phaseInterface, interface_, iter)
    {
        if (0b01 << iter.index() & set)
        {
            talpha.ref() +=
                protect
              ? max(iter().residualAlpha(), alphas[iter().index()])()
              : alphas[iter().index()];
        }
    }

    return talpha;
}


Foam::tmp<Foam::volScalarField> Foam::blendingMethod::parameter
(
    const UPtrList<const volScalarField>& alphas,
    const label set,
    const Pair<blendingParameter>& parameters
) const
{
    tmp<volScalarField> talphaParameter = constant(alphas, 0);

    forAllConstIter(phaseInterface, interface_, iter)
    {
        if (0b01 << iter.index() & set)
        {
            talphaParameter.ref() +=
                max(iter().residualAlpha(), alphas[iter().index()])
               *parameters[iter.index()].value;
        }
    }

    return talphaParameter/alpha(alphas, set, true);
}


Foam::tmp<Foam::volScalarField> Foam::blendingMethod::x
(
    const UPtrList<const volScalarField>& alphas,
    const label phaseSet,
    const label systemSet
) const
{
    return
        systemSet == 0b00
      ? alpha(alphas, phaseSet, false)
      : alpha(alphas, phaseSet, true)/alpha(alphas, systemSet, true);
}


Foam::tmp<Foam::volScalarField> Foam::blendingMethod::f
(
    const UPtrList<const volScalarField>& alphas,
    const label phaseSet,
    const label systemSet
) const
{
    label canBeContinuousPhaseSet = 0b00;
    label canBeContinuousSystemSet = 0b00;
    forAllConstIter(phaseInterface, interface_, iter)
    {
        if (canBeContinuous(iter.index()))
        {
            canBeContinuousPhaseSet += 0b01 << iter.index() & phaseSet;
            canBeContinuousSystemSet += 0b01 << iter.index() & systemSet;
        }
    }

    if (canBeContinuousPhaseSet == 0)
    {
        return constant(alphas, 0);
    }

    if (canBeContinuousPhaseSet == canBeContinuousSystemSet)
    {
        return constant(alphas, 1);
    }

    return
        fContinuous
        (
            alphas,
            canBeContinuousPhaseSet,
            canBeContinuousSystemSet
        );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::blendingMethod::blendingMethod
(
    const dictionary& dict,
    const phaseInterface& interface
)
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
    return f(alphas, 0b10, 0b11);
}


Foam::tmp<Foam::volScalarField> Foam::blendingMethod::f2DispersedIn1
(
    const UPtrList<const volScalarField>& alphas
) const
{
    return f(alphas, 0b01, 0b11);
}


Foam::tmp<Foam::volScalarField> Foam::blendingMethod::fDisplaced
(
    const UPtrList<const volScalarField>& alphas
) const
{
    return 1 - f(alphas, 0b11, 0b00);
}


// ************************************************************************* //
