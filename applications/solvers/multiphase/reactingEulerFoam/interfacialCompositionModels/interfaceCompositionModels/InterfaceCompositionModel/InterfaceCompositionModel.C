/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2018 OpenFOAM Foundation
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

#include "InterfaceCompositionModel.H"
#include "phaseModel.H"
#include "phasePair.H"
#include "pureMixture.H"
#include "multiComponentMixture.H"
#include "rhoThermo.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Thermo, class OtherThermo>
template<class ThermoType>
const typename Foam::multiComponentMixture<ThermoType>::thermoType&
Foam::InterfaceCompositionModel<Thermo, OtherThermo>::getLocalThermo
(
    const word& speciesName,
    const multiComponentMixture<ThermoType>& globalThermo
) const
{
    return
        globalThermo.getLocalThermo
        (
            globalThermo.species()
            [
                speciesName
            ]
        );
}


template<class Thermo, class OtherThermo>
template<class ThermoType>
const typename Foam::pureMixture<ThermoType>::thermoType&
Foam::InterfaceCompositionModel<Thermo, OtherThermo>::getLocalThermo
(
    const word& speciesName,
    const pureMixture<ThermoType>& globalThermo
) const
{
    return globalThermo.cellMixture(0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo, class OtherThermo>
Foam::InterfaceCompositionModel<Thermo, OtherThermo>::InterfaceCompositionModel
(
    const dictionary& dict,
    const phasePair& pair
)
:
    interfaceCompositionModel(dict, pair),
    thermo_
    (
        pair.phase1().mesh().lookupObject<Thermo>
        (
            IOobject::groupName(basicThermo::dictName, pair.phase1().name())
        )
    ),
    otherThermo_
    (
        pair.phase2().mesh().lookupObject<OtherThermo>
        (
            IOobject::groupName(basicThermo::dictName, pair.phase2().name())
        )
    ),
    Le_("Le", dimless, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Thermo, class OtherThermo>
Foam::InterfaceCompositionModel<Thermo, OtherThermo>::
~InterfaceCompositionModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Thermo, class OtherThermo>
Foam::tmp<Foam::volScalarField>
Foam::InterfaceCompositionModel<Thermo, OtherThermo>::dY
(
    const word& speciesName,
    const volScalarField& Tf
) const
{
    return
        Yf(speciesName, Tf)
      - thermo_.composition().Y()
        [
            thermo_.composition().species()[speciesName]
        ];
}


template<class Thermo, class OtherThermo>
Foam::tmp<Foam::volScalarField>
Foam::InterfaceCompositionModel<Thermo, OtherThermo>::D
(
    const word& speciesName
) const
{
    const typename Thermo::thermoType& localThermo =
        getLocalThermo
        (
            speciesName,
            thermo_
        );

    const volScalarField& p(thermo_.p());

    const volScalarField& T(thermo_.T());

    tmp<volScalarField> tmpD
    (
        new volScalarField
        (
            IOobject
            (
                IOobject::groupName("D", pair_.name()),
                p.time().timeName(),
                p.mesh()
            ),
            p.mesh(),
            dimensionedScalar("zero", dimArea/dimTime, 0)
        )
    );

    volScalarField& D(tmpD.ref());

    forAll(p, celli)
    {
        D[celli] =
            localThermo.alphah(p[celli], T[celli])
           /localThermo.rho(p[celli], T[celli]);
    }

    D /= Le_;

    return tmpD;
}


template<class Thermo, class OtherThermo>
Foam::tmp<Foam::volScalarField>
Foam::InterfaceCompositionModel<Thermo, OtherThermo>::L
(
    const word& speciesName,
    const volScalarField& Tf
) const
{
    const typename Thermo::thermoType& localThermo =
        getLocalThermo
        (
            speciesName,
            thermo_
        );
    const typename OtherThermo::thermoType& otherLocalThermo =
        getLocalThermo
        (
            speciesName,
            otherThermo_
        );

    const volScalarField& p(thermo_.p());
    const volScalarField& otherP(otherThermo_.p());

    tmp<volScalarField> tmpL
    (
        new volScalarField
        (
            IOobject
            (
                IOobject::groupName("L", pair_.name()),
                p.time().timeName(),
                p.mesh()
            ),
            p.mesh(),
            dimensionedScalar("zero", dimEnergy/dimMass, 0)
        )
    );

    volScalarField& L(tmpL.ref());

    forAll(p, celli)
    {
        L[celli] =
            localThermo.Ha(p[celli], Tf[celli])
          - otherLocalThermo.Ha(otherP[celli], Tf[celli]);
    }

    return tmpL;
}


template<class Thermo, class OtherThermo>
void Foam::InterfaceCompositionModel<Thermo, OtherThermo>::addMDotL
(
    const volScalarField& K,
    const volScalarField& Tf,
    volScalarField& mDotL,
    volScalarField& mDotLPrime
) const
{
    forAllConstIter(hashedWordList, this->speciesNames_, iter)
    {
        volScalarField rhoKDL
        (
            thermo_.rhoThermo::rho()
           *K
           *D(*iter)
           *L(*iter, Tf)
        );

        mDotL += rhoKDL*dY(*iter, Tf);
        mDotLPrime += rhoKDL*YfPrime(*iter, Tf);
    }
}


// ************************************************************************* //
