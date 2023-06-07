/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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

#include "Reaction.H"


// * * * * * * * * * * * * * * * * Static Data * * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::scalar Foam::Reaction<ThermoType>::TlowDefault(0);

template<class ThermoType>
Foam::scalar Foam::Reaction<ThermoType>::ThighDefault(great);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class ThermoType>
void Foam::Reaction<ThermoType>::setThermo
(
    const PtrList<ThermoType>& speciesThermo
)
{
    typename ThermoType::thermoType rhsThermo
    (
        rhs()[0].stoichCoeff
       *speciesThermo[rhs()[0].index].W()
       *speciesThermo[rhs()[0].index]
    );

    for (label i=1; i<rhs().size(); ++i)
    {
        rhsThermo +=
            rhs()[i].stoichCoeff
           *speciesThermo[rhs()[i].index].W()
           *speciesThermo[rhs()[i].index];
    }

    typename ThermoType::thermoType lhsThermo
    (
        lhs()[0].stoichCoeff
       *speciesThermo[lhs()[0].index].W()
       *speciesThermo[lhs()[0].index]
    );

    for (label i=1; i<lhs().size(); ++i)
    {
        lhsThermo +=
            lhs()[i].stoichCoeff
           *speciesThermo[lhs()[i].index].W()
           *speciesThermo[lhs()[i].index];
    }

    // Check for mass imbalance in the reaction. A value of 1 corresponds to an
    // error of 1 H atom in the reaction; i.e. 1 kg/kmol.
    if (mag(lhsThermo.Y() - rhsThermo.Y()) > 0.1)
    {
        FatalErrorInFunction
            << "Mass imbalance for reaction " << name() << ": "
            << mag(lhsThermo.Y() - rhsThermo.Y()) << " kg/kmol"
            << exit(FatalError);
    }

    ThermoType::thermoType::operator=(lhsThermo == rhsThermo);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::Reaction<ThermoType>::Reaction
(
    const speciesTable& species,
    const PtrList<ThermoType>& speciesThermo,
    const List<specieCoeffs>& lhs,
    const List<specieCoeffs>& rhs
)
:
    reaction(species, lhs, rhs),
    ThermoType::thermoType(speciesThermo[0]),
    Tlow_(TlowDefault),
    Thigh_(ThighDefault)
{
    setThermo(speciesThermo);
}


template<class ThermoType>
Foam::Reaction<ThermoType>::Reaction
(
    const Reaction<ThermoType>& r,
    const speciesTable& species
)
:
    reaction(r, species),
    ThermoType::thermoType(r),
    Tlow_(r.Tlow()),
    Thigh_(r.Thigh())
{}


template<class ThermoType>
Foam::Reaction<ThermoType>::Reaction
(
    const speciesTable& species,
    const PtrList<ThermoType>& speciesThermo,
    const dictionary& dict
)
:
    reaction(species, dict),
    ThermoType::thermoType(speciesThermo[0]),
    Tlow_(dict.lookupOrDefault<scalar>("Tlow", TlowDefault)),
    Thigh_(dict.lookupOrDefault<scalar>("Thigh", ThighDefault))
{
    setThermo(speciesThermo);
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::autoPtr<Foam::Reaction<ThermoType>>
Foam::Reaction<ThermoType>::New
(
    const speciesTable& species,
    const PtrList<ThermoType>& speciesThermo,
    const dictionary& dict
)
{
    const word& reactionTypeName = dict.lookup("type");

    typename dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(reactionTypeName);

    // Backwards compatibility check. Reaction names used to have "Reaction"
    // (Reaction<ThermoType>::typeName_()) appended. This was
    // removed as it is unnecessary given the context in which the reaction is
    // specified. If this reaction name was not found, search also for the old
    // name.
    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        cstrIter = dictionaryConstructorTablePtr_->find
        (
            reactionTypeName.removeTrailing(typeName_())
        );
    }

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown reaction type "
            << reactionTypeName << nl << nl
            << "Valid reaction types are :" << nl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<Reaction<ThermoType>>
    (
        cstrIter()(species, speciesThermo, dict)
    );
}


template<class ThermoType>
Foam::autoPtr<Foam::Reaction<ThermoType>>
Foam::Reaction<ThermoType>::New
(
    const speciesTable& species,
    const PtrList<ThermoType>& speciesThermo,
    const objectRegistry& ob,
    const dictionary& dict
)
{
    // If the objectRegistry constructor table is empty
    // use the dictionary constructor table only
    if (!objectRegistryConstructorTablePtr_)
    {
        return New(species, speciesThermo, dict);
    }

    const word& reactionTypeName = dict.lookup("type");

    typename objectRegistryConstructorTable::iterator cstrIter =
        objectRegistryConstructorTablePtr_->find(reactionTypeName);

    // Backwards compatibility check. See above.
    if (cstrIter == objectRegistryConstructorTablePtr_->end())
    {
        cstrIter = objectRegistryConstructorTablePtr_->find
        (
            reactionTypeName.removeTrailing(typeName_())
        );
    }

    if (cstrIter == objectRegistryConstructorTablePtr_->end())
    {
        typename dictionaryConstructorTable::iterator cstrIter =
            dictionaryConstructorTablePtr_->find(reactionTypeName);

        // Backwards compatibility check. See above.
        if (cstrIter == dictionaryConstructorTablePtr_->end())
        {
            cstrIter = dictionaryConstructorTablePtr_->find
            (
                reactionTypeName.removeTrailing(typeName_())
            );
        }

        if (cstrIter == dictionaryConstructorTablePtr_->end())
        {
            FatalErrorInFunction
                << "Unknown reaction type "
                << reactionTypeName << nl << nl
                << "Valid reaction types are :" << nl
                << dictionaryConstructorTablePtr_->sortedToc()
                << objectRegistryConstructorTablePtr_->sortedToc()
                << exit(FatalError);
        }

        return autoPtr<Reaction<ThermoType>>
        (
            cstrIter()(species, speciesThermo, dict)
        );
    }

    return autoPtr<Reaction<ThermoType>>
    (
        cstrIter()(species, speciesThermo, ob, dict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ThermoType>
void Foam::Reaction<ThermoType>::write(Ostream& os) const
{
    reaction::write(os);
}


template<class ThermoType>
void Foam::Reaction<ThermoType>::C
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li,
    scalar& Cf,
    scalar& Cr
) const
{
    Cf = Cr = 1;

    forAll(lhs(), i)
    {
        const label si = lhs()[i].index;
        const specieExponent& el = lhs()[i].exponent;
        Cf *= c[si] >= small || el >= 1 ? pow(max(c[si], 0), el) : 0;
    }

    forAll(rhs(), i)
    {
        const label si = rhs()[i].index;
        const specieExponent& er = rhs()[i].exponent;
        Cr *= c[si] >= small || er >= 1 ? pow(max(c[si], 0), er) : 0;
    }
}


template<class ThermoType>
Foam::scalar Foam::Reaction<ThermoType>::omega
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li,
    scalar& omegaf,
    scalar& omegar
) const
{
    const scalar clippedT = min(max(T, this->Tlow()), this->Thigh());

    // Rate constants
    const scalar kf = this->kf(p, clippedT, c, li);
    const scalar kr = this->kr(kf, p, clippedT, c, li);

    // Concentration products
    scalar Cf, Cr;
    this->C(p, T, c, li, Cf, Cr);

    omegaf = kf*Cf;
    omegar = kr*Cr;

    return omegaf - omegar;
}


template<class ThermoType>
void Foam::Reaction<ThermoType>::dNdtByV
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li,
    scalarField& dNdtByV,
    const bool reduced,
    const List<label>& c2s,
    const label Nsi0
) const
{
    scalar omegaf, omegar;
    const scalar omega = this->omega(p, T, c, li, omegaf, omegar);

    forAll(lhs(), i)
    {
        const label si = reduced ? c2s[lhs()[i].index] : lhs()[i].index;
        const scalar sl = lhs()[i].stoichCoeff;
        dNdtByV[Nsi0 + si] -= sl*omega;
    }
    forAll(rhs(), i)
    {
        const label si = reduced ? c2s[rhs()[i].index] : rhs()[i].index;
        const scalar sr = rhs()[i].stoichCoeff;
        dNdtByV[Nsi0 + si] += sr*omega;
    }
}


template<class ThermoType>
void Foam::Reaction<ThermoType>::ddNdtByVdcTp
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li,
    scalarField& dNdtByV,
    scalarSquareMatrix& ddNdtByVdcTp,
    const bool reduced,
    const List<label>& c2s,
    const label Nsi0,
    const label Tsi,
    scalarField& cTpWork0,
    scalarField& cTpWork1
) const
{
    // Rate constants
    const scalar kf = this->kf(p, T, c, li);
    const scalar kr = this->kr(kf, p, T, c, li);

    // Concentration products
    scalar Cf, Cr;
    this->C(p, T, c, li, Cf, Cr);

    // Overall reaction rate
    const scalar omega = kf*Cf - kr*Cr;

    // Specie reaction rates
    forAll(lhs(), i)
    {
        const label si = reduced ? c2s[lhs()[i].index] : lhs()[i].index;
        const scalar sl = lhs()[i].stoichCoeff;
        dNdtByV[Nsi0 + si] -= sl*omega;
    }
    forAll(rhs(), i)
    {
        const label si = reduced ? c2s[rhs()[i].index] : rhs()[i].index;
        const scalar sr = rhs()[i].stoichCoeff;
        dNdtByV[Nsi0 + si] += sr*omega;
    }

    // Jacobian contributions from the derivative of the concentration products
    // w.r.t. concentration
    {
        forAll(lhs(), j)
        {
            const label sj = reduced ? c2s[lhs()[j].index] : lhs()[j].index;

            scalar dCfdcj = 1;
            forAll(lhs(), i)
            {
                const label si = lhs()[i].index;
                const specieExponent& el = lhs()[i].exponent;
                if (i == j)
                {
                    dCfdcj *=
                        c[si] >= small || el >= 1
                      ? el*pow(max(c[si], 0), el - specieExponent(label(1)))
                      : 0;
                }
                else
                {
                    dCfdcj *=
                        c[si] >= small || el >= 1
                      ? pow(max(c[si], 0), el)
                      : 0;
                }
            }

            forAll(lhs(), i)
            {
                const label si = reduced ? c2s[lhs()[i].index] : lhs()[i].index;
                const scalar sl = lhs()[i].stoichCoeff;
                ddNdtByVdcTp(Nsi0 + si, Nsi0 + sj) -= sl*kf*dCfdcj;
            }
            forAll(rhs(), i)
            {
                const label si = reduced ? c2s[rhs()[i].index] : rhs()[i].index;
                const scalar sr = rhs()[i].stoichCoeff;
                ddNdtByVdcTp(Nsi0 + si, Nsi0 + sj) += sr*kf*dCfdcj;
            }
        }

        forAll(rhs(), j)
        {
            const label sj = reduced ? c2s[rhs()[j].index] : rhs()[j].index;

            scalar dCrcj = 1;
            forAll(rhs(), i)
            {
                const label si = rhs()[i].index;
                const specieExponent& er = rhs()[i].exponent;
                if (i == j)
                {
                    dCrcj *=
                        c[si] >= small || er >= 1
                      ? er*pow(max(c[si], 0), er - specieExponent(label(1)))
                      : 0;
                }
                else
                {
                    dCrcj *=
                        c[si] >= small || er >= 1
                      ? pow(max(c[si], 0), er)
                      : 0;
                }
            }

            forAll(lhs(), i)
            {
                const label si = reduced ? c2s[lhs()[i].index] : lhs()[i].index;
                const scalar sl = lhs()[i].stoichCoeff;
                ddNdtByVdcTp(Nsi0 + si, Nsi0 + sj) += sl*kr*dCrcj;
            }
            forAll(rhs(), i)
            {
                const label si = reduced ? c2s[rhs()[i].index] : rhs()[i].index;
                const scalar sr = rhs()[i].stoichCoeff;
                ddNdtByVdcTp(Nsi0 + si, Nsi0 + sj) -= sr*kr*dCrcj;
            }
        }
    }

    // Jacobian contributions from the derivative of the rate constants
    // w.r.t. temperature
    {
        const scalar dkfdT = this->dkfdT(p, T, c, li);
        const scalar dkrdT = this->dkrdT(p, T, c, li, dkfdT, kr);

        const scalar dwdT = dkfdT*Cf - dkrdT*Cr;
        forAll(lhs(), i)
        {
            const label si = reduced ? c2s[lhs()[i].index] : lhs()[i].index;
            const scalar sl = lhs()[i].stoichCoeff;
            ddNdtByVdcTp(Nsi0 + si, Tsi) -= sl*dwdT;
        }
        forAll(rhs(), i)
        {
            const label si = reduced ? c2s[rhs()[i].index] : rhs()[i].index;
            const scalar sr = rhs()[i].stoichCoeff;
            ddNdtByVdcTp(Nsi0 + si, Tsi) += sr*dwdT;
        }
    }

    // Jacobian contributions from the derivative of the rate constants
    // w.r.t. concentration
    if (hasDkdc())
    {
        scalarField& dkfdc = cTpWork0;
        scalarField& dkrdc = cTpWork1;

        this->dkfdc(p, T, c, li, dkfdc);
        this->dkrdc(p, T, c, li, dkfdc, kr, dkrdc);

        forAll(c, j)
        {
            const label sj = reduced ? c2s[j] : j;

            if (sj == -1) continue;

            const scalar dwdc = dkfdc[j]*Cf - dkrdc[j]*Cr;
            forAll(lhs(), i)
            {
                const label si = reduced ? c2s[lhs()[i].index] : lhs()[i].index;
                const scalar sl = lhs()[i].stoichCoeff;
                ddNdtByVdcTp(Nsi0 + si, Nsi0 + sj) -= sl*dwdc;
            }
            forAll(rhs(), i)
            {
                const label si = reduced ? c2s[rhs()[i].index] : rhs()[i].index;
                const scalar sr = rhs()[i].stoichCoeff;
                ddNdtByVdcTp(Nsi0 + si, Nsi0 + sj) += sr*dwdc;
            }
        }
    }
}


// ************************************************************************* //
