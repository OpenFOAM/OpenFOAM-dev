/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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

template<class ReactionThermo>
Foam::scalar Foam::Reaction<ReactionThermo>::TlowDefault(0);

template<class ReactionThermo>
Foam::scalar Foam::Reaction<ReactionThermo>::ThighDefault(great);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class ReactionThermo>
void Foam::Reaction<ReactionThermo>::setThermo
(
    const HashPtrTable<ReactionThermo>& thermoDatabase
)
{
    typename ReactionThermo::thermoType rhsThermo
    (
        rhs()[0].stoichCoeff
       *(*thermoDatabase[species()[rhs()[0].index]]).W()
       *(*thermoDatabase[species()[rhs()[0].index]])
    );

    for (label i=1; i<rhs().size(); ++i)
    {
        rhsThermo +=
            rhs()[i].stoichCoeff
           *(*thermoDatabase[species()[rhs()[i].index]]).W()
           *(*thermoDatabase[species()[rhs()[i].index]]);
    }

    typename ReactionThermo::thermoType lhsThermo
    (
        lhs()[0].stoichCoeff
       *(*thermoDatabase[species()[lhs()[0].index]]).W()
       *(*thermoDatabase[species()[lhs()[0].index]])
    );

    for (label i=1; i<lhs().size(); ++i)
    {
        lhsThermo +=
            lhs()[i].stoichCoeff
           *(*thermoDatabase[species()[lhs()[i].index]]).W()
           *(*thermoDatabase[species()[lhs()[i].index]]);
    }

    // Check for mass imbalance in the reaction
    // A value of 1 corresponds to an error of 1 H atom in the reaction,
    // i.e. 1 kg/kmol
    if (mag(lhsThermo.Y() - rhsThermo.Y()) > 0.1)
    {
        FatalErrorInFunction
            << "Mass imbalance for reaction " << name() << ": "
            << mag(lhsThermo.Y() - rhsThermo.Y()) << " kg/kmol"
            << exit(FatalError);
    }

    ReactionThermo::thermoType::operator=(lhsThermo == rhsThermo);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ReactionThermo>
Foam::Reaction<ReactionThermo>::Reaction
(
    const speciesTable& species,
    const List<specieCoeffs>& lhs,
    const List<specieCoeffs>& rhs,
    const HashPtrTable<ReactionThermo>& thermoDatabase
)
:
    reaction(species, lhs, rhs),
    ReactionThermo::thermoType(*thermoDatabase[species[0]]),
    Tlow_(TlowDefault),
    Thigh_(ThighDefault)
{
    setThermo(thermoDatabase);
}


template<class ReactionThermo>
Foam::Reaction<ReactionThermo>::Reaction
(
    const Reaction<ReactionThermo>& r,
    const speciesTable& species
)
:
    reaction(r, species),
    ReactionThermo::thermoType(r),
    Tlow_(r.Tlow()),
    Thigh_(r.Thigh())
{}


template<class ReactionThermo>
Foam::Reaction<ReactionThermo>::Reaction
(
    const speciesTable& species,
    const HashPtrTable<ReactionThermo>& thermoDatabase,
    const dictionary& dict
)
:
    reaction(species, dict),
    ReactionThermo::thermoType(*thermoDatabase[species[0]]),
    Tlow_(dict.lookupOrDefault<scalar>("Tlow", TlowDefault)),
    Thigh_(dict.lookupOrDefault<scalar>("Thigh", ThighDefault))
{
    setThermo(thermoDatabase);
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

template<class ReactionThermo>
Foam::autoPtr<Foam::Reaction<ReactionThermo>>
Foam::Reaction<ReactionThermo>::New
(
    const speciesTable& species,
    const HashPtrTable<ReactionThermo>& thermoDatabase,
    const dictionary& dict
)
{
    const word& reactionTypeName = dict.lookup("type");

    typename dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(reactionTypeName);

    // Backwards compatibility check. Reaction names used to have "Reaction"
    // (Reaction<ReactionThermo>::typeName_()) appended. This was removed as it
    // is unnecessary given the context in which the reaction is specified. If
    // this reaction name was not found, search also for the old name.
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

    return autoPtr<Reaction<ReactionThermo>>
    (
        cstrIter()(species, thermoDatabase, dict)
    );
}


template<class ReactionThermo>
Foam::autoPtr<Foam::Reaction<ReactionThermo>>
Foam::Reaction<ReactionThermo>::New
(
    const speciesTable& species,
    const HashPtrTable<ReactionThermo>& thermoDatabase,
    const objectRegistry& ob,
    const dictionary& dict
)
{
    // If the objectRegistry constructor table is empty
    // use the dictionary constructor table only
    if (!objectRegistryConstructorTablePtr_)
    {
        return New(species, thermoDatabase, dict);
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

        return autoPtr<Reaction<ReactionThermo>>
        (
            cstrIter()(species, thermoDatabase, dict)
        );
    }

    return autoPtr<Reaction<ReactionThermo>>
    (
        cstrIter()(species, thermoDatabase, ob, dict)
    );
}


template<class ReactionThermo>
Foam::autoPtr<Foam::Reaction<ReactionThermo>>
Foam::Reaction<ReactionThermo>::New
(
    const speciesTable& species,
    const PtrList<ReactionThermo>& speciesThermo,
    const dictionary& dict
)
{
    HashPtrTable<ReactionThermo> thermoDatabase;
    forAll(speciesThermo, i)
    {
        thermoDatabase.insert
        (
            speciesThermo[i].name(),
            speciesThermo[i].clone().ptr()
        );
    }

    return New(species, thermoDatabase, dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ReactionThermo>
void Foam::Reaction<ReactionThermo>::write(Ostream& os) const
{
    reaction::write(os);
}


template<class ReactionThermo>
void Foam::Reaction<ReactionThermo>::omega
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li,
    scalarField& dcdt
) const
{
    scalar omegaf, omegar;

    const scalar omegaI = omega(p, T, c, li, omegaf, omegar);

    forAll(lhs(), i)
    {
        const label si = lhs()[i].index;
        const scalar sl = lhs()[i].stoichCoeff;
        dcdt[si] -= sl*omegaI;
    }
    forAll(rhs(), i)
    {
        const label si = rhs()[i].index;
        const scalar sr = rhs()[i].stoichCoeff;
        dcdt[si] += sr*omegaI;
    }
}


template<class ReactionThermo>
Foam::scalar Foam::Reaction<ReactionThermo>::omega
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li,
    scalar& omegaf,
    scalar& omegar
) const
{
    scalar clippedT = min(max(T, this->Tlow()), this->Thigh());

    omegaf = this->kf(p, clippedT, c, li);
    omegar = this->kr(omegaf, p, clippedT, c, li);

    forAll(lhs(), i)
    {
        const label si = lhs()[i].index;
        const scalar el = lhs()[i].exponent;
        omegaf *= c[si] >= small || el >= 1 ? pow(max(c[si], 0), el) : 0;
    }
    forAll(rhs(), i)
    {
        const label si = rhs()[i].index;
        const scalar er = rhs()[i].exponent;
        omegar *= c[si] >= small || er >= 1 ? pow(max(c[si], 0), er) : 0;
    }

    return omegaf - omegar;
}


template<class ReactionThermo>
void Foam::Reaction<ReactionThermo>::dwdc
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li,
    scalarSquareMatrix& J,
    scalarField& dcdt,
    scalar& omegaI,
    scalar& kfwd,
    scalar& kbwd,
    const bool reduced,
    const List<label>& c2s
) const
{
    scalar omegaf, omegar;

    omegaI = omega(p, T, c, li, omegaf, omegar);

    forAll(lhs(), i)
    {
        const label si = reduced ? c2s[lhs()[i].index] : lhs()[i].index;
        const scalar sl = lhs()[i].stoichCoeff;
        dcdt[si] -= sl*omegaI;
    }
    forAll(rhs(), i)
    {
        const label si = reduced ? c2s[rhs()[i].index] : rhs()[i].index;
        const scalar sr = rhs()[i].stoichCoeff;
        dcdt[si] += sr*omegaI;
    }

    kfwd = this->kf(p, T, c, li);
    kbwd = this->kr(kfwd, p, T, c, li);

    forAll(lhs(), j)
    {
        const label sj = reduced ? c2s[lhs()[j].index] : lhs()[j].index;
        scalar kf = kfwd;
        forAll(lhs(), i)
        {
            const label si = lhs()[i].index;
            const scalar el = lhs()[i].exponent;
            if (i == j)
            {
                kf *=
                    c[si] >= small || el >= 1
                  ? el*pow(max(c[si], 0), el - 1)
                  : 0;
            }
            else
            {
                kf *=
                    c[si] >= small || el >= 1
                  ? pow(max(c[si], 0), el)
                  : 0;
            }
        }

        forAll(lhs(), i)
        {
            const label si = reduced ? c2s[lhs()[i].index] : lhs()[i].index;
            const scalar sl = lhs()[i].stoichCoeff;
            J(si, sj) -= sl*kf;
        }
        forAll(rhs(), i)
        {
            const label si = reduced ? c2s[rhs()[i].index] : rhs()[i].index;
            const scalar sr = rhs()[i].stoichCoeff;
            J(si, sj) += sr*kf;
        }
    }

    forAll(rhs(), j)
    {
        const label sj = reduced ? c2s[rhs()[j].index] : rhs()[j].index;
        scalar kr = kbwd;
        forAll(rhs(), i)
        {
            const label si = rhs()[i].index;
            const scalar er = rhs()[i].exponent;
            if (i == j)
            {
                kr *=
                    c[si] >= small || er >= 1
                  ? er*pow(max(c[si], 0), er - 1)
                  : 0;
            }
            else
            {
                kr *=
                    c[si] >= small || er >= 1
                  ? pow(max(c[si], 0), er)
                  : 0;
            }
        }

        forAll(lhs(), i)
        {
            const label si = reduced ? c2s[lhs()[i].index] : lhs()[i].index;
            const scalar sl = lhs()[i].stoichCoeff;
            J(si, sj) += sl*kr;
        }
        forAll(rhs(), i)
        {
            const label si = reduced ? c2s[rhs()[i].index] : rhs()[i].index;
            const scalar sr = rhs()[i].stoichCoeff;
            J(si, sj) -= sr*kr;
        }
    }

    // When third-body species are involved, additional terms are added
    // beta function returns an empty list when third-body are not involved
    const List<Tuple2<label, scalar>>& beta = this->beta();
    if (notNull(beta))
    {
        // This temporary array needs to be cached for efficiency
        scalarField dcidc(beta.size());
        this->dcidc(p, T, c, li, dcidc);

        forAll(beta, j)
        {
            label sj = beta[j].first();
            sj = reduced ? c2s[sj] : sj;
            if (sj != -1)
            {
                forAll(lhs(), i)
                {
                    const label si =
                        reduced ? c2s[lhs()[i].index] : lhs()[i].index;
                    const scalar sl = lhs()[i].stoichCoeff;
                    J(si, sj) -= sl*dcidc[j]*omegaI;
                }
                forAll(rhs(), i)
                {
                    const label si =
                        reduced ? c2s[rhs()[i].index] : rhs()[i].index;
                    const scalar sr = rhs()[i].stoichCoeff;
                    J(si, sj) += sr*dcidc[j]*omegaI;
                }
            }
        }
    }
}


template<class ReactionThermo>
void Foam::Reaction<ReactionThermo>::dwdT
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li,
    const scalar omegaI,
    const scalar kfwd,
    const scalar kbwd,
    scalarSquareMatrix& J,
    const bool reduced,
    const List<label>& c2s,
    const label indexT
) const
{
    scalar dkfdT = this->dkfdT(p, T, c, li);
    scalar dkrdT = this->dkrdT(p, T, c, li, dkfdT, kbwd);

    forAll(lhs(), i)
    {
        const label si = lhs()[i].index;
        const scalar el = lhs()[i].exponent;
        dkfdT *= c[si] >= small || el >= 1 ? pow(max(c[si], 0), el) : 0;
    }
    forAll(rhs(), i)
    {
        const label si = rhs()[i].index;
        const scalar er = rhs()[i].exponent;
        dkrdT *= c[si] >= small || er >= 1 ? pow(max(c[si], 0), er) : 0;
    }

    const scalar dqidT = dkfdT - dkrdT;

    // For reactions including third-body efficiencies or pressure dependent
    // reaction, an additional term is needed
    const scalar dcidT = omegaI*this->dcidT(p, T, c, li);

    // J(i, indexT) = sum_reactions nu_i dqdT
    forAll(lhs(), i)
    {
        const label si = reduced ? c2s[lhs()[i].index] : lhs()[i].index;
        const scalar sl = lhs()[i].stoichCoeff;
        J(si, indexT) -= sl*(dqidT + dcidT);
    }
    forAll(rhs(), i)
    {
        const label si = reduced ? c2s[rhs()[i].index] : rhs()[i].index;
        const scalar sr = rhs()[i].stoichCoeff;
        J(si, indexT) += sr*(dqidT + dcidT);
    }
}


// ************************************************************************* //
