/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
#include "DynamicList.H"

// * * * * * * * * * * * * * * * * Static Data * * * * * * * * * * * * * * * //

template<class ReactionThermo>
Foam::label Foam::Reaction<ReactionThermo>::nUnNamedReactions(0);

template<class ReactionThermo>
Foam::scalar Foam::Reaction<ReactionThermo>::TlowDefault(0);

template<class ReactionThermo>
Foam::scalar Foam::Reaction<ReactionThermo>::ThighDefault(great);

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class ReactionThermo>
void Foam::Reaction<ReactionThermo>::reactionStrLeft
(
    OStringStream& reaction
) const
{
    for (label i = 0; i < lhs_.size(); ++i)
    {
        if (i > 0)
        {
            reaction << " + ";
        }
        if (mag(lhs_[i].stoichCoeff - 1) > small)
        {
            reaction << lhs_[i].stoichCoeff;
        }
        reaction << species_[lhs_[i].index];
        if (mag(lhs_[i].exponent - lhs_[i].stoichCoeff) > small)
        {
            reaction << "^" << lhs_[i].exponent;
        }
    }
}


template<class ReactionThermo>
void Foam::Reaction<ReactionThermo>::reactionStrRight
(
    OStringStream& reaction
) const
{
    for (label i = 0; i < rhs_.size(); ++i)
    {
        if (i > 0)
        {
            reaction << " + ";
        }
        if (mag(rhs_[i].stoichCoeff - 1) > small)
        {
            reaction << rhs_[i].stoichCoeff;
        }
        reaction << species_[rhs_[i].index];
        if (mag(rhs_[i].exponent - rhs_[i].stoichCoeff) > small)
        {
            reaction << "^" << rhs_[i].exponent;
        }
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class ReactionThermo>
Foam::label Foam::Reaction<ReactionThermo>::getNewReactionID()
{
    return nUnNamedReactions++;
}


template<class ReactionThermo>
Foam::string Foam::Reaction<ReactionThermo>::reactionStr
(
    OStringStream& reaction
) const
{
    reactionStrLeft(reaction);
    reaction << " = ";
    reactionStrRight(reaction);
    return reaction.str();
}


template<class ReactionThermo>
void Foam::Reaction<ReactionThermo>::setThermo
(
    const HashPtrTable<ReactionThermo>& thermoDatabase
)
{
    typename ReactionThermo::thermoType rhsThermo
    (
        rhs_[0].stoichCoeff
       *(*thermoDatabase[species_[rhs_[0].index]]).W()
       *(*thermoDatabase[species_[rhs_[0].index]])
    );

    for (label i=1; i<rhs_.size(); ++i)
    {
        rhsThermo +=
            rhs_[i].stoichCoeff
           *(*thermoDatabase[species_[rhs_[i].index]]).W()
           *(*thermoDatabase[species_[rhs_[i].index]]);
    }

    typename ReactionThermo::thermoType lhsThermo
    (
        lhs_[0].stoichCoeff
       *(*thermoDatabase[species_[lhs_[0].index]]).W()
       *(*thermoDatabase[species_[lhs_[0].index]])
    );

    for (label i=1; i<lhs_.size(); ++i)
    {
        lhsThermo +=
            lhs_[i].stoichCoeff
           *(*thermoDatabase[species_[lhs_[i].index]]).W()
           *(*thermoDatabase[species_[lhs_[i].index]]);
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
    ReactionThermo::thermoType(*thermoDatabase[species[0]]),
    name_("un-named-reaction-" + Foam::name(getNewReactionID())),
    species_(species),
    Tlow_(TlowDefault),
    Thigh_(ThighDefault),
    lhs_(lhs),
    rhs_(rhs)
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
    ReactionThermo::thermoType(r),
    name_(r.name() + "Copy"),
    species_(species),
    Tlow_(r.Tlow()),
    Thigh_(r.Thigh()),
    lhs_(r.lhs_),
    rhs_(r.rhs_)
{}


template<class ReactionThermo>
Foam::Reaction<ReactionThermo>::specieCoeffs::specieCoeffs
(
    const speciesTable& species,
    Istream& is
)
{
    token t(is);
    if (t.isNumber())
    {
        stoichCoeff = t.number();
        is >> t;
    }
    else
    {
        stoichCoeff = 1;
    }

    exponent = stoichCoeff;

    if (t.isWord())
    {
        word specieName = t.wordToken();

        size_t i = specieName.find('^');

        if (i != word::npos)
        {
            string exponentStr = specieName
            (
                i + 1,
                specieName.size() - i - 1
            );
            exponent = atof(exponentStr.c_str());
            specieName = specieName(0, i);
        }

        if (species.contains(specieName))
        {
            index = species[specieName];
        }
        else
        {
            index = -1;
        }
    }
    else
    {
        FatalIOErrorInFunction(is)
            << "Expected a word but found " << t.info()
            << exit(FatalIOError);
    }
}


template<class ReactionThermo>
void Foam::Reaction<ReactionThermo>::setLRhs
(
    Istream& is,
    const speciesTable& species,
    List<specieCoeffs>& lhs,
    List<specieCoeffs>& rhs
)
{
    DynamicList<specieCoeffs> dlrhs;

    while (is.good())
    {
        dlrhs.append(specieCoeffs(species, is));

        if (dlrhs.last().index != -1)
        {
            token t(is);
            if (t.isPunctuation())
            {
                if (t == token::ADD)
                {
                }
                else if (t == token::ASSIGN)
                {
                    lhs = dlrhs.shrink();
                    dlrhs.clear();
                }
                else
                {
                    rhs = dlrhs.shrink();
                    is.putBack(t);
                    return;
                }
            }
            else
            {
                rhs = dlrhs.shrink();
                is.putBack(t);
                return;
            }
        }
        else
        {
            dlrhs.remove();
            if (is.good())
            {
                token t(is);
                if (t.isPunctuation())
                {
                    if (t == token::ADD)
                    {
                    }
                    else if (t == token::ASSIGN)
                    {
                        lhs = dlrhs.shrink();
                        dlrhs.clear();
                    }
                    else
                    {
                        rhs = dlrhs.shrink();
                        is.putBack(t);
                        return;
                    }
                }
            }
            else
            {
                if (!dlrhs.empty())
                {
                    rhs = dlrhs.shrink();
                }
                return;
            }
        }
    }

    FatalIOErrorInFunction(is)
        << "Cannot continue reading reaction data from stream"
        << exit(FatalIOError);
}


template<class ReactionThermo>
Foam::Reaction<ReactionThermo>::Reaction
(
    const speciesTable& species,
    const HashPtrTable<ReactionThermo>& thermoDatabase,
    const dictionary& dict
)
:
    ReactionThermo::thermoType(*thermoDatabase[species[0]]),
    name_(dict.dictName()),
    species_(species),
    Tlow_(dict.lookupOrDefault<scalar>("Tlow", TlowDefault)),
    Thigh_(dict.lookupOrDefault<scalar>("Thigh", ThighDefault))
{
    setLRhs
    (
        IStringStream(dict.lookup("reaction"))(),
        species_,
        lhs_,
        rhs_
    );
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

    typename dictionaryConstructorTable::iterator cstrIter
        = dictionaryConstructorTablePtr_->find(reactionTypeName);

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


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ReactionThermo>
void Foam::Reaction<ReactionThermo>::write(Ostream& os) const
{
    OStringStream reaction;
    os.writeKeyword("reaction") << reactionStr(reaction)
        << token::END_STATEMENT << nl;
}


template<class ReactionThermo>
void Foam::Reaction<ReactionThermo>::ddot
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    scalarField& d
) const
{
}


template<class ReactionThermo>
void Foam::Reaction<ReactionThermo>::fdot
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    scalarField& f
) const
{
}


template<class ReactionThermo>
void Foam::Reaction<ReactionThermo>::omega
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    scalarField& dcdt
) const
{
    scalar pf, cf, pr, cr;
    label lRef, rRef;

    scalar omegaI = omega
    (
        p, T, c, pf, cf, lRef, pr, cr, rRef
    );

    forAll(lhs_, i)
    {
        const label si = lhs_[i].index;
        const scalar sl = lhs_[i].stoichCoeff;
        dcdt[si] -= sl*omegaI;
    }
    forAll(rhs_, i)
    {
        const label si = rhs_[i].index;
        const scalar sr = rhs_[i].stoichCoeff;
        dcdt[si] += sr*omegaI;
    }
}


template<class ReactionThermo>
Foam::scalar Foam::Reaction<ReactionThermo>::omega
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    scalar& pf,
    scalar& cf,
    label& lRef,
    scalar& pr,
    scalar& cr,
    label& rRef
) const
{

    scalar clippedT = min(max(T, this->Tlow()), this->Thigh());

    const scalar kf = this->kf(p, clippedT, c);
    const scalar kr = this->kr(kf, p, clippedT, c);

    pf = 1;
    pr = 1;

    const label Nl = lhs_.size();
    const label Nr = rhs_.size();

    label slRef = 0;
    lRef = lhs_[slRef].index;

    pf = kf;
    for (label s = 1; s < Nl; s++)
    {
        const label si = lhs_[s].index;

        if (c[si] < c[lRef])
        {
            const scalar exp = lhs_[slRef].exponent;
            pf *= pow(max(c[lRef], 0), exp);
            lRef = si;
            slRef = s;
        }
        else
        {
            const scalar exp = lhs_[s].exponent;
            pf *= pow(max(c[si], 0), exp);
        }
    }
    cf = max(c[lRef], 0);

    {
        const scalar exp = lhs_[slRef].exponent;
        if (exp < 1)
        {
            if (cf > small)
            {
                pf *= pow(cf, exp - 1);
            }
            else
            {
                pf = 0;
            }
        }
        else
        {
            pf *= pow(cf, exp - 1);
        }
    }

    label srRef = 0;
    rRef = rhs_[srRef].index;

    // Find the matrix element and element position for the rhs
    pr = kr;
    for (label s = 1; s < Nr; s++)
    {
        const label si = rhs_[s].index;
        if (c[si] < c[rRef])
        {
            const scalar exp = rhs_[srRef].exponent;
            pr *= pow(max(c[rRef], 0), exp);
            rRef = si;
            srRef = s;
        }
        else
        {
            const scalar exp = rhs_[s].exponent;
            pr *= pow(max(c[si], 0), exp);
        }
    }
    cr = max(c[rRef], 0);

    {
        const scalar exp = rhs_[srRef].exponent;
        if (exp < 1)
        {
            if (cr > small)
            {
                pr *= pow(cr, exp - 1);
            }
            else
            {
                pr = 0;
            }
        }
        else
        {
            pr *= pow(cr, exp - 1);
        }
    }

    return pf*cf - pr*cr;
}


template<class ReactionThermo>
void Foam::Reaction<ReactionThermo>::dwdc
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    scalarSquareMatrix& J,
    scalarField& dcdt,
    scalar& omegaI,
    scalar& kfwd,
    scalar& kbwd,
    const bool reduced,
    const List<label>& c2s
) const
{
    scalar pf, cf, pr, cr;
    label lRef, rRef;

    omegaI = omega(p, T, c, pf, cf, lRef, pr, cr, rRef);

    forAll(lhs_, i)
    {
        const label si = reduced ? c2s[lhs_[i].index] : lhs_[i].index;
        const scalar sl = lhs_[i].stoichCoeff;
        dcdt[si] -= sl*omegaI;
    }
    forAll(rhs_, i)
    {
        const label si = reduced ? c2s[rhs_[i].index] : rhs_[i].index;
        const scalar sr = rhs_[i].stoichCoeff;
        dcdt[si] += sr*omegaI;
    }

    kfwd = this->kf(p, T, c);
    kbwd = this->kr(kfwd, p, T, c);

    forAll(lhs_, j)
    {
        const label sj = reduced ? c2s[lhs_[j].index] : lhs_[j].index;
        scalar kf = kfwd;
        forAll(lhs_, i)
        {
            const label si = lhs_[i].index;
            const scalar el = lhs_[i].exponent;
            if (i == j)
            {
                if (el < 1)
                {
                    if (c[si] > SMALL)
                    {
                        kf *= el*pow(c[si] + VSMALL, el - 1);
                    }
                    else
                    {
                        kf = 0;
                    }
                }
                else
                {
                    kf *= el*pow(c[si], el - 1);
                }
            }
            else
            {
                kf *= pow(c[si], el);
            }
        }

        forAll(lhs_, i)
        {
            const label si = reduced ? c2s[lhs_[i].index] : lhs_[i].index;
            const scalar sl = lhs_[i].stoichCoeff;
            J(si, sj) -= sl*kf;
        }
        forAll(rhs_, i)
        {
            const label si = reduced ? c2s[rhs_[i].index] : rhs_[i].index;
            const scalar sr = rhs_[i].stoichCoeff;
            J(si, sj) += sr*kf;
        }
    }

    forAll(rhs_, j)
    {
        const label sj = reduced ? c2s[rhs_[j].index] : rhs_[j].index;
        scalar kr = kbwd;
        forAll(rhs_, i)
        {
            const label si = rhs_[i].index;
            const scalar er = rhs_[i].exponent;
            if (i == j)
            {
                if (er < 1)
                {
                    if (c[si] > SMALL)
                    {
                        kr *= er*pow(c[si] + VSMALL, er - 1);
                    }
                    else
                    {
                        kr = 0;
                    }
                }
                else
                {
                    kr *= er*pow(c[si], er - 1);
                }
            }
            else
            {
                kr *= pow(c[si], er);
            }
        }

        forAll(lhs_, i)
        {
            const label si = reduced ? c2s[lhs_[i].index] : lhs_[i].index;
            const scalar sl = lhs_[i].stoichCoeff;
            J(si, sj) += sl*kr;
        }
        forAll(rhs_, i)
        {
            const label si = reduced ? c2s[rhs_[i].index] : rhs_[i].index;
            const scalar sr = rhs_[i].stoichCoeff;
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
        this->dcidc(p, T, c, dcidc);

        forAll(beta, j)
        {
            label sj = beta[j].first();
            sj = reduced ? c2s[sj] : sj;
            if (sj != -1)
            {
                forAll(lhs_, i)
                {
                    const label si =
                        reduced ? c2s[lhs_[i].index] : lhs_[i].index;
                    const scalar sl = lhs_[i].stoichCoeff;
                    J(si, sj) -= sl*dcidc[j]*omegaI;
                }
                forAll(rhs_, i)
                {
                    const label si =
                        reduced ? c2s[rhs_[i].index] : rhs_[i].index;
                    const scalar sr = rhs_[i].stoichCoeff;
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
    const scalar omegaI,
    const scalar kfwd,
    const scalar kbwd,
    scalarSquareMatrix& J,
    const bool reduced,
    const List<label>& c2s,
    const label indexT
) const
{
    scalar kf = kfwd;
    scalar kr = kbwd;

    scalar dkfdT = this->dkfdT(p, T, c);
    scalar dkrdT = this->dkrdT(p, T, c, dkfdT, kr);

    scalar sumExp = 0.0;
    forAll(lhs_, i)
    {
        const label si = lhs_[i].index;
        const scalar el = lhs_[i].exponent;
        const scalar cExp = pow(c[si], el);
        dkfdT *= cExp;
        kf *= cExp;
        sumExp += el;
    }
    kf *= -sumExp/T;

    sumExp = 0.0;
    forAll(rhs_, i)
    {
        const label si = rhs_[i].index;
        const scalar er = rhs_[i].exponent;
        const scalar cExp = pow(c[si], er);
        dkrdT *= cExp;
        kr *= cExp;
        sumExp += er;
    }
    kr *= -sumExp/T;

    // dqidT includes the third-body (or pressure dependent) effect
    scalar dqidT = dkfdT - dkrdT + kf - kr;

    // For reactions including third-body efficiencies or pressure dependent
    // reaction, an additional term is needed
    scalar dcidT = this->dcidT(p, T, c);
    dcidT *= omegaI;

    // J(i, indexT) = sum_reactions nu_i dqdT
    forAll(lhs_, i)
    {
        const label si = reduced ? c2s[lhs_[i].index] : lhs_[i].index;
        const scalar sl = lhs_[i].stoichCoeff;
        J(si, indexT) -= sl*(dqidT + dcidT);
    }
    forAll(rhs_, i)
    {
        const label si = reduced ? c2s[rhs_[i].index] : rhs_[i].index;
        const scalar sr = rhs_[i].stoichCoeff;
        J(si, indexT) += sr*(dqidT + dcidT);
    }
}


template<class ReactionThermo>
const Foam::speciesTable& Foam::Reaction<ReactionThermo>::species() const
{
    return species_;
}


// ************************************************************************* //
