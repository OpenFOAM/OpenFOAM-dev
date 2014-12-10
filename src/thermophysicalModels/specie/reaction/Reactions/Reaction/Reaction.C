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

#include "Reaction.H"
#include "DynamicList.H"

// * * * * * * * * * * * * * * * * Static Data * * * * * * * * * * * * * * * //

template<class ReactionThermo>
Foam::label Foam::Reaction<ReactionThermo>::nUnNamedReactions = 0;

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
        if (mag(lhs_[i].stoichCoeff - 1) > SMALL)
        {
            reaction << lhs_[i].stoichCoeff;
        }
        reaction << species_[lhs_[i].index];
        if (mag(lhs_[i].exponent - lhs_[i].stoichCoeff) > SMALL)
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
        if (mag(rhs_[i].stoichCoeff - 1) > SMALL)
        {
            reaction << rhs_[i].stoichCoeff;
        }
        reaction << species_[rhs_[i].index];
        if (mag(rhs_[i].exponent - rhs_[i].stoichCoeff) > SMALL)
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
    if (rhs_.size() > 0)
    {
        ReactionThermo::thermoType::operator=
        (
            rhs_[0].stoichCoeff*(*thermoDatabase[species_[rhs_[0].index]])
        );

        for (label i=1; i<rhs_.size(); ++i)
        {
            this->operator+=
            (
                rhs_[i].stoichCoeff*(*thermoDatabase[species_[rhs_[i].index]])
            );
        }
    }

    forAll(lhs_, i)
    {
        this->operator-=
        (
            lhs_[i].stoichCoeff*(*thermoDatabase[species_[lhs_[i].index]])
        );
    }
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
        stoichCoeff = 1.0;
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
        FatalIOErrorIn("Reaction<ReactionThermo>::lrhs(Istream& is)", is)
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

    FatalIOErrorIn("Reaction<ReactionThermo>::setLRhs(Istream& is)", is)
        << "Cannot continue reading reaction data from stream"
        << exit(FatalIOError);
}


template<class ReactionThermo>
Foam::Reaction<ReactionThermo>::Reaction
(
    const speciesTable& species,
    const HashPtrTable<ReactionThermo>& thermoDatabase,
    Istream& is
)
:
    ReactionThermo::thermoType(*thermoDatabase[species[0]]),
    name_("un-named-reaction" + Foam::name(getNewReactionID())),
    species_(species)
{
    setLRhs(is, species, lhs_, rhs_);
    setThermo(thermoDatabase);
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
    species_(species)
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
Foam::autoPtr<Foam::Reaction<ReactionThermo> >
Foam::Reaction<ReactionThermo>::New
(
    const speciesTable& species,
    const HashPtrTable<ReactionThermo>& thermoDatabase,
    Istream& is
)
{
    if (is.eof())
    {
        FatalIOErrorIn
        (
            "Reaction<ReactionThermo>::New(const speciesTable&, "
            " const HashPtrTable<ReactionThermo>&, Istream&)",
            is
        )   << "Reaction type not specified" << nl << nl
            << "Valid Reaction types are :" << nl
            << IstreamConstructorTablePtr_->sortedToc()
            << exit(FatalIOError);
    }

    const word reactionTypeName(is);

    typename IstreamConstructorTable::iterator cstrIter
        = IstreamConstructorTablePtr_->find(reactionTypeName);

    if (cstrIter == IstreamConstructorTablePtr_->end())
    {
        FatalIOErrorIn
        (
            "Reaction<ReactionThermo>::New(const speciesTable&, "
            " const HashPtrTable<ReactionThermo>&, Istream&)",
            is
        )   << "Unknown reaction type "
            << reactionTypeName << nl << nl
            << "Valid reaction types are :" << nl
            << IstreamConstructorTablePtr_->sortedToc()
            << exit(FatalIOError);
    }

    return autoPtr<Reaction<ReactionThermo> >
    (
        cstrIter()(species, thermoDatabase, is)
    );
}


template<class ReactionThermo>
Foam::autoPtr<Foam::Reaction<ReactionThermo> >
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
        FatalErrorIn
        (
            "Reaction<ReactionThermo>::New"
            "("
                "const speciesTable&, "
                "const HashPtrTable<ReactionThermo>&, "
                "const dictionary&"
            ")"
        )   << "Unknown reaction type "
            << reactionTypeName << nl << nl
            << "Valid reaction types are :" << nl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<Reaction<ReactionThermo> >
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
Foam::scalar Foam::Reaction<ReactionThermo>::kf
(
    const scalar p,
    const scalar T,
    const scalarField& c
) const
{
    return 0.0;
}


template<class ReactionThermo>
Foam::scalar Foam::Reaction<ReactionThermo>::kr
(
    const scalar kfwd,
    const scalar p,
    const scalar T,
    const scalarField& c
) const
{
    return 0.0;
}


template<class ReactionThermo>
Foam::scalar Foam::Reaction<ReactionThermo>::kr
(
    const scalar p,
    const scalar T,
    const scalarField& c
) const
{
    return 0.0;
}


template<class ReactionThermo>
const Foam::speciesTable& Foam::Reaction<ReactionThermo>::species() const
{
    return species_;
}


template<class ReactionThermo>
const Foam::speciesTable& Foam::Reaction<ReactionThermo>::gasSpecies() const
{
    notImplemented
    (
        "const speciesTable& gasSpecies() const"
        " for this reaction"
    );
    return *reinterpret_cast<speciesTable*>(0);
}


template<class ReactionThermo>
const Foam::List<typename Foam::Reaction<ReactionThermo>::specieCoeffs>&
Foam::Reaction<ReactionThermo>::glhs() const
{
    notImplemented
    (
        "inline const List<typename Reaction<ReactionThermo>::specieCoeffs>&"
        "Reaction<ReactionThermo>::glhs()"
    );
    return *reinterpret_cast<List<specieCoeffs>*>(0);
}


template<class ReactionThermo>
const Foam::List<typename Foam::Reaction<ReactionThermo>::specieCoeffs>&
Foam::Reaction<ReactionThermo>::grhs() const
{
    notImplemented
    (
        "inline const List<typename Reaction<ReactionThermo>::specieCoeffs>&"
        "Reaction<ReactionThermo>::grhs()"
    );
    return *reinterpret_cast<List<specieCoeffs>*>(0);
}

// ************************************************************************* //
