/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "phaseProperties.H"
#include "dictionaryEntry.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseProperties::phaseProperties(Istream& is)
:
    phase_(UNKNOWN),
    stateLabel_("(unknown)"),
    names_(0),
    Y_(0),
    globalIds_(0),
    globalCarrierIds_(0)
{
    is.check("Foam::phaseProperties::phaseProperties(Istream& is)");

    dictionaryEntry phaseInfo(dictionary::null, is);

    phase_ = phaseTypeNames_[phaseInfo.keyword()];
    stateLabel_ = phaseToStateLabel(phase_);

    if (phaseInfo.size() > 0)
    {
        label nComponents = phaseInfo.size();
        names_.setSize(nComponents, "unknownSpecie");
        Y_.setSize(nComponents, 0.0);
        globalIds_.setSize(nComponents, -1);
        globalCarrierIds_.setSize(nComponents, -1);

        label cmptI = 0;
        forAllConstIter(IDLList<entry>, phaseInfo, iter)
        {
            names_[cmptI] = iter().keyword();
            Y_[cmptI] = readScalar(phaseInfo.lookup(names_[cmptI]));
            cmptI++;
        }

        checkTotalMassFraction();
    }
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, phaseProperties& pp)
{
    is.check
    (
        "Foam::Istream& Foam::operator>>(Istream&, phaseProperties&)"
    );

    dictionaryEntry phaseInfo(dictionary::null, is);

    pp.phase_ = pp.phaseTypeNames_[phaseInfo.keyword()];
    pp.stateLabel_ = pp.phaseToStateLabel(pp.phase_);

    if (phaseInfo.size() > 0)
    {
        label nComponents = phaseInfo.size();

        pp.names_.setSize(nComponents, "unknownSpecie");
        pp.Y_.setSize(nComponents, 0.0);
        pp.globalIds_.setSize(nComponents, -1);
        pp.globalCarrierIds_.setSize(nComponents, -1);

        label cmptI = 0;
        forAllConstIter(IDLList<entry>, phaseInfo, iter)
        {
            pp.names_[cmptI] = iter().keyword();
            pp.Y_[cmptI] = readScalar(phaseInfo.lookup(pp.names_[cmptI]));
            cmptI++;
        }

        pp.checkTotalMassFraction();
    }

    return is;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const phaseProperties& pp)
{
    os.check
    (
        "Foam::Ostream& Foam::operator<<(Ostream&, const phaseProperties&)"
    );

    os  << pp.phaseTypeNames_[pp.phase_] << nl << token::BEGIN_BLOCK << nl
        << incrIndent;

    forAll(pp.names_, cmptI)
    {
        os.writeKeyword(pp.names_[cmptI]) << pp.Y_[cmptI]
            << token::END_STATEMENT << nl;
    }

    os  << decrIndent << token::END_BLOCK << nl;

    os.check
    (
        "Foam::Ostream& Foam::operator<<(Ostream&, const phaseProperties&)"
    );

    return os;
}


// ************************************************************************* //
