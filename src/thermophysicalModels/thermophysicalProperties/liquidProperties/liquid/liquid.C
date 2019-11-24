/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019 OpenFOAM Foundation
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

#include "liquid.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(liquid, 0);
    addToRunTimeSelectionTable(liquidProperties, liquid, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::liquid::liquid(const dictionary& dict)
:
    liquidProperties(dict),
    rho_(thermophysicalFunction::New(dict, "rho")),
    pv_(thermophysicalFunction::New(dict, "pv")),
    hl_(thermophysicalFunction::New(dict, "hl")),
    Cp_(thermophysicalFunction::New(dict, "Cp")),
    h_(thermophysicalFunction::New(dict, "h")),
    Cpg_(thermophysicalFunction::New(dict, "Cpg")),
    B_(thermophysicalFunction::New(dict, "B")),
    mu_(thermophysicalFunction::New(dict, "mu")),
    mug_(thermophysicalFunction::New(dict, "mug")),
    kappa_(thermophysicalFunction::New(dict, "kappa")),
    kappag_(thermophysicalFunction::New(dict, "kappag")),
    sigma_(thermophysicalFunction::New(dict, "sigma")),
    D_(thermophysicalFunction::New(dict, "D"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::liquid::writeData(Ostream& os) const
{
    liquidProperties::writeData(os); os << nl;
    rho_->writeData(os); os << nl;
    pv_->writeData(os); os << nl;
    hl_->writeData(os); os << nl;
    Cp_->writeData(os); os << nl;
    h_->writeData(os); os << nl;
    Cpg_->writeData(os); os << nl;
    B_->writeData(os); os << nl;
    mu_->writeData(os); os << nl;
    mug_->writeData(os); os << nl;
    kappa_->writeData(os); os << nl;
    kappag_->writeData(os); os << nl;
    sigma_->writeData(os); os << nl;
    D_->writeData(os); os << endl;
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const liquid& l)
{
    l.writeData(os);
    return os;
}


// ************************************************************************* //
