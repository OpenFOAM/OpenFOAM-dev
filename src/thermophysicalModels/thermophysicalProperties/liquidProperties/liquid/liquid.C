/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019-2020 OpenFOAM Foundation
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
    rho_(thermophysicalFunction::New("rho", dict)),
    pv_(thermophysicalFunction::New("pv", dict)),
    hl_(thermophysicalFunction::New("hl", dict)),
    Cp_(thermophysicalFunction::New("Cp", dict)),
    h_(thermophysicalFunction::New("h", dict)),
    Cpg_(thermophysicalFunction::New("Cpg", dict)),
    B_(thermophysicalFunction::New("B", dict)),
    mu_(thermophysicalFunction::New("mu", dict)),
    mug_(thermophysicalFunction::New("mug", dict)),
    kappa_(thermophysicalFunction::New("kappa", dict)),
    kappag_(thermophysicalFunction::New("kappag", dict)),
    sigma_(thermophysicalFunction::New("sigma", dict)),
    D_(thermophysicalFunction::New("D", dict))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::liquid::write(Ostream& os) const
{
    liquidProperties::write(os); os << nl;
    rho_->write(os); os << nl;
    pv_->write(os); os << nl;
    hl_->write(os); os << nl;
    Cp_->write(os); os << nl;
    h_->write(os); os << nl;
    Cpg_->write(os); os << nl;
    B_->write(os); os << nl;
    mu_->write(os); os << nl;
    mug_->write(os); os << nl;
    kappa_->write(os); os << nl;
    kappag_->write(os); os << nl;
    sigma_->write(os); os << nl;
    D_->write(os); os << endl;
}


// ************************************************************************* //
