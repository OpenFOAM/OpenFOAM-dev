/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019-2021 OpenFOAM Foundation
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
#include "None.H"
#include "addToRunTimeSelectionTable.H"

#include "thermodynamicConstants.H"
using namespace Foam::constant::thermodynamic;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(liquid, 0);
    addToRunTimeSelectionTable(liquidProperties, liquid, dictionary);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::autoPtr<Foam::Function1<Foam::scalar>> Foam::liquid::New
(
    const word& name,
    const dictionary& dict
)
{
    if (dict.isDict(name))
    {
        return Function1<scalar>::New(name, dict);
    }
    else
    {
        return autoPtr<Function1<scalar>>
        (
            new Function1s::None<scalar>(name, dict)
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::liquid::liquid(const dictionary& dict)
:
    liquidProperties(dict),
    rho_(New("rho", dict)),
    pv_(New("pv", dict)),
    hl_(New("hl", dict)),
    Cp_(New("Cp", dict)),
    h_(New("h", dict)),
    Cpg_(New("Cpg", dict)),
    B_(New("B", dict)),
    mu_(New("mu", dict)),
    mug_(New("mug", dict)),
    kappa_(New("kappa", dict)),
    kappag_(New("kappag", dict)),
    sigma_(New("sigma", dict)),
    D_(New("D", dict)),
    Hf_(h_->value(Tstd))
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
