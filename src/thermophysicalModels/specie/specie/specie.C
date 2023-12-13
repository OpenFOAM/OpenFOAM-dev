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

#include "specie.H"

/* * * * * * * * * * * * * * * public constants  * * * * * * * * * * * * * * */

namespace Foam
{
    defineTypeNameAndDebug(specie, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::specie::specie(const word& name, const dictionary& dict)
:
    name_(name),
    Y_(dict.subDict("specie").lookupOrDefault("massFraction", 1.0)),
    molWeight_(dict.subDict("specie").lookup<scalar>("molWeight"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::specie::write(Ostream& os) const
{
    dictionary dict("specie");
    if (Y_ != 1)
    {
        dict.add("massFraction", Y_);
    }
    dict.add("molWeight", molWeight_);
    os  << indent << dict.dictName() << dict;
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const specie& st)
{
    st.write(os);
    os.check("Ostream& operator<<(Ostream& os, const specie& st)");
    return os;
}


// ************************************************************************* //
