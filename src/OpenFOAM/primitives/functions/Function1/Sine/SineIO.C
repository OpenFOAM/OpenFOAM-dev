/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenFOAM Foundation
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

#include "Sine.H"

// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class Type>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const Function1Types::Sine<Type>& se
)
{
    os  << static_cast<const Function1<Type>& >(se)
        << token::SPACE << se.t0_
        << token::SPACE << se.amplitude_
        << token::SPACE << se.frequency_
        << token::SPACE << se.scale_
        << token::SPACE << se.level_;

    // Check state of Ostream
    os.check
    (
        "Ostream& operator<<(Ostream&, const Sine<Type>&)"
    );

    return os;
}


template<class Type>
void Foam::Function1Types::Sine<Type>::writeData(Ostream& os) const
{
    Function1<Type>::writeData(os);
    os  << token::END_STATEMENT << nl;
    os  << indent << word(this->name() + "Coeffs") << nl;
    os  << indent << token::BEGIN_BLOCK << incrIndent << nl;
    os.writeKeyword("t0") << t0_ << token::END_STATEMENT << nl;
    os.writeKeyword("amplitude") << amplitude_ << token::END_STATEMENT << nl;
    os.writeKeyword("frequency") << frequency_ << token::END_STATEMENT << nl;
    os.writeKeyword("scale") << scale_ << token::END_STATEMENT << nl;
    os.writeKeyword("level") << level_ << token::END_STATEMENT << nl;
    os  << decrIndent << indent << token::END_BLOCK << endl;
}


// ************************************************************************* //
