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

#include "correlationFunction.H"
#include "IOstreams.H"

template<class Type>
bool Foam::correlationFunction<Type>::writeAveraged(Ostream& os) const
{
    Field<scalar> averageCF(averaged());

    forAll(averageCF, v)
    {
        os  << v*sampleInterval()
            << token::SPACE
            << averageCF[v]
            << nl;
    }

    return os.good();
}


template<class Type>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const correlationFunction<Type>& cF
)
{
    os  << cF.duration()
        << nl << cF.sampleInterval()
        << nl << cF.averagingInterval()
        << nl << cF.sampleSteps()
        << nl << cF.tZeroBuffers()
        << nl << static_cast<const bufferedAccumulator<scalar>&>(cF);

    // Check state of Ostream
    os.check
    (
        "Foam::Ostream& Foam::operator<<"
        "(Ostream&, const correlationFunction<Type>&)"
    );

    return os;
}


// ************************************************************************* //
