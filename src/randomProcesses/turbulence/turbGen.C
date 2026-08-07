/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2026 OpenFOAM Foundation
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

#include "fft.H"
#include "turbGen.H"
#include "Kmesh.H"
#include "primitiveFields.H"
#include "Ek.H"
#include "mathematicalConstants.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::turbGen::turbGen(const Kmesh& k, const scalar EA, const scalar K0)
:
    K_(k),
    Ea_(EA),
    k0_(K0),
    rndGen_(label(0))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::vectorField Foam::turbGen::U()
{
    const vectorField& K(K_);

    vectorField s(K.size());
    scalarField rndPhases(K.size());

    forAll(K, i)
    {
        s[i] = rndGen_.sample01<vector>();
        rndPhases[i] = rndGen_.scalar01();
    }

    s = K ^ s;
    s = s/(mag(s) + 1.0e-20);

    s = Ek(Ea_, k0_, mag(K))*s;

    complexVectorField up
    (
        fft::reverseTransform
        (
            ComplexField
            (
                eval(cos(constant::mathematical::twoPi*rndPhases)*s),
                eval(sin(constant::mathematical::twoPi*rndPhases)*s)
            ),
            K_.nn()
        )
    );

    return ReImSum(up);
}


// ************************************************************************* //
