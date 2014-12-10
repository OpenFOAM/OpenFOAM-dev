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

#include "geomDecomp.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::geomDecomp::geomDecomp
(
    const dictionary& decompositionDict,
    const word& derivedType
)
:
    decompositionMethod(decompositionDict),
    geomDecomDict_(decompositionDict.subDict(derivedType + "Coeffs")),
    n_(geomDecomDict_.lookup("n")),
    delta_(readScalar(geomDecomDict_.lookup("delta"))),
    rotDelta_(I)
{
    // check that the case makes sense :

    if (nProcessors_ != n_.x()*n_.y()*n_.z())
    {
        FatalErrorIn
        (
            "geomDecomp::geomDecomp"
            "(const dictionary& decompositionDict)"
        )   << "Wrong number of processor divisions in geomDecomp:" << nl
            << "Number of domains    : " << nProcessors_ << nl
            << "Wanted decomposition : " << n_
            << exit(FatalError);
    }

    scalar d = 1 - 0.5*delta_*delta_;
    scalar d2 = sqr(d);

    scalar a = delta_;
    scalar a2 = sqr(a);

    rotDelta_ = tensor
    (
        d2,         -a*d,         a,
        a*d - a2*d,  a*a2 + d2,  -2*a*d,
        a*d2 + a2,   a*d - a2*d,  d2 - a2
    );
}


// ************************************************************************* //
