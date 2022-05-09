/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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

#include "coupledFvPatch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(coupledFvPatch, 0);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::coupledFvPatch::~coupledFvPatch()
{}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::coupledFvPatch::makeWeights
(
    scalarField& w,
    const vectorField& nbrSf,
    const vectorField& nbrDelta
) const
{
    const vectorField delta(coupledFvPatch::delta());

    const scalarField nfDelta(nf() & delta);

    const scalarField nbrNfDelta((nbrSf/(mag(nbrSf) + vSmall)) & nbrDelta);

    forAll(delta, facei)
    {
        const scalar ndoi = nfDelta[facei];
        const scalar ndni = nbrNfDelta[facei];
        const scalar ndi = ndoi + ndni;

        if (ndni/vGreat < ndi)
        {
            w[facei] = ndni/ndi;
        }
        else
        {
            const scalar doi = mag(delta[facei]);
            const scalar dni = mag(nbrDelta[facei]);
            const scalar di = doi + dni;

            w[facei] = dni/di;
        }
    }
}


Foam::tmp<Foam::vectorField> Foam::coupledFvPatch::delta
(
    const vectorField& nbrDelta
) const
{
    if (transform().transforms())
    {
        return
            coupledFvPatch::delta()
          - transform().transform(nbrDelta);
    }
    else
    {
        return coupledFvPatch::delta() - nbrDelta;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::vectorField> Foam::coupledFvPatch::delta() const
{
    return Cf() - Cn();
}


// ************************************************************************* //
