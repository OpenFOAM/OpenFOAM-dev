/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "kShellIntegration.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::graph Foam::kShellIntegration
(
    const complexVectorField& Ek,
    const Kmesh& K
)
{
    // evaluate the radial component of the spectra as an average
    // over the shells of thickness dk

    graph kShellMeanEk = kShellMean(Ek, K);
    const scalarField& x = kShellMeanEk.x();
    scalarField& y = *kShellMeanEk.begin()();

    // now multiply by 4pi k^2 (the volume of each shell) to get the
    // spectra E(k). int E(k) dk is now the total energy in a box
    // of side 2pi

    y *= sqr(x)*4.0*constant::mathematical::pi;

    // now scale this to get the energy in a box of side l0

    scalar l0(K.sizeOfBox()[0]*(scalar(K.nn()[0])/(scalar(K.nn()[0])-1.0)));
    scalar factor = pow((l0/(2.0*constant::mathematical::pi)),3.0);

    y *= factor;

    // and divide by the number of points in the box, to give the
    // energy density.

    y /= scalar(K.size());

    return kShellMeanEk;
}


// kShellMean : average over the points in a k-shell to evaluate the
// radial part of the energy spectrum.

Foam::graph Foam::kShellMean
(
    const complexVectorField& Ek,
    const Kmesh& K
)
{
    const label tnp = Ek.size();
    const label NoSubintervals = label
    (
        pow(scalar(tnp), 1.0/vector::dim)*pow(1.0/vector::dim, 0.5) - 0.5
    );

    scalarField k1D(NoSubintervals);
    scalarField Ek1D(NoSubintervals);
    scalarField EWeight(NoSubintervals);

    scalar kmax = K.max()*pow(1.0/vector::dim,0.5);
    scalar delta_k = kmax/(NoSubintervals);

    forAll(Ek1D, a)
    {
        k1D[a] = (a + 1)*delta_k;
        Ek1D[a] = 0.0;
        EWeight[a] = 0;
    }

    forAll(K, l)
    {
        scalar kmag = mag(K[l]);

        for (label a=0; a<NoSubintervals; a++)
        {
            if
            (
                kmag <= ((a + 1)*delta_k + delta_k/2.0)
             && kmag > ((a + 1)*delta_k - delta_k/2.0)
            )
            {
                scalar dist = delta_k/2.0 - mag((a + 1)*delta_k - kmag);

                Ek1D[a] += dist*
                magSqr
                (
                    vector
                    (
                        mag(Ek[l].x()),
                        mag(Ek[l].y()),
                        mag(Ek[l].z())
                    )
                 );

                EWeight[a] += dist;
            }
        }
    }

    for (label a=0; a<NoSubintervals; a++)
    {
        if (EWeight[a] > 0)
        {
            Ek1D[a] /= EWeight[a];
        }
    }

    return graph("E(k)", "k", "E(k)", k1D, Ek1D);
}


// ************************************************************************* //
