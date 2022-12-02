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

#include "pointLinear.H"
#include "fvMesh.H"
#include "volPointInterpolation.H"
#include "triangle.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::SurfaceField<Type>>
Foam::pointLinear<Type>::
correction
(
    const VolField<Type>& vf
) const
{
    const fvMesh& mesh = this->mesh();

    PointField<Type> pvf
    (
        volPointInterpolation::New(mesh).interpolate(vf)
    );

    tmp<SurfaceField<Type>> tsfCorr =
        linearInterpolate(vf);

    Field<Type>& sfCorr = tsfCorr.ref().primitiveFieldRef();

    const pointField& points = mesh.points();
    const pointField& C = mesh.C();
    const faceList& faces = mesh.faces();
    const scalarField& w = mesh.weights();
    const labelList& owner = mesh.owner();
    const labelList& neighbour = mesh.neighbour();

    forAll(sfCorr, facei)
    {
        point pi =
            w[owner[facei]]*C[owner[facei]]
          + (1.0 - w[owner[facei]])*C[neighbour[facei]];

        const face& f = faces[facei];

        scalar at = triangle<point, const point&>
        (
            pi,
            points[f[0]],
            points[f[f.size()-1]]
        ).mag();

        scalar sumAt = at;
        Type sumPsip = at*(1.0/3.0)*
        (
            sfCorr[facei]
          + pvf[f[0]]
          + pvf[f[f.size()-1]]
        );

        for (label pointi=1; pointi<f.size(); pointi++)
        {
            at = triangle<point, const point&>
            (
                pi,
                points[f[pointi]],
                points[f[pointi-1]]
            ).mag();

            sumAt += at;
            sumPsip += at*(1.0/3.0)*
            (
                sfCorr[facei]
              + pvf[f[pointi]]
              + pvf[f[pointi-1]]
            );

        }

        sfCorr[facei] = sumPsip/sumAt - sfCorr[facei];
    }


    typename SurfaceField<Type>::
        Boundary& bSfCorr = tsfCorr.ref().boundaryFieldRef();

    forAll(bSfCorr, patchi)
    {
        fvsPatchField<Type>& pSfCorr = bSfCorr[patchi];

        if (pSfCorr.coupled())
        {
            const fvPatch& fvp = mesh.boundary()[patchi];
            const scalarField& pWghts = mesh.weights().boundaryField()[patchi];
            const labelUList& pOwner = fvp.faceCells();
            const vectorField& pNbrC = mesh.C().boundaryField()[patchi];

            forAll(pOwner, facei)
            {
                label own = pOwner[facei];

                point pi =
                    pWghts[facei]*C[own]
                  + (1.0 - pWghts[facei])*pNbrC[facei];

                const face& f = faces[facei+fvp.start()];

                scalar at = triangle<point, const point&>
                (
                    pi,
                    points[f[0]],
                    points[f[f.size()-1]]
                ).mag();

                scalar sumAt = at;
                Type sumPsip = at*(1.0/3.0)*
                (
                    pSfCorr[facei]
                  + pvf[f[0]]
                  + pvf[f[f.size()-1]]
                );

                for (label pointi=1; pointi<f.size(); pointi++)
                {
                    at = triangle<point, const point&>
                    (
                        pi,
                        points[f[pointi]],
                        points[f[pointi-1]]
                    ).mag();

                    sumAt += at;
                    sumPsip += at*(1.0/3.0)*
                    (
                        pSfCorr[facei]
                      + pvf[f[pointi]]
                      + pvf[f[pointi-1]]
                    );

                }

                pSfCorr[facei] = sumPsip/sumAt - pSfCorr[facei];
            }
        }
        else
        {
            pSfCorr = Zero;
        }
    }

    return tsfCorr;
}


namespace Foam
{
    makeSurfaceInterpolationScheme(pointLinear);
}

// ************************************************************************* //
