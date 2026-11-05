/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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

#include "zonal_surfaceInterpolationScheme.H"
#include "UCompactListList.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::zonal<Type>::zonal
(
    const fvMesh& mesh,
    const surfaceScalarField& faceFlux,
    Istream& schemeData
)
:
    surfaceInterpolationScheme<Type>(mesh),
    faceFlux_(faceFlux),
    tScheme_
    (
        surfaceInterpolationScheme<Type>::New
        (
            mesh,
            faceFlux,
            schemeData
        )
    ),
    zoneSchemesDict_(schemeData),
    zoneSchemes_(zoneSchemesDict_.size())
{
    label i = 0;
    forAllConstIter(dictionary, zoneSchemesDict_, iter)
    {
        zoneSchemes_.set
        (
            i++,
            surfaceInterpolationScheme<Type>::New
            (
                mesh,
                faceFlux,
                iter().stream()
            )
        );
    }
}


template<class Type>
Foam::zonal<Type>::zonal
(
    const fvMesh& mesh,
    Istream& schemeData
)
:
    zonal
    (
        mesh,
        mesh.lookupObject<surfaceScalarField>(word(schemeData)),
        schemeData
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::surfaceScalarField>
Foam::zonal<Type>::weights(const VolField<Type>& vf) const
{
    tmp<surfaceScalarField> tw = tScheme_().weights(vf);
    surfaceScalarField& w = tw.ref();

    const fvMesh& mesh = this->mesh();
    const faceZoneList& faceZones = mesh.faceZones();

    label zsi = 0;
    forAllConstIter(dictionary, zoneSchemesDict_, iter)
    {
        tmp<surfaceScalarField> tzw = zoneSchemes_[zsi].weights(vf);
        const surfaceScalarField& zw = tzw.ref();

        forAll(faceZones, fzi)
        {
            const faceZone& fZone = faceZones[fzi];

            if (iter().keyword().match(fZone.name()))
            {
                forAll(fZone, fi)
                {
                    const label facei = fZone[fi];

                    if (mesh.isInternalFace(facei))
                    {
                        w[facei] = zw[facei];
                    }
                    else
                    {
                        const label bFacei = facei - mesh.nInternalFaces();
                        const labelUList patches
                        (
                            mesh.polyBFacePatches()[bFacei]
                        );
                        const labelUList patchFaces
                        (
                            mesh.polyBFacePatchFaces()[bFacei]
                        );

                        forAll(patches, pi)
                        {
                            w.boundaryFieldRef()[patches[pi]][patchFaces[pi]] =
                                zw.boundaryField()[patches[pi]][patchFaces[pi]];
                        }
                    }
                }
            }
        }

        zsi++;
    }

    return tw;
}


template<class Type>
bool Foam::zonal<Type>::corrected() const
{
    bool corrected = tScheme_().corrected();
    forAll(zoneSchemes_, i)
    {
        if (corrected) break;
        corrected = corrected || zoneSchemes_[i].corrected();
    }
    return corrected;
}


template<class Type>
Foam::tmp<Foam::SurfaceField<Type>>
Foam::zonal<Type>::correction(const VolField<Type>& vf) const
{
    if (corrected())
    {
        const fvMesh& mesh = this->mesh();

        tmp<SurfaceField<Type>> tcorr;

        if (tScheme_().corrected())
        {
            tcorr = tScheme_().correction(vf);
        }
        else
        {
            tcorr = SurfaceField<Type>::New
            (
                word("zonal::correction(", vf.name(), ')'),
                mesh,
                dimensioned<Type>(vf.name(), vf.dimensions(), Zero)
            );
        }

        SurfaceField<Type>& corr = tcorr.ref();

        const faceZoneList& faceZones = mesh.faceZones();

        label zsi = 0;
        forAllConstIter(dictionary, zoneSchemesDict_, iter)
        {
            if (zoneSchemes_[zsi].corrected())
            {
                tmp<SurfaceField<Type>> tzcorr
                (
                    zoneSchemes_[zsi].correction(vf)
                );
                const SurfaceField<Type>& zcorr = tzcorr.ref();

                forAll(faceZones, fzi)
                {
                    const faceZone& fZone = faceZones[fzi];

                    if (iter().keyword().match(fZone.name()))
                    {
                        forAll(fZone, fi)
                        {
                            const label facei = fZone[fi];

                            if (mesh.isInternalFace(facei))
                            {
                                corr[facei] = zcorr[facei];
                            }
                            else
                            {
                                const label bFacei =
                                    facei - mesh.nInternalFaces();
                                const labelUList patches
                                (
                                    mesh.polyBFacePatches()[bFacei]
                                );
                                const labelUList patchFaces
                                (
                                    mesh.polyBFacePatchFaces()[bFacei]
                                );

                                forAll(patches, pi)
                                {
                                    corr.boundaryFieldRef()
                                    [patches[pi]][patchFaces[pi]] =
                                        zcorr.boundaryField()
                                        [patches[pi]][patchFaces[pi]];
                                }
                            }
                        }
                    }
                }
            }
            else
            {
                forAll(faceZones, fzi)
                {
                    const faceZone& fZone = faceZones[fzi];

                    if (iter().keyword().match(fZone.name()))
                    {
                        forAll(fZone, fi)
                        {
                            const label facei = fZone[fi];

                            if (mesh.isInternalFace(facei))
                            {
                                corr[facei] = Zero;
                            }
                            else
                            {
                                const label bFacei =
                                    facei - mesh.nInternalFaces();
                                const labelUList patches
                                (
                                    mesh.polyBFacePatches()[bFacei]
                                );
                                const labelUList patchFaces
                                (
                                    mesh.polyBFacePatchFaces()[bFacei]
                                );

                                forAll(patches, pi)
                                {
                                    corr.boundaryFieldRef()
                                    [patches[pi]][patchFaces[pi]] = Zero;
                                }
                            }
                        }
                    }
                }
            }

            zsi++;
        }

        return tcorr;
    }
    else
    {
        return tmp<SurfaceField<Type>>(nullptr);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makeSurfaceInterpolationScheme(zonal);
}

// ************************************************************************* //
