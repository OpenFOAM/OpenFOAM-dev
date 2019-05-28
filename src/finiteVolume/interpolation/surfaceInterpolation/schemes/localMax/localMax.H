/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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

Class
    Foam::localMax

Description
    LocalMax-mean differencing scheme class.

    This scheme interpolates 1/field using a scheme specified at run-time
    and return the reciprocal of the interpolate.

SourceFiles
    localMax.C

\*---------------------------------------------------------------------------*/

#ifndef localMax_H
#define localMax_H

#include "surfaceInterpolationScheme.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class localMax Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class localMax
:
    public surfaceInterpolationScheme<Type>
{
    // Private Member Functions

        //- Disallow default bitwise assignment
        void operator=(const localMax&) = delete;


public:

    //- Runtime type information
    TypeName("localMax");


    // Constructors

        //- Construct from mesh
        localMax(const fvMesh& mesh)
        :
            surfaceInterpolationScheme<Type>(mesh)
        {}

        //- Construct from Istream.
        //  The name of the flux field is read from the Istream and looked-up
        //  from the mesh objectRegistry
        localMax
        (
            const fvMesh& mesh,
            Istream& is
        )
        :
            surfaceInterpolationScheme<Type>(mesh)
        {}

        //- Construct from faceFlux and Istream
        localMax
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            Istream& is
        )
        :
            surfaceInterpolationScheme<Type>(mesh)
        {}


    // Member Functions

        //- Return the interpolation weighting factors
        virtual tmp<surfaceScalarField> weights
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const
        {
            NotImplemented;

            return tmp<surfaceScalarField>(nullptr);
        }

        //- Return the face-interpolate of the given cell field
        virtual tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
        interpolate
        (
            const GeometricField<Type, fvPatchField, volMesh>& vf
        ) const
        {
            const fvMesh& mesh = vf.mesh();

            tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> tvff
            (
                new GeometricField<Type, fvsPatchField, surfaceMesh>
                (
                    IOobject
                    (
                        "localMax::interpolate(" + vf.name() + ')',
                        mesh.time().timeName(),
                        mesh
                    ),
                    mesh,
                    vf.dimensions()
                )
            );
            GeometricField<Type, fvsPatchField, surfaceMesh>& vff = tvff.ref();

            const labelUList& own = mesh.owner();
            const labelUList& nei = mesh.neighbour();

            forAll(vff, facei)
            {
                vff[facei] = max(vf[own[facei]], vf[nei[facei]]);
            }

            typename GeometricField<Type, fvsPatchField, surfaceMesh>::
                Boundary& bff = vff.boundaryFieldRef();

            forAll(bff, patchi)
            {
                const fvPatchField<Type>& pf = vf.boundaryField()[patchi];
                Field<Type>& pff = bff[patchi];

                if (pf.coupled())
                {
                    tmp<Field<Type>> tpif(pf.patchInternalField());
                    const Field<Type>& pif = tpif();

                    tmp<Field<Type>> tpnf(pf.patchNeighbourField());
                    const Field<Type>& pnf = tpnf();

                    forAll(pff, i)
                    {
                        pff[i] = max(pif[i], pnf[i]);
                    }
                }
                else
                {
                    pff = pf;
                }
            }

            return tvff;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
