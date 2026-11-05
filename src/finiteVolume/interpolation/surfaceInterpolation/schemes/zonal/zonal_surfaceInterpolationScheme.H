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

Class
    Foam::zonal

Description
    Zonal discretisation scheme which interpolates the given \c volField to the
    mesh faces using a default scheme and then applies different schemes to the
    faces of the given set of \c faceZones.

    Examples

    Setting the default scheme to \c limitedLinearV and upwind in the \c
    faceZones fz1 and fz2:
    \verbatim
    div(phi,U)      Gauss zonal limitedLinearV 1
    {
        fz1     upwind;
        fz2     upwind;
    };
    \endverbatim

    The usual wildcards and selection sets are supported, for example to
    select the upwind scheme in all \c faceZones:
    \verbatim
    div(phi,U)      Gauss zonal linearUpwind grad(U)
    {
        ".*"    upwind;
    };
    \endverbatim

SourceFiles
    zonal_surfaceInterpolationScheme.C

\*---------------------------------------------------------------------------*/

#ifndef zonal_surfaceInterpolationScheme_H
#define zonal_surfaceInterpolationScheme_H

#include "surfaceInterpolationScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class zonal Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class zonal
:
    public surfaceInterpolationScheme<Type>
{
    // Private Data

        //- Face-flux field used by the base and zone schemes
        const surfaceScalarField& faceFlux_;

        //- The base interpolation scheme
        tmp<surfaceInterpolationScheme<Type>> tScheme_;

        //- Dictionary of faceZone schemes
        dictionary zoneSchemesDict_;

        //- List of instantiated faceZone schemes
        PtrList<surfaceInterpolationScheme<Type>> zoneSchemes_;


public:

    //- Runtime type information
    TypeName("zonal");


    // Constructors

        //- Construct from mesh, faceFlux and Istream
        zonal
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            Istream& schemeData
        );

        //- Construct from mesh and Istream
        zonal(const fvMesh& mesh, Istream& schemeData);

        //- Disallow default bitwise copy construction
        zonal(const zonal&) = delete;


    // Member Functions

        //- Return the interpolation weighting factors
        tmp<surfaceScalarField> weights(const VolField<Type>& vf) const;

        //- Return true if this scheme uses an explicit correction
        virtual bool corrected() const;

        //- Return the explicit correction to the face-interpolate
        //  set to zero on the near-boundary faces where upwind is applied
        virtual tmp<SurfaceField<Type>>
        correction(const VolField<Type>& vf) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const zonal&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
