/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

Description
    Abstract base class for surface interpolation schemes.

\*---------------------------------------------------------------------------*/

#include "surfaceInterpolationScheme.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "coupledFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

// Return weighting factors for scheme given by name in dictionary
template<class Type>
tmp<surfaceInterpolationScheme<Type> > surfaceInterpolationScheme<Type>::New
(
    const fvMesh& mesh,
    Istream& schemeData
)
{
    if (schemeData.eof())
    {
        FatalIOErrorIn
        (
            "surfaceInterpolationScheme<Type>::New(const fvMesh&, Istream&)",
            schemeData
        )   << "Discretisation scheme not specified"
            << endl << endl
            << "Valid schemes are :" << endl
            << MeshConstructorTablePtr_->sortedToc()
            << exit(FatalIOError);
    }

    const word schemeName(schemeData);

    if (surfaceInterpolation::debug || surfaceInterpolationScheme<Type>::debug)
    {
        Info<< "surfaceInterpolationScheme<Type>::New"
               "(const fvMesh&, Istream&)"
               " : discretisation scheme = "
            << schemeName
            << endl;
    }

    typename MeshConstructorTable::iterator constructorIter =
        MeshConstructorTablePtr_->find(schemeName);

    if (constructorIter == MeshConstructorTablePtr_->end())
    {
        FatalIOErrorIn
        (
            "surfaceInterpolationScheme<Type>::New(const fvMesh&, Istream&)",
            schemeData
        )   << "Unknown discretisation scheme "
            << schemeName << nl << nl
            << "Valid schemes are :" << endl
            << MeshConstructorTablePtr_->sortedToc()
            << exit(FatalIOError);
    }

    return constructorIter()(mesh, schemeData);
}


// Return weighting factors for scheme given by name in dictionary
template<class Type>
tmp<surfaceInterpolationScheme<Type> > surfaceInterpolationScheme<Type>::New
(
    const fvMesh& mesh,
    const surfaceScalarField& faceFlux,
    Istream& schemeData
)
{
    if (schemeData.eof())
    {
        FatalIOErrorIn
        (
            "surfaceInterpolationScheme<Type>::New"
            "(const fvMesh&, const surfaceScalarField&, Istream&)",
            schemeData
        )   << "Discretisation scheme not specified"
            << endl << endl
            << "Valid schemes are :" << endl
            << MeshConstructorTablePtr_->sortedToc()
            << exit(FatalIOError);
    }

    const word schemeName(schemeData);

    if (surfaceInterpolation::debug || surfaceInterpolationScheme<Type>::debug)
    {
        Info<< "surfaceInterpolationScheme<Type>::New"
               "(const fvMesh&, const surfaceScalarField&, Istream&)"
               " : discretisation scheme = "
            << schemeName
            << endl;
    }

    typename MeshFluxConstructorTable::iterator constructorIter =
        MeshFluxConstructorTablePtr_->find(schemeName);

    if (constructorIter == MeshFluxConstructorTablePtr_->end())
    {
        FatalIOErrorIn
        (
            "surfaceInterpolationScheme<Type>::New"
            "(const fvMesh&, const surfaceScalarField&, Istream&)",
            schemeData
        )   << "Unknown discretisation scheme "
            << schemeName << nl << nl
            << "Valid schemes are :" << endl
            << MeshFluxConstructorTablePtr_->sortedToc()
            << exit(FatalIOError);
    }

    return constructorIter()(mesh, faceFlux, schemeData);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
surfaceInterpolationScheme<Type>::~surfaceInterpolationScheme()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- Return the face-interpolate of the given cell field
//  with the given owner and neighbour weighting factors
template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >
surfaceInterpolationScheme<Type>::interpolate
(
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const tmp<surfaceScalarField>& tlambdas,
    const tmp<surfaceScalarField>& tys
)
{
    if (surfaceInterpolation::debug)
    {
        Info<< "surfaceInterpolationScheme<Type>::uncorrectedInterpolate"
               "(const GeometricField<Type, fvPatchField, volMesh>&, "
               "const tmp<surfaceScalarField>&, "
               "const tmp<surfaceScalarField>&) : "
               "interpolating "
            << vf.type() << " "
            << vf.name()
            << " from cells to faces "
               "without explicit correction"
            << endl;
    }

    const surfaceScalarField& lambdas = tlambdas();
    const surfaceScalarField& ys = tys();

    const Field<Type>& vfi = vf.internalField();
    const scalarField& lambda = lambdas.internalField();
    const scalarField& y = ys.internalField();

    const fvMesh& mesh = vf.mesh();
    const labelUList& P = mesh.owner();
    const labelUList& N = mesh.neighbour();

    tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > tsf
    (
        new GeometricField<Type, fvsPatchField, surfaceMesh>
        (
            IOobject
            (
                "interpolate("+vf.name()+')',
                vf.instance(),
                vf.db()
            ),
            mesh,
            vf.dimensions()
        )
    );
    GeometricField<Type, fvsPatchField, surfaceMesh>& sf = tsf();

    Field<Type>& sfi = sf.internalField();

    for (register label fi=0; fi<P.size(); fi++)
    {
        sfi[fi] = lambda[fi]*vfi[P[fi]] + y[fi]*vfi[N[fi]];
    }


    // Interpolate across coupled patches using given lambdas and ys

    forAll(lambdas.boundaryField(), pi)
    {
        const fvsPatchScalarField& pLambda = lambdas.boundaryField()[pi];
        const fvsPatchScalarField& pY = ys.boundaryField()[pi];

        if (vf.boundaryField()[pi].coupled())
        {
            sf.boundaryField()[pi] =
                pLambda*vf.boundaryField()[pi].patchInternalField()
              + pY*vf.boundaryField()[pi].patchNeighbourField();
        }
        else
        {
            sf.boundaryField()[pi] = vf.boundaryField()[pi];
        }
    }

    tlambdas.clear();
    tys.clear();

    return tsf;
}


//- Return the face-interpolate of the given cell field
//  with the given weighting factors
template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >
surfaceInterpolationScheme<Type>::interpolate
(
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const tmp<surfaceScalarField>& tlambdas
)
{
    if (surfaceInterpolation::debug)
    {
        Info<< "surfaceInterpolationScheme<Type>::interpolate"
               "(const GeometricField<Type, fvPatchField, volMesh>&, "
               "const tmp<surfaceScalarField>&) : "
               "interpolating "
            << vf.type() << " "
            << vf.name()
            << " from cells to faces "
               "without explicit correction"
            << endl;
    }

    const surfaceScalarField& lambdas = tlambdas();

    const Field<Type>& vfi = vf.internalField();
    const scalarField& lambda = lambdas.internalField();

    const fvMesh& mesh = vf.mesh();
    const labelUList& P = mesh.owner();
    const labelUList& N = mesh.neighbour();

    tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > tsf
    (
        new GeometricField<Type, fvsPatchField, surfaceMesh>
        (
            IOobject
            (
                "interpolate("+vf.name()+')',
                vf.instance(),
                vf.db()
            ),
            mesh,
            vf.dimensions()
        )
    );
    GeometricField<Type, fvsPatchField, surfaceMesh>& sf = tsf();

    Field<Type>& sfi = sf.internalField();

    for (register label fi=0; fi<P.size(); fi++)
    {
        sfi[fi] = lambda[fi]*(vfi[P[fi]] - vfi[N[fi]]) + vfi[N[fi]];
    }

    // Interpolate across coupled patches using given lambdas

    forAll(lambdas.boundaryField(), pi)
    {
        const fvsPatchScalarField& pLambda = lambdas.boundaryField()[pi];

        if (vf.boundaryField()[pi].coupled())
        {
            tsf().boundaryField()[pi] =
                pLambda*vf.boundaryField()[pi].patchInternalField()
             + (1.0 - pLambda)*vf.boundaryField()[pi].patchNeighbourField();
        }
        else
        {
            sf.boundaryField()[pi] = vf.boundaryField()[pi];
        }
    }

    tlambdas.clear();

    return tsf;
}


//- Return the face-interpolate of the given cell field
//  with explicit correction
template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >
surfaceInterpolationScheme<Type>::interpolate
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    if (surfaceInterpolation::debug)
    {
        Info<< "surfaceInterpolationScheme<Type>::interpolate"
               "(const GeometricField<Type, fvPatchField, volMesh>&) : "
               "interpolating "
            << vf.type() << " "
            << vf.name()
            << " from cells to faces"
            << endl;
    }

    tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > tsf
        = interpolate(vf, weights(vf));

    if (corrected())
    {
        tsf() += correction(vf);
    }

    return tsf;
}


//- Return the face-interpolate of the given cell field
//  with explicit correction
template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >
surfaceInterpolationScheme<Type>::interpolate
(
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tvf
) const
{
    tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > tinterpVf
        = interpolate(tvf());
    tvf.clear();
    return tinterpVf;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
