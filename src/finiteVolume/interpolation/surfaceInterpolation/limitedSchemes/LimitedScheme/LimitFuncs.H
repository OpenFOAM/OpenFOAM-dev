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

Namespace
    Foam::limitFuncs

Description
    Namespace for limiting functions


Class
    Foam::limitFuncs::LimitFuncs

Description
    Class to create NVD/TVD limited weighting-factors.

    The particular differencing scheme class is supplied as a template
    argument, the weight function of which is called by the weight function
    of this class for the internal faces as well as faces of coupled
    patches (e.g. processor-processor patches). The weight function is
    supplied the central-differencing weighting factor, the face-flux, the
    cell and face gradients (from which the normalised variable
    distribution may be created) and the cell centre distance.

    This code organisation is both neat and efficient, allowing for
    convenient implementation of new schemes to run on parallelised cases.

SourceFiles
    LimitFuncs.C

\*---------------------------------------------------------------------------*/

#ifndef LimitFuncs_H
#define LimitFuncs_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace limitFuncs
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
class null
{
public:

    null()
    {}

    inline tmp<GeometricField<Type, fvPatchField, volMesh>> operator()
    (
        const GeometricField<Type, fvPatchField, volMesh>& phi
    ) const
    {
        return phi;
    }
};


template<class Type>
class magSqr
{
public:

    magSqr()
    {}

    inline tmp<volScalarField> operator()
    (
        const GeometricField<Type, fvPatchField, volMesh>&
    ) const;
};

template<>
inline tmp<volScalarField> magSqr<scalar>::operator()
(
    const volScalarField& phi
) const;

template<>
inline tmp<volScalarField> magSqr<symmTensor>::operator()
(
    const volSymmTensorField& phi
) const;

template<>
inline tmp<volScalarField> magSqr<tensor>::operator()
(
    const volTensorField& phi
) const;


template<class Type>
class rhoMagSqr
{
public:

    rhoMagSqr()
    {}

    inline tmp<volScalarField> operator()
    (
        const GeometricField<Type, fvPatchField, volMesh>&
    ) const;
};

template<>
inline tmp<volScalarField> rhoMagSqr<scalar>::operator()
(
    const volScalarField& phi
) const;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace limitFuncs
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "LimitFuncs.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
