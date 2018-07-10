/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2018 OpenFOAM Foundation
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

#include "PhaseLimitStabilization.H"
#include "fvMatrices.H"
#include "fvmSup.H"
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::fv::PhaseLimitStabilization<Type>::PhaseLimitStabilization
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    option(name, modelType, dict, mesh),
    fieldName_(coeffs_.lookup("field")),
    rateName_(coeffs_.lookup("rate")),
    residualAlpha_(readScalar(coeffs_.lookup("residualAlpha")))
{
    fieldNames_.setSize(1, fieldName_);
    applied_.setSize(1, false);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::fv::PhaseLimitStabilization<Type>::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    fvMatrix<Type>& eqn,
    const label fieldi
)
{
    const GeometricField<Type, fvPatchField, volMesh>& psi = eqn.psi();

    uniformDimensionedScalarField& rate =
        mesh_.lookupObjectRef<uniformDimensionedScalarField>(rateName_);

    eqn -= fvm::Sp(max(residualAlpha_ - alpha, scalar(0))*rho*rate, psi);
}


template<class Type>
bool Foam::fv::PhaseLimitStabilization<Type>::read(const dictionary& dict)
{
    if (option::read(dict))
    {
        coeffs_.lookup("residualAlpha") >> residualAlpha_;

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
