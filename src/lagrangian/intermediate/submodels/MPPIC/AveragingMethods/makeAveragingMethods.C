/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2014 OpenFOAM Foundation
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

#include "Field.H"
#include "fvcGrad.H"
#include "polyMeshTetDecomposition.H"

#include "Basic.H"
#include "Dual.H"
#include "Moment.H"

namespace Foam
{
    // Scalar interpolation
    defineNamedTemplateTypeNameAndDebug(AveragingMethod<scalar>, 0);
    defineTemplateRunTimeSelectionTable
    (
        AveragingMethod<scalar>,
        dictionary
    );

    // Vector interpolation
    defineNamedTemplateTypeNameAndDebug(AveragingMethod<vector>, 0);
    defineTemplateRunTimeSelectionTable
    (
        AveragingMethod<vector>,
        dictionary
    );

    namespace AveragingMethods
    {
        // Basic interpolation
        defineNamedTemplateTypeNameAndDebug(Basic<scalar>, 0);
        AveragingMethod<scalar>::
            adddictionaryConstructorToTable<Basic<scalar> >
            addBasicscalarConstructorToTable_;

        defineNamedTemplateTypeNameAndDebug(Basic<vector>, 0);
        AveragingMethod<vector>::
            adddictionaryConstructorToTable<Basic<vector> >
            addBasicvectorConstructorToTable_;

        // Dual interpolation
        defineNamedTemplateTypeNameAndDebug(Dual<scalar>, 0);
        AveragingMethod<scalar>::
            adddictionaryConstructorToTable<Dual<scalar> >
            addDualscalarConstructorToTable_;

        defineNamedTemplateTypeNameAndDebug(Dual<vector>, 0);
        AveragingMethod<vector>::
            adddictionaryConstructorToTable<Dual<vector> >
            addDualvectorConstructorToTable_;

        // Moment interpolation
        defineNamedTemplateTypeNameAndDebug(Moment<scalar>, 0);
        AveragingMethod<scalar>::
            adddictionaryConstructorToTable<Moment<scalar> >
            addMomentscalarConstructorToTable_;

        defineNamedTemplateTypeNameAndDebug(Moment<vector>, 0);
        AveragingMethod<vector>::
            adddictionaryConstructorToTable<Moment<vector> >
            addMomentvectorConstructorToTable_;
    }
}


// ************************************************************************* //
