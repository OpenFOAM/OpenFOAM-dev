/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2019 OpenFOAM Foundation
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

#include "residuals.H"
#include "volFields.H"
#include "Residuals.H"
#include "ListOps.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void Foam::functionObjects::residuals::writeFileHeader(const word& fieldName)
{
    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;

    if (obr_.foundObject<fieldType>(fieldName))
    {
        typename pTraits<Type>::labelType validComponents
        (
            mesh_.validComponents<Type>()
        );

        for (direction cmpt=0; cmpt<pTraits<Type>::nComponents; cmpt++)
        {
            if (component(validComponents, cmpt) != -1)
            {
                writeTabbed
                (
                    file(),
                    fieldName + word(pTraits<Type>::componentNames[cmpt])
                );
            }
        }
    }
}


template<class Type>
void Foam::functionObjects::residuals::writeResidual(const word& fieldName)
{
    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;

    if (obr_.foundObject<fieldType>(fieldName))
    {
        if (Residuals<Type>::found(mesh_, fieldName))
        {
            const DynamicList<SolverPerformance<Type>>& sp
            (
                Residuals<Type>::field(mesh_, fieldName)
            );

            const Type& residual = sp.first().initialResidual();

            typename pTraits<Type>::labelType validComponents
            (
                mesh_.validComponents<Type>()
            );

            for (direction cmpt=0; cmpt<pTraits<Type>::nComponents; cmpt++)
            {
                if (component(validComponents, cmpt) != -1)
                {
                    file() << tab << component(residual, cmpt);
                }
            }
        }
        else
        {
            file() << tab << "N/A";
        }
    }
}


// ************************************************************************* //
