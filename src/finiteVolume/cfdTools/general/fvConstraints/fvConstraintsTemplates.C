/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021 OpenFOAM Foundation
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

// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

template<class Type>
bool Foam::fvConstraints::constrain(fvMatrix<Type>& eqn) const
{
    checkApplied();

    const PtrListDictionary<fvConstraint>& constraintList(*this);

    bool constrained = false;

    forAll(constraintList, i)
    {
        const fvConstraint& constraint = constraintList[i];

        if (constraint.constrainsField(eqn.psi().name()))
        {
            constrainedFields_[i].insert(eqn.psi().name());

            if (debug)
            {
                Info<< "Applying constraint " << constraint.name()
                    << " to field " << eqn.psi().name() << endl;
            }

            constrained =
                constraint.constrain(eqn, eqn.psi().name()) || constrained;
        }
    }

    return constrained;
}


template<class Type>
bool Foam::fvConstraints::constrain
(
    GeometricField<Type, fvPatchField, volMesh>& field
) const
{
    const word& fieldName = field.name();

    const PtrListDictionary<fvConstraint>& constraintList(*this);

    bool constrained = false;

    forAll(constraintList, i)
    {
        const fvConstraint& constraint = constraintList[i];

        if (constraint.constrainsField(fieldName))
        {
            constrainedFields_[i].insert(fieldName);

            if (debug)
            {
                Info<< "Applying constraint " << constraint.name()
                    << " for field " << fieldName << endl;
            }

            constrained =
                constraint.constrain(field) || constrained;
        }
    }

    return constrained;
}


// ************************************************************************* //
