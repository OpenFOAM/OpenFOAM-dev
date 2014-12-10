/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

#include "LduMatrix.H"
#include "lduInterfaceField.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type, class DType, class LUType>
void Foam::LduMatrix<Type, DType, LUType>::initMatrixInterfaces
(
    const FieldField<Field, LUType>& interfaceCoeffs,
    const Field<Type>& psiif,
    Field<Type>& result
) const
{
    if
    (
        Pstream::defaultCommsType == Pstream::blocking
     || Pstream::defaultCommsType == Pstream::nonBlocking
    )
    {
        forAll (interfaces_, interfaceI)
        {
            if (interfaces_.set(interfaceI))
            {
                interfaces_[interfaceI].initInterfaceMatrixUpdate
                (
                    result,
                    psiif,
                    interfaceCoeffs[interfaceI],
                    //Amultiplier<Type, LUType>(interfaceCoeffs[interfaceI]),
                    Pstream::defaultCommsType
                );
            }
        }
    }
    else if (Pstream::defaultCommsType == Pstream::scheduled)
    {
        const lduSchedule& patchSchedule = this->patchSchedule();

        // Loop over the "global" patches are on the list of interfaces but
        // beyond the end of the schedule which only handles "normal" patches
        for
        (
            label interfaceI=patchSchedule.size()/2;
            interfaceI<interfaces_.size();
            interfaceI++
        )
        {
            if (interfaces_.set(interfaceI))
            {
                interfaces_[interfaceI].initInterfaceMatrixUpdate
                (
                    result,
                    psiif,
                    interfaceCoeffs[interfaceI],
                    //Amultiplier<Type, LUType>(interfaceCoeffs[interfaceI]),
                    Pstream::blocking
                );
            }
        }
    }
    else
    {
        FatalErrorIn("LduMatrix<Type, DType, LUType>::initMatrixInterfaces")
            << "Unsuported communications type "
            << Pstream::commsTypeNames[Pstream::defaultCommsType]
            << exit(FatalError);
    }
}


template<class Type, class DType, class LUType>
void Foam::LduMatrix<Type, DType, LUType>::updateMatrixInterfaces
(
    const FieldField<Field, LUType>& interfaceCoeffs,
    const Field<Type>& psiif,
    Field<Type>& result
) const
{
    if
    (
        Pstream::defaultCommsType == Pstream::blocking
     || Pstream::defaultCommsType == Pstream::nonBlocking
    )
    {
        // Block until all sends/receives have been finished
        if (Pstream::defaultCommsType == Pstream::nonBlocking)
        {
            IPstream::waitRequests();
            OPstream::waitRequests();
        }

        forAll (interfaces_, interfaceI)
        {
            if (interfaces_.set(interfaceI))
            {
                interfaces_[interfaceI].updateInterfaceMatrix
                (
                    result,
                    psiif,
                    interfaceCoeffs[interfaceI],
                    //Amultiplier<Type, LUType>(interfaceCoeffs[interfaceI]),
                    Pstream::defaultCommsType
                );
            }
        }
    }
    else if (Pstream::defaultCommsType == Pstream::scheduled)
    {
        const lduSchedule& patchSchedule = this->patchSchedule();

        // Loop over all the "normal" interfaces relating to standard patches
        forAll (patchSchedule, i)
        {
            label interfaceI = patchSchedule[i].patch;

            if (interfaces_.set(interfaceI))
            {
                if (patchSchedule[i].init)
                {
                    interfaces_[interfaceI].initInterfaceMatrixUpdate
                    (
                        result,
                        psiif,
                        interfaceCoeffs[interfaceI],
                      //Amultiplier<Type, LUType>(interfaceCoeffs[interfaceI]),
                        Pstream::scheduled
                    );
                }
                else
                {
                    interfaces_[interfaceI].updateInterfaceMatrix
                    (
                        result,
                        psiif,
                        interfaceCoeffs[interfaceI],
                      //Amultiplier<Type, LUType>(interfaceCoeffs[interfaceI]),
                        Pstream::scheduled
                    );
                }
            }
        }

        // Loop over the "global" patches are on the list of interfaces but
        // beyond the end of the schedule which only handles "normal" patches
        for
        (
            label interfaceI=patchSchedule.size()/2;
            interfaceI<interfaces_.size();
            interfaceI++
        )
        {
            if (interfaces_.set(interfaceI))
            {
                interfaces_[interfaceI].updateInterfaceMatrix
                (
                    result,
                    psiif,
                    interfaceCoeffs[interfaceI],
                    //Amultiplier<Type, LUType>(interfaceCoeffs[interfaceI]),
                    Pstream::blocking
                );
            }
        }
    }
    else
    {
        FatalErrorIn("LduMatrix<Type, DType, LUType>::updateMatrixInterfaces")
            << "Unsuported communications type "
            << Pstream::commsTypeNames[Pstream::defaultCommsType]
            << exit(FatalError);
    }
}


// ************************************************************************* //
