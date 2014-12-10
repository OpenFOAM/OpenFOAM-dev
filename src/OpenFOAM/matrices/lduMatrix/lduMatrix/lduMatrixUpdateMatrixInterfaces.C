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

\*---------------------------------------------------------------------------*/

#include "lduMatrix.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::lduMatrix::initMatrixInterfaces
(
    const FieldField<Field, scalar>& coupleCoeffs,
    const lduInterfaceFieldPtrsList& interfaces,
    const scalarField& psiif,
    scalarField& result,
    const direction cmpt
) const
{
    if
    (
        Pstream::defaultCommsType == Pstream::blocking
     || Pstream::defaultCommsType == Pstream::nonBlocking
    )
    {
        forAll(interfaces, interfaceI)
        {
            if (interfaces.set(interfaceI))
            {
                interfaces[interfaceI].initInterfaceMatrixUpdate
                (
                    result,
                    psiif,
                    coupleCoeffs[interfaceI],
                    cmpt,
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
            interfaceI<interfaces.size();
            interfaceI++
        )
        {
            if (interfaces.set(interfaceI))
            {
                interfaces[interfaceI].initInterfaceMatrixUpdate
                (
                    result,
                    psiif,
                    coupleCoeffs[interfaceI],
                    cmpt,
                    Pstream::blocking
                );
            }
        }
    }
    else
    {
        FatalErrorIn("lduMatrix::initMatrixInterfaces(..)")
            << "Unsuported communications type "
            << Pstream::commsTypeNames[Pstream::defaultCommsType]
            << exit(FatalError);
    }
}


void Foam::lduMatrix::updateMatrixInterfaces
(
    const FieldField<Field, scalar>& coupleCoeffs,
    const lduInterfaceFieldPtrsList& interfaces,
    const scalarField& psiif,
    scalarField& result,
    const direction cmpt
) const
{
    if (Pstream::defaultCommsType == Pstream::blocking)
    {
        forAll(interfaces, interfaceI)
        {
            if (interfaces.set(interfaceI))
            {
                interfaces[interfaceI].updateInterfaceMatrix
                (
                    result,
                    psiif,
                    coupleCoeffs[interfaceI],
                    cmpt,
                    Pstream::defaultCommsType
                );
            }
        }
    }
    else if (Pstream::defaultCommsType == Pstream::nonBlocking)
    {
        // Try and consume interfaces as they become available
        bool allUpdated = false;

        for (label i = 0; i < UPstream::nPollProcInterfaces; i++)
        {
            allUpdated = true;

            forAll(interfaces, interfaceI)
            {
                if (interfaces.set(interfaceI))
                {
                    if (!interfaces[interfaceI].updatedMatrix())
                    {
                        if (interfaces[interfaceI].ready())
                        {
                            interfaces[interfaceI].updateInterfaceMatrix
                            (
                                result,
                                psiif,
                                coupleCoeffs[interfaceI],
                                cmpt,
                                Pstream::defaultCommsType
                            );
                        }
                        else
                        {
                            allUpdated = false;
                        }
                    }
                }
            }

            if (allUpdated)
            {
                break;
            }
        }

        // Block for everything
        if (Pstream::parRun())
        {
            if (allUpdated)
            {
                // All received. Just remove all storage of requests
                // Note that we don't know what starting number of requests
                // was before start of sends and receives (since set from
                // initMatrixInterfaces) so set to 0 and loose any in-flight
                // requests.
                UPstream::resetRequests(0);
            }
            else
            {
                // Block for all requests and remove storage
                UPstream::waitRequests();
            }
        }

        // Consume
        forAll(interfaces, interfaceI)
        {
            if
            (
                interfaces.set(interfaceI)
            && !interfaces[interfaceI].updatedMatrix()
            )
            {
                interfaces[interfaceI].updateInterfaceMatrix
                (
                    result,
                    psiif,
                    coupleCoeffs[interfaceI],
                    cmpt,
                    Pstream::defaultCommsType
                );
            }
        }
    }
    else if (Pstream::defaultCommsType == Pstream::scheduled)
    {
        const lduSchedule& patchSchedule = this->patchSchedule();

        // Loop over all the "normal" interfaces relating to standard patches
        forAll(patchSchedule, i)
        {
            label interfaceI = patchSchedule[i].patch;

            if (interfaces.set(interfaceI))
            {
                if (patchSchedule[i].init)
                {
                    interfaces[interfaceI].initInterfaceMatrixUpdate
                    (
                        result,
                        psiif,
                        coupleCoeffs[interfaceI],
                        cmpt,
                        Pstream::scheduled
                    );
                }
                else
                {
                    interfaces[interfaceI].updateInterfaceMatrix
                    (
                        result,
                        psiif,
                        coupleCoeffs[interfaceI],
                        cmpt,
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
            interfaceI<interfaces.size();
            interfaceI++
        )
        {
            if (interfaces.set(interfaceI))
            {
                interfaces[interfaceI].updateInterfaceMatrix
                (
                    result,
                    psiif,
                    coupleCoeffs[interfaceI],
                    cmpt,
                    Pstream::blocking
                );
            }
        }
    }
    else
    {
        FatalErrorIn("lduMatrix::updateMatrixInterfaces(..)")
            << "Unsuported communications type "
            << Pstream::commsTypeNames[Pstream::defaultCommsType]
            << exit(FatalError);
    }
}


// ************************************************************************* //
