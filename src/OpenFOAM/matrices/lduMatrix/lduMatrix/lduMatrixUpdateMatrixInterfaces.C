/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
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
        Pstream::defaultCommsType == Pstream::commsTypes::blocking
     || Pstream::defaultCommsType == Pstream::commsTypes::nonBlocking
    )
    {
        forAll(interfaces, interfacei)
        {
            if (interfaces.set(interfacei))
            {
                interfaces[interfacei].initInterfaceMatrixUpdate
                (
                    result,
                    psiif,
                    coupleCoeffs[interfacei],
                    cmpt,
                    Pstream::defaultCommsType
                );
            }
        }
    }
    else if (Pstream::defaultCommsType == Pstream::commsTypes::scheduled)
    {
        const lduSchedule& patchSchedule = this->patchSchedule();

        // Loop over the "global" patches are on the list of interfaces but
        // beyond the end of the schedule which only handles "normal" patches
        for
        (
            label interfacei=patchSchedule.size()/2;
            interfacei<interfaces.size();
            interfacei++
        )
        {
            if (interfaces.set(interfacei))
            {
                interfaces[interfacei].initInterfaceMatrixUpdate
                (
                    result,
                    psiif,
                    coupleCoeffs[interfacei],
                    cmpt,
                    Pstream::commsTypes::blocking
                );
            }
        }
    }
    else
    {
        FatalErrorInFunction
            << "Unsupported communications type "
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
    if (Pstream::defaultCommsType == Pstream::commsTypes::blocking)
    {
        forAll(interfaces, interfacei)
        {
            if (interfaces.set(interfacei))
            {
                interfaces[interfacei].updateInterfaceMatrix
                (
                    result,
                    psiif,
                    coupleCoeffs[interfacei],
                    cmpt,
                    Pstream::defaultCommsType
                );
            }
        }
    }
    else if (Pstream::defaultCommsType == Pstream::commsTypes::nonBlocking)
    {
        // Try and consume interfaces as they become available
        bool allUpdated = false;

        for (label i=0; i<UPstream::nPollProcInterfaces; i++)
        {
            allUpdated = true;

            forAll(interfaces, interfacei)
            {
                if (interfaces.set(interfacei))
                {
                    if (!interfaces[interfacei].updatedMatrix())
                    {
                        if (interfaces[interfacei].ready())
                        {
                            interfaces[interfacei].updateInterfaceMatrix
                            (
                                result,
                                psiif,
                                coupleCoeffs[interfacei],
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
        forAll(interfaces, interfacei)
        {
            if
            (
                interfaces.set(interfacei)
            && !interfaces[interfacei].updatedMatrix()
            )
            {
                interfaces[interfacei].updateInterfaceMatrix
                (
                    result,
                    psiif,
                    coupleCoeffs[interfacei],
                    cmpt,
                    Pstream::defaultCommsType
                );
            }
        }
    }
    else if (Pstream::defaultCommsType == Pstream::commsTypes::scheduled)
    {
        const lduSchedule& patchSchedule = this->patchSchedule();

        // Loop over all the "normal" interfaces relating to standard patches
        forAll(patchSchedule, i)
        {
            label interfacei = patchSchedule[i].patch;

            if (interfaces.set(interfacei))
            {
                if (patchSchedule[i].init)
                {
                    interfaces[interfacei].initInterfaceMatrixUpdate
                    (
                        result,
                        psiif,
                        coupleCoeffs[interfacei],
                        cmpt,
                        Pstream::commsTypes::scheduled
                    );
                }
                else
                {
                    interfaces[interfacei].updateInterfaceMatrix
                    (
                        result,
                        psiif,
                        coupleCoeffs[interfacei],
                        cmpt,
                        Pstream::commsTypes::scheduled
                    );
                }
            }
        }

        // Loop over the "global" patches are on the list of interfaces but
        // beyond the end of the schedule which only handles "normal" patches
        for
        (
            label interfacei=patchSchedule.size()/2;
            interfacei<interfaces.size();
            interfacei++
        )
        {
            if (interfaces.set(interfacei))
            {
                interfaces[interfacei].updateInterfaceMatrix
                (
                    result,
                    psiif,
                    coupleCoeffs[interfacei],
                    cmpt,
                    Pstream::commsTypes::blocking
                );
            }
        }
    }
    else
    {
        FatalErrorInFunction
            << "Unsupported communications type "
            << Pstream::commsTypeNames[Pstream::defaultCommsType]
            << exit(FatalError);
    }
}


// ************************************************************************* //
