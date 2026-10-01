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

\*---------------------------------------------------------------------------*/

#include "interfaceCompressionCoefficient.H"
#include "generateInterfacialModels.H"

// * * * * * * * * * * * * * * * * Selector  * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::interfaceCompressionCoefficient>
Foam::interfaceCompressionCoefficient::New
(
    const UPtrList<const entry>& entries,
    const phaseInterface& interface
)
{
    const Foam::entry& entry =
        modelEntry<interfaceCompressionCoefficient>(entries);

    return
        autoPtr<interfaceCompressionCoefficient>
        (
            new interfaceCompressionCoefficient
            (
                IOobject::groupName(typeName, interface.name()),
                dimless,
                entry.stream()
            )
        );
}


// ************************************************************************* //
