/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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

#include "carrierLagrangianFieldSource.H"
#include "coupled.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::carrierLagrangianFieldSource<Type>::carrierLagrangianFieldSource
(
    const regIOobject& iIo,
    const dictionary& dict
)
:
    LagrangianFieldSource<Type>(iIo, dict),
    CloudLagrangianFieldSource<Type>(*this),
    fieldcName_(dict.lookupOrDefault<word>("fieldc", iIo.name()))
{}


template<class Type>
Foam::carrierLagrangianFieldSource<Type>::carrierLagrangianFieldSource
(
    const carrierLagrangianFieldSource& field,
    const regIOobject& iIo
)
:
    LagrangianFieldSource<Type>(field, iIo),
    CloudLagrangianFieldSource<Type>(*this),
    fieldcName_(field.fieldcName_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::carrierLagrangianFieldSource<Type>::~carrierLagrangianFieldSource()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::LagrangianSubField<Type>>
Foam::carrierLagrangianFieldSource<Type>::value
(
    const LagrangianInjection& injection,
    const LagrangianSubMesh& subMesh
) const
{
    return
        this->template cloud<clouds::coupled>(injection, subMesh).carrierField
        (
            subMesh.mesh().mesh().lookupObject<VolField<Type>>(fieldcName_)
        )(subMesh);
}


template<class Type>
void Foam::carrierLagrangianFieldSource<Type>::write
(
    Ostream& os
) const
{
    LagrangianFieldSource<Type>::write(os);

    writeEntryIfDifferent<word>
    (
        os,
        "fieldc",
        this->internalField().name(),
        fieldcName_
    );
}


// ************************************************************************* //
