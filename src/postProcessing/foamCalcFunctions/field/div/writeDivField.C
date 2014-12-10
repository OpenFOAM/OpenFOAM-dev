/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

template<class Type>
void Foam::calcTypes::div::writeDivField
(
    const IOobject& header,
    const fvMesh& mesh,
    bool& processed
)
{
    if (header.headerClassName() == Type::typeName)
    {
        Info<< "    Reading " << header.name() << endl;
        Type field(header, mesh);

        Info<< "    Calculating div" << header.name() << endl;
        volScalarField divField
        (
            IOobject
            (
                "div" + header.name(),
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ
            ),
            fvc::div(field)
        );
        divField.write();

        processed = true;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
