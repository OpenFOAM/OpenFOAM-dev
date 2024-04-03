/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2024 OpenFOAM Foundation
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

#include "vtkUnstructuredReader.H"
#include "labelIOField.H"
#include "scalarIOField.H"
#include "stringIOList.H"
#include "cellModeller.H"
#include "vectorIOField.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class T>
void Foam::vtkUnstructuredReader::readBlock
(
    Istream& inFile,
    const label n,
    List<T>& lst
) const
{
    lst.setSize(n);
    forAll(lst, i)
    {
        inFile >> lst[i];
    }
}


template<class Type>
void Foam::vtkUnstructuredReader::printFieldStats
(
    const objectRegistry& obj
) const
{
    wordList fieldNames(obj.toc(Type::typeName));

    if (fieldNames.size() > 0)
    {
        Info<< "Read " << fieldNames.size() << " " << Type::typeName
            << " fields:" << endl;
        Info<< "Size\tName" << nl
            << "----\t----" << endl;
        forAll(fieldNames, i)
        {
            Info<< obj.lookupObject<Type>(fieldNames[i]).size()
                << "\t" << fieldNames[i]
                << endl;
        }
        Info<< endl;
    }
}


// ************************************************************************* //
