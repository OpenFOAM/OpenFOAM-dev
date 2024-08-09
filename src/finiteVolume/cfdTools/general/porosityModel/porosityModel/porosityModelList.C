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

#include "porosityModelList.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::porosityModelList::porosityModelList
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    PtrList<porosityModel>(),
    mesh_(mesh)
{
    reset(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::porosityModelList::~porosityModelList()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::porosityModelList::reset(const dictionary& dict)
{
    label count = 0;
    forAllConstIter(dictionary, dict, iter)
    {
        if (iter().isDict())
        {
            count++;
        }
    }

    this->setSize(count);
    label i = 0;
    forAllConstIter(dictionary, dict, iter)
    {
        if (iter().isDict())
        {
            const word& name = iter().keyword();
            const dictionary& modelDict = iter().dict();

            this->set
            (
                i++,
                porosityModel::New(name, mesh_, modelDict)
            );
        }
    }
}


void Foam::porosityModelList::addResistance
(
    fvVectorMatrix& UEqn
)
{
    forAll(*this, i)
    {
        this->operator[](i).addResistance(UEqn);
    }
}


void Foam::porosityModelList::addResistance
(
    const fvVectorMatrix& UEqn,
    volTensorField& AU,
    bool correctAUprocBC
)
{
    forAll(*this, i)
    {
        this->operator[](i).addResistance(UEqn, AU, correctAUprocBC);
    }
}


bool Foam::porosityModelList::movePoints()
{
    forAll(*this, i)
    {
        this->operator[](i).movePoints();
    }

    return true;
}


void Foam::porosityModelList::topoChange(const polyTopoChangeMap& map)
{
    forAll(*this, i)
    {
        this->operator[](i).topoChange(map);
    }
}


void Foam::porosityModelList::mapMesh(const polyMeshMap& map)
{
    forAll(*this, i)
    {
        this->operator[](i).mapMesh(map);
    }
}


void Foam::porosityModelList::distribute
(
    const polyDistributionMap& map
)
{
    forAll(*this, i)
    {
        this->operator[](i).distribute(map);
    }
}


bool Foam::porosityModelList::read(const dictionary& dict)
{
    bool allOk = true;
    forAll(*this, i)
    {
        porosityModel& pm = this->operator[](i);
        bool ok = pm.read(dict.subDict(pm.name()));
        allOk = (allOk && ok);
    }
    return allOk;
}


// ************************************************************************* //
