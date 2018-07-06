/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
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

Description
    Constructor of cellModeller: just sets the cellModeller's params.

\*---------------------------------------------------------------------------*/

#include "cellModeller.H"
#include "IFstream.H"
#include "etcFiles.H"

// * * * * * * * * * * * * * * * Static data * * * * * * * * * * * * * * * * //

Foam::PtrList<Foam::cellModel> Foam::cellModeller::models_;

Foam::List<Foam::cellModel*> Foam::cellModeller::modelPtrs_;

Foam::HashTable<const Foam::cellModel*> Foam::cellModeller::modelDictionary_;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cellModeller::cellModeller()
{
    if (modelPtrs_.size())
    {
        FatalErrorInFunction
            << "attempt to re-construct cellModeller when it already exists :"
            << modelPtrs_.size()
            << exit(FatalError);
    }

    label maxIndex = 0;
    forAll(models_, i)
    {
        if (models_[i].index() > maxIndex) maxIndex = models_[i].index();
    }

    modelPtrs_.setSize(maxIndex + 1);
    modelPtrs_ = nullptr;

    // For all the words in the wordlist, set the details of the model
    // to those specified by the word name and the other parameters
    // given. This should result in an automatic 'read' of the model
    // from its File (see cellModel class).
    forAll(models_, i)
    {
        if (modelPtrs_[models_[i].index()])
        {
            FatalErrorInFunction
                << "more than one model share the index "
                << models_[i].index()
                << exit(FatalError);
        }

        modelPtrs_[models_[i].index()] = &models_[i];

        if (modelDictionary_.found(models_[i].name()))
        {
            FatalErrorInFunction
                << "more than one model share the name "
                << models_[i].name()
                << exit(FatalError);
        }

        modelDictionary_.insert(models_[i].name(), &models_[i]);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::cellModeller::~cellModeller()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const Foam::cellModel* Foam::cellModeller::lookup(const word& name)
{
    if (models_.empty())
    {
        IFstream is(findEtcFile("cellModels", true));
        is  >> models_;
        cellModeller m;
    }

    HashTable<const cellModel*>::iterator iter = modelDictionary_.find(name);

    if (iter != modelDictionary_.end())
    {
        return iter();
    }
    else
    {
        return nullptr;
    }
}


const Foam::cellModel* Foam::cellModeller::lookup(const label i)
{
    if (models_.empty())
    {
        IFstream is(findEtcFile("cellModels", true));
        is  >> models_;
        cellModeller m;
    }

    return modelPtrs_[i];
}


// ************************************************************************* //
