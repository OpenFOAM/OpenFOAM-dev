/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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

#include "ejectionModelList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace surfaceFilmModels
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

ejectionModelList::ejectionModelList(surfaceFilm& film)
:
    PtrList<ejectionModel>(),
    filmSubModelBase(film)
{}


ejectionModelList::ejectionModelList
(
    surfaceFilm& film,
    const dictionary& dict
)
:
    PtrList<ejectionModel>(),
    filmSubModelBase
    (
        "ejectionModelList",
        film,
        dict,
        "ejectionModelList",
        "ejectionModelList"
    ),
    massEjected_(film.intCoupledPatchIDs().size(), 0.0)
{
    Info<< "    Selecting film ejection" << endl;

    if (dict.isDict("ejection"))
    {
        const dictionary& ejectionDict(dict.subDict("ejection"));
        this->setSize(ejectionDict.size());

        label i = 0;
        forAllConstIter(dictionary, ejectionDict, iter)
        {
            set
            (
                i++,
                ejectionModel::New
                (
                    film,
                    ejectionDict.isDict(iter().keyword())
                  ? ejectionDict.subDict(iter().keyword())
                  : dictionary::null,
                    iter().keyword()
                )
            );
        }
    }
    else if (dict.found("ejectionModels"))
    {
        const wordList models(dict.lookup("ejectionModels"));
        this->setSize(models.size());

        forAll(models, i)
        {
            set(i, ejectionModel::New(film, dict, models[i]));
        }
    }

    if (!size())
    {
        Info<< "        none" << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

ejectionModelList::~ejectionModelList()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void ejectionModelList::correct
(
    scalarField& availableMass,
    volScalarField& massToEject,
    volScalarField& diameterToEject
)
{
    // Correct models that accumulate mass and diameter transfers
    forAll(*this, i)
    {
        ejectionModel& im = operator[](i);
        im.correct(availableMass, massToEject, diameterToEject);
    }

    // Push values to boundaries ready for transfer to the primary region
    massToEject.correctBoundaryConditions();
    diameterToEject.correctBoundaryConditions();

    const labelList& patchIDs = film().intCoupledPatchIDs();

    forAll(patchIDs, i)
    {
        label patchi = patchIDs[i];
        massEjected_[i] =
            massEjected_[i] + sum(massToEject.boundaryField()[patchi]);
    }
}


void ejectionModelList::info(Ostream& os)
{
    const polyBoundaryMesh& pbm = film().mesh().boundaryMesh();

    scalar ejectedMass = 0;
    scalarField patchEjectedMasses
    (
        pbm.size() - film().mesh().globalData().processorPatches().size(),
        0
    );

    forAll(*this, i)
    {
        const ejectionModel& im = operator[](i);
        ejectedMass += im.ejectedMassTotal();
        im.patchEjectedMassTotals(patchEjectedMasses);
    }

    os  << indent << "ejected mass      = " << ejectedMass << nl;

    forAll(patchEjectedMasses, patchi)
    {
        if (mag(patchEjectedMasses[patchi]) > vSmall)
        {
            os  << indent << indent << "from patch " << pbm[patchi].name()
                << " = " << patchEjectedMasses[patchi] << nl;
        }
    }

    scalarField mass0(massEjected_.size(), 0);
    this->getBaseProperty("massEjected", mass0);

    scalarField mass(massEjected_);
    Pstream::listCombineGather(mass, plusEqOp<scalar>());
    mass += mass0;

    const labelList& patchIDs = film().intCoupledPatchIDs();

    forAll(patchIDs, i)
    {
        label patchi = patchIDs[i];
        Info<< indent << "  - patch: " << pbm[patchi].name() << ": "
            << mass[i] << endl;
    }

    if (film().time().writeTime())
    {
        setBaseProperty("massEjected", mass);
        massEjected_ = 0.0;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace Foam

// ************************************************************************* //
