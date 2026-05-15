/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2026 OpenFOAM Foundation
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

#include "solution.H"
#include "Time.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(solution, 0);
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void Foam::solution::readDict()
{
    if (found("cache"))
    {
        printDictionary print(subDict("cache"));
        cache_ = subDict("cache");
        caching_ = cache_.lookupOrDefault("active", true);
    }

    fieldRelaxDict_ = &dictionary::null;
    eqnRelaxDict_ = &dictionary::null;
    if (found("relaxationFactors"))
    {
        const dictionary& relaxDict(subDict("relaxationFactors"));
        printDictionary print(relaxDict);

        if (relaxDict.found("fields") || relaxDict.found("equations"))
        {
            if (relaxDict.found("fields"))
            {
                fieldRelaxDict_ = &relaxDict.subDict("fields");
            }

            if (relaxDict.found("equations"))
            {
                eqnRelaxDict_ = &relaxDict.subDict("equations");
            }
        }
        else
        {
            IOWarningInFunction(*this)
                << "Neither fields nor equations specified" << endl;
        }
    }

    if (found("solvers"))
    {
        solvers_ = &subDict("solvers");
        printDictionary print(*solvers_);
    }
    else
    {
        solvers_ = &dictionary::null;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solution::solution
(
    const objectRegistry& obr,
    const fileName& dictName
)
:
    IOdictionary
    (
        IOobject
        (
            dictName,
            obr.time().system(),
            obr,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    cache_("cache", *this),
    caching_(false),
    fieldRelaxDict_(nullptr),
    eqnRelaxDict_(nullptr),
    fieldRelaxDefault_(0),
    eqnRelaxDefault_(0),
    solvers_(nullptr)
{
    readDict();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::solution::cache(const word& name) const
{
    if (caching_)
    {
        if (debug)
        {
            Info<< "Cache: find entry for " << name << endl;
        }

        return cache_.found(name);
    }
    else
    {
        return false;
    }
}


void Foam::solution::enableCache(const word& name) const
{
    caching_ = true;

    if (debug)
    {
        Info<< "Enable cache for " << name << endl;
    }

    cache_.add(name, true);
}


bool Foam::solution::relaxField(const word& name) const
{
    if (debug)
    {
        Info<< "Field relaxation factor for " << name
            << " is " << (fieldRelaxDict_->found(name) ? "set" : "unset")
            << endl;
    }

    return fieldRelaxDict_->found(name) || fieldRelaxDict_->found("default");
}


bool Foam::solution::relaxEquation(const word& name) const
{
    if (debug)
    {
        Info<< "Find equation relaxation factor for " << name << endl;
    }

    return eqnRelaxDict_->found(name) || eqnRelaxDict_->found("default");
}


Foam::scalar Foam::solution::fieldRelaxationFactor(const word& name) const
{
    if (debug)
    {
        Info<< "Lookup variable relaxation factor for " << name << endl;
    }

    if (fieldRelaxDict_->found(name))
    {
        return fieldRelaxDict_->lookup<scalar>(name);
    }
    else if (fieldRelaxDefault_ > small)
    {
        return fieldRelaxDefault_;
    }
    else
    {
        FatalIOErrorInFunction
        (
            *fieldRelaxDict_
        )   << "Cannot find variable relaxation factor for '" << name
            << "' or a suitable default value."
            << exit(FatalIOError);

        return 0;
    }
}


Foam::scalar Foam::solution::equationRelaxationFactor(const word& name) const
{
    if (debug)
    {
        Info<< "Lookup equation relaxation factor for " << name << endl;
    }

    if (eqnRelaxDict_->found(name))
    {
        return eqnRelaxDict_->lookup<scalar>(name);
    }
    else if (eqnRelaxDefault_ > small)
    {
        return eqnRelaxDefault_;
    }
    else
    {
        FatalIOErrorInFunction
        (
            *eqnRelaxDict_
        )   << "Cannot find equation relaxation factor for '" << name
            << "' or a suitable default value."
            << exit(FatalIOError);

        return 0;
    }
}


const Foam::dictionary& Foam::solution::solversDict() const
{
    return *solvers_;
}


const Foam::dictionary& Foam::solution::solverDict(const word& name) const
{
    if (debug)
    {
        Info<< "Lookup solver for " << name << endl;
    }

    return solvers_->subDict(name);
}


bool Foam::solution::read()
{
    if (regIOobject::read())
    {
        readDict();

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
