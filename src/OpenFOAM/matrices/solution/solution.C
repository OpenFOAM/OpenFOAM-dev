/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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

// These are for old syntax compatibility:
#include "BICCG.H"
#include "ICCG.H"
#include "IStringStream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineDebugSwitchWithName(solution, "solution", 0);
}

// List of sub-dictionaries to rewrite
//! \cond localScope
static const Foam::List<Foam::word> subDictNames
(
    Foam::IStringStream("(preconditioner smoother)")()
);
//! \endcond


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void Foam::solution::read(const dictionary& dict)
{
    if (dict.found("cache"))
    {
        cache_ = dict.subDict("cache");
        caching_ = cache_.lookupOrDefault("active", true);
    }

    if (dict.found("relaxationFactors"))
    {
        const dictionary& relaxDict(dict.subDict("relaxationFactors"));
        if (relaxDict.found("fields") || relaxDict.found("equations"))
        {
            if (relaxDict.found("fields"))
            {
                fieldRelaxDict_ = relaxDict.subDict("fields");
            }

            if (relaxDict.found("equations"))
            {
                eqnRelaxDict_ = relaxDict.subDict("equations");
            }
        }
        else
        {
            // backwards compatibility
            fieldRelaxDict_.clear();

            const wordList entryNames(relaxDict.toc());
            forAll(entryNames, i)
            {
                const word& e = entryNames[i];
                scalar value = readScalar(relaxDict.lookup(e));

                if (e(0, 1) == "p")
                {
                    fieldRelaxDict_.add(e, value);
                }
                else if (e.length() >= 3)
                {
                    if (e(0, 3) == "rho")
                    {
                        fieldRelaxDict_.add(e, value);
                    }
                }

            }

            eqnRelaxDict_ = relaxDict;
        }

        fieldRelaxDefault_ =
            fieldRelaxDict_.lookupOrDefault<scalar>("default", 0.0);

        eqnRelaxDefault_ =
            eqnRelaxDict_.lookupOrDefault<scalar>("default", 0.0);

        if (debug)
        {
            Info<< "relaxation factors:" << nl
                << "fields: " << fieldRelaxDict_ << nl
                << "equations: " << eqnRelaxDict_ << endl;
        }
    }


    if (dict.found("solvers"))
    {
        solvers_ = dict.subDict("solvers");
        upgradeSolverDict(solvers_);
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
            (
                obr.readOpt() == IOobject::MUST_READ
             || obr.readOpt() == IOobject::READ_IF_PRESENT
              ? IOobject::MUST_READ_IF_MODIFIED
              : obr.readOpt()
            ),
            IOobject::NO_WRITE
        )
    ),
    cache_(dictionary::null),
    caching_(false),
    fieldRelaxDict_(dictionary::null),
    eqnRelaxDict_(dictionary::null),
    fieldRelaxDefault_(0),
    eqnRelaxDefault_(0),
    solvers_(dictionary::null)
{
    if
    (
        readOpt() == IOobject::MUST_READ
     || readOpt() == IOobject::MUST_READ_IF_MODIFIED
     || (readOpt() == IOobject::READ_IF_PRESENT && headerOk())
    )
    {
        read(solutionDict());
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::solution::upgradeSolverDict
(
    dictionary& dict,
    const bool verbose
)
{
    label nChanged = 0;

    // backward compatibility:
    // recast primitive entries into dictionary entries
    forAllIter(dictionary, dict, iter)
    {
        if (!iter().isDict())
        {
            Istream& is = iter().stream();
            word name(is);
            dictionary subdict;

            if (name == "BICCG")
            {
                // special treatment for very old syntax
                subdict = BICCG::solverDict(is);
            }
            else if (name == "ICCG")
            {
                // special treatment for very old syntax
                subdict = ICCG::solverDict(is);
            }
            else
            {
                subdict.add("solver", name);
                subdict <<= dictionary(is);

                // preconditioner and smoother entries can be
                // 1) primitiveEntry w/o settings,
                // 2) or a dictionaryEntry.
                // transform primitiveEntry with settings -> dictionaryEntry
                forAll(subDictNames, dictI)
                {
                    const word& dictName = subDictNames[dictI];
                    entry* ePtr = subdict.lookupEntryPtr(dictName,false,false);

                    if (ePtr && !ePtr->isDict())
                    {
                        Istream& is = ePtr->stream();
                        is >> name;

                        if (!is.eof())
                        {
                            dictionary newDict;
                            newDict.add(dictName, name);
                            newDict <<= dictionary(is);

                            subdict.set(dictName, newDict);
                        }
                    }
                }
            }


            // write out information to help people adjust to the new syntax
            if (verbose && Pstream::master())
            {
                Info<< "// using new solver syntax:\n"
                    << iter().keyword() << subdict << endl;
            }

            // overwrite with dictionary entry
            dict.set(iter().keyword(), subdict);

            nChanged++;
        }
    }

    return nChanged;
}


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


bool Foam::solution::relaxField(const word& name) const
{
    if (debug)
    {
        Info<< "Field relaxation factor for " << name
            << " is " << (fieldRelaxDict_.found(name) ? "set" : "unset")
            << endl;
    }

    return fieldRelaxDict_.found(name) || fieldRelaxDict_.found("default");
}


bool Foam::solution::relaxEquation(const word& name) const
{
    if (debug)
    {
        Info<< "Find equation relaxation factor for " << name << endl;
    }

    return eqnRelaxDict_.found(name) || eqnRelaxDict_.found("default");
}


Foam::scalar Foam::solution::fieldRelaxationFactor(const word& name) const
{
    if (debug)
    {
        Info<< "Lookup variable relaxation factor for " << name << endl;
    }

    if (fieldRelaxDict_.found(name))
    {
        return readScalar(fieldRelaxDict_.lookup(name));
    }
    else if (fieldRelaxDefault_ > SMALL)
    {
        return fieldRelaxDefault_;
    }
    else
    {
        FatalIOErrorIn
        (
            "Foam::solution::fieldRelaxationFactor(const word&)",
            fieldRelaxDict_
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

    if (eqnRelaxDict_.found(name))
    {
        return readScalar(eqnRelaxDict_.lookup(name));
    }
    else if (eqnRelaxDefault_ > SMALL)
    {
        return eqnRelaxDefault_;
    }
    else
    {
        FatalIOErrorIn
        (
            "Foam::solution::eqnRelaxationFactor(const word&)",
            eqnRelaxDict_
        )   << "Cannot find equation relaxation factor for '" << name
            << "' or a suitable default value."
            << exit(FatalIOError);

        return 0;
    }
}


const Foam::dictionary& Foam::solution::solutionDict() const
{
    if (found("select"))
    {
        return subDict(word(lookup("select")));
    }
    else
    {
        return *this;
    }
}


const Foam::dictionary& Foam::solution::solverDict(const word& name) const
{
    if (debug)
    {
        InfoIn("solution::solverDict(const word&)")
            << "Lookup solver for " << name << endl;
    }

    return solvers_.subDict(name);
}


const Foam::dictionary& Foam::solution::solver(const word& name) const
{
    if (debug)
    {
        InfoIn("solution::solver(const word&)")
            << "Lookup solver for " << name << endl;
    }

    return solvers_.subDict(name);
}


bool Foam::solution::read()
{
    if (regIOobject::read())
    {
        read(solutionDict());

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
