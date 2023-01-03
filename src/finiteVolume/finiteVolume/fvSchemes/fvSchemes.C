/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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

#include "fvSchemes.H"
#include "Time.H"
#include "steadyStateDdtScheme.H"
#include "steadyStateD2dt2Scheme.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(fvSchemes, 0);
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void Foam::fvSchemes::clear()
{
    ddtSchemes_.clear();
    defaultDdtScheme_.clear();
    d2dt2Schemes_.clear();
    defaultD2dt2Scheme_.clear();
    interpolationSchemes_.clear();
    defaultInterpolationScheme_.clear();
    divSchemes_.clear();
    defaultDivScheme_.clear();
    gradSchemes_.clear();
    defaultGradScheme_.clear();
    snGradSchemes_.clear();
    defaultSnGradScheme_.clear();
    laplacianSchemes_.clear();
    defaultLaplacianScheme_.clear();
    // Do not clear fluxRequired settings
}


void Foam::fvSchemes::read(const dictionary& dict)
{
    if (dict.found("ddtSchemes"))
    {
        ddtSchemes_ = dict.subDict("ddtSchemes");
    }
    else
    {
        ddtSchemes_.set("default", "none");
    }

    if
    (
        ddtSchemes_.found("default")
     && word(ddtSchemes_.lookup("default")) != "none"
    )
    {
        defaultDdtScheme_ = ddtSchemes_.lookup("default");
        steady_ =
        (
            word(defaultDdtScheme_)
         == fv::steadyStateDdtScheme<scalar>::typeName
        );
    }


    if (dict.found("d2dt2Schemes"))
    {
        d2dt2Schemes_ = dict.subDict("d2dt2Schemes");
    }
    else
    {
        d2dt2Schemes_.set("default", "none");
    }

    if
    (
        d2dt2Schemes_.found("default")
     && word(d2dt2Schemes_.lookup("default")) != "none"
    )
    {
        defaultD2dt2Scheme_ = d2dt2Schemes_.lookup("default");
        steady_ =
        (
            word(defaultD2dt2Scheme_)
         == fv::steadyStateD2dt2Scheme<scalar>::typeName
        );
    }


    if (dict.found("interpolationSchemes"))
    {
        interpolationSchemes_ = dict.subDict("interpolationSchemes");
    }
    else if (!interpolationSchemes_.found("default"))
    {
        interpolationSchemes_.add("default", "linear");
    }

    if
    (
        interpolationSchemes_.found("default")
     && word(interpolationSchemes_.lookup("default")) != "none"
    )
    {
        defaultInterpolationScheme_ =
            interpolationSchemes_.lookup("default");
    }


    divSchemes_ = dict.subDict("divSchemes");

    if
    (
        divSchemes_.found("default")
     && word(divSchemes_.lookup("default")) != "none"
    )
    {
        defaultDivScheme_ = divSchemes_.lookup("default");
    }


    gradSchemes_ = dict.subDict("gradSchemes");

    if
    (
        gradSchemes_.found("default")
     && word(gradSchemes_.lookup("default")) != "none"
    )
    {
        defaultGradScheme_ = gradSchemes_.lookup("default");
    }


    if (dict.found("snGradSchemes"))
    {
        snGradSchemes_ = dict.subDict("snGradSchemes");
    }
    else if (!snGradSchemes_.found("default"))
    {
        snGradSchemes_.add("default", "corrected");
    }

    if
    (
        snGradSchemes_.found("default")
     && word(snGradSchemes_.lookup("default")) != "none"
    )
    {
        defaultSnGradScheme_ = snGradSchemes_.lookup("default");
    }


    laplacianSchemes_ = dict.subDict("laplacianSchemes");

    if
    (
        laplacianSchemes_.found("default")
     && word(laplacianSchemes_.lookup("default")) != "none"
    )
    {
        defaultLaplacianScheme_ = laplacianSchemes_.lookup("default");
    }


    if (dict.found("fluxRequired"))
    {
        fluxRequired_.merge(dict.subDict("fluxRequired"));

        if
        (
            fluxRequired_.found("default")
         && word(fluxRequired_.lookup("default")) != "none"
        )
        {
            defaultFluxRequired_ = Switch(fluxRequired_.lookup("default"));
        }
    }
}


Foam::word Foam::fvSchemes::filterGroup(const word& name) const
{
    word filteredName(name);

    word::size_type n = 0;
    word::iterator iter2 = filteredName.begin();

    bool found = false;

    for
    (
        word::const_iterator iter1 = filteredName.begin();
        iter1 != filteredName.end();
        ++iter1
    )
    {
        const char c = *iter1;

        if (c == '.')
        {
            found = true;
        }
        else if (!found || !isalnum(c))
        {
            found = false;
            *iter2++ = c;
            n++;
        }
    }

    filteredName.resize(n);

    return filteredName;
}


Foam::ITstream& Foam::fvSchemes::lookupScheme
(
    const word& name,
    const dictionary& schemes,
    const ITstream& defaultScheme
) const
{
    if (debug)
    {
        Info<< "Lookup scheme for " << name
            << " in dictionary " << schemes.name().caseName() << endl;
    }

    // Lookup scheme with optional wildcards
    const entry* schemePtr = schemes.lookupEntryPtr(name, false, true);

    if (schemePtr)
    {
        return schemePtr->stream();
    }
    else
    {
        // If scheme not found check if it contains group names
        bool filtered = (name.find('.') != word::npos);

        if (filtered)
        {
            // Filter-out the group names and lookup
            const entry* schemePtr =
                schemes.lookupEntryPtr(filterGroup(name), false, true);

            if (schemePtr)
            {
                return schemePtr->stream();
            }
        }

        // If scheme still not found check if a default scheme is provided
        if (!defaultScheme.empty())
        {
            const_cast<ITstream&>(defaultScheme).rewind();
            return const_cast<ITstream&>(defaultScheme);
        }
        else
        {
            // Cannot find scheme

            FatalIOErrorInFunction(schemes)
                << "Cannot find scheme for " << name;

            if (filtered)
            {
                FatalIOError << " or " << filterGroup(name);
            }

            FatalIOError
                << " in dictionary " << schemes.name().caseName()
                << exit(FatalIOError);
            return const_cast<ITstream&>(defaultScheme);
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvSchemes::fvSchemes(const objectRegistry& obr)
:
    IOdictionary
    (
        IOobject
        (
            "fvSchemes",
            obr.time().system(),
            obr,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    ddtSchemes_
    (
        ITstream
        (
            objectPath() + ".ddtSchemes",
            tokenList()
        )()
    ),
    defaultDdtScheme_
    (
        ddtSchemes_.name() + ".default",
        tokenList()
    ),
    d2dt2Schemes_
    (
        ITstream
        (
            objectPath() + ".d2dt2Schemes",
            tokenList()
        )()
    ),
    defaultD2dt2Scheme_
    (
        d2dt2Schemes_.name() + ".default",
        tokenList()
    ),
    interpolationSchemes_
    (
        ITstream
        (
            objectPath() + ".interpolationSchemes",
            tokenList()
        )()
    ),
    defaultInterpolationScheme_
    (
        interpolationSchemes_.name() + ".default",
        tokenList()
    ),
    divSchemes_
    (
        ITstream
        (
            objectPath() + ".divSchemes",
            tokenList()
        )()
    ),
    defaultDivScheme_
    (
        divSchemes_.name() + ".default",
        tokenList()
    ),
    gradSchemes_
    (
        ITstream
        (
            objectPath() + ".gradSchemes",
            tokenList()
        )()
    ),
    defaultGradScheme_
    (
        gradSchemes_.name() + ".default",
        tokenList()
    ),
    snGradSchemes_
    (
        ITstream
        (
            objectPath() + ".snGradSchemes",
            tokenList()
        )()
    ),
    defaultSnGradScheme_
    (
        snGradSchemes_.name() + ".default",
        tokenList()
    ),
    laplacianSchemes_
    (
        ITstream
        (
            objectPath() + ".laplacianSchemes",
            tokenList()
        )()
    ),
    defaultLaplacianScheme_
    (
        laplacianSchemes_.name() + ".default",
        tokenList()
    ),
    fluxRequired_
    (
        ITstream
        (
            objectPath() + ".fluxRequired",
            tokenList()
        )()
    ),
    defaultFluxRequired_(false),
    steady_(false)
{
    read(dict());
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fvSchemes::read()
{
    if (regIOobject::read())
    {
        // Clear current settings except fluxRequired
        clear();

        read(dict());

        return true;
    }
    else
    {
        return false;
    }
}


const Foam::dictionary& Foam::fvSchemes::dict() const
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


Foam::ITstream& Foam::fvSchemes::ddt(const word& name) const
{
    return lookupScheme(name, ddtSchemes_, defaultDdtScheme_);
}


Foam::ITstream& Foam::fvSchemes::d2dt2(const word& name) const
{
    return lookupScheme(name, d2dt2Schemes_, defaultD2dt2Scheme_);
}


Foam::ITstream& Foam::fvSchemes::interpolation(const word& name) const
{
    return lookupScheme
    (
        name,
        interpolationSchemes_,
        defaultInterpolationScheme_
    );
}


Foam::ITstream& Foam::fvSchemes::div(const word& name) const
{
    return lookupScheme(name, divSchemes_, defaultDivScheme_);
}


Foam::ITstream& Foam::fvSchemes::grad(const word& name) const
{
    return lookupScheme(name, gradSchemes_, defaultGradScheme_);
}


Foam::ITstream& Foam::fvSchemes::snGrad(const word& name) const
{
    return lookupScheme(name, snGradSchemes_, defaultSnGradScheme_);
}


Foam::ITstream& Foam::fvSchemes::laplacian(const word& name) const
{
    return lookupScheme(name, laplacianSchemes_, defaultLaplacianScheme_);
}


void Foam::fvSchemes::setFluxRequired(const word& name) const
{
    if (debug)
    {
        Info<< "Setting fluxRequired for " << name << endl;
    }

    fluxRequired_.add(name, true, true);
}


bool Foam::fvSchemes::fluxRequired(const word& name) const
{
    if (debug)
    {
        Info<< "Lookup fluxRequired for " << name << endl;
    }

    if (fluxRequired_.found(name))
    {
        return true;
    }
    else
    {
        return defaultFluxRequired_;
    }
}


// ************************************************************************* //
