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

\*---------------------------------------------------------------------------*/

#include "fvSchemes.H"
#include "Time.H"
#include "steadyStateDdtScheme.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

int Foam::fvSchemes::debug(Foam::debug::debugSwitch("fvSchemes", false));


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
            (
                obr.readOpt() == IOobject::MUST_READ
             || obr.readOpt() == IOobject::READ_IF_PRESENT
              ? IOobject::MUST_READ_IF_MODIFIED
              : obr.readOpt()
            ),
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
    if
    (
        readOpt() == IOobject::MUST_READ
     || readOpt() == IOobject::MUST_READ_IF_MODIFIED
     || (readOpt() == IOobject::READ_IF_PRESENT && headerOk())
    )
    {
        read(schemesDict());
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fvSchemes::read()
{
    if (regIOobject::read())
    {
        // Clear current settings except fluxRequired
        clear();

        read(schemesDict());

        return true;
    }
    else
    {
        return false;
    }
}


const Foam::dictionary& Foam::fvSchemes::schemesDict() const
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


Foam::ITstream& Foam::fvSchemes::ddtScheme(const word& name) const
{
    if (debug)
    {
        Info<< "Lookup ddtScheme for " << name << endl;
    }

    if (ddtSchemes_.found(name) || defaultDdtScheme_.empty())
    {
        return ddtSchemes_.lookup(name);
    }
    else
    {
        const_cast<ITstream&>(defaultDdtScheme_).rewind();
        return const_cast<ITstream&>(defaultDdtScheme_);
    }
}


Foam::ITstream& Foam::fvSchemes::d2dt2Scheme(const word& name) const
{
    if (debug)
    {
        Info<< "Lookup d2dt2Scheme for " << name << endl;
    }

    if (d2dt2Schemes_.found(name) || defaultD2dt2Scheme_.empty())
    {
        return d2dt2Schemes_.lookup(name);
    }
    else
    {
        const_cast<ITstream&>(defaultD2dt2Scheme_).rewind();
        return const_cast<ITstream&>(defaultD2dt2Scheme_);
    }
}


Foam::ITstream& Foam::fvSchemes::interpolationScheme(const word& name) const
{
    if (debug)
    {
        Info<< "Lookup interpolationScheme for " << name << endl;
    }

    if
    (
        interpolationSchemes_.found(name)
     || defaultInterpolationScheme_.empty()
    )
    {
        return interpolationSchemes_.lookup(name);
    }
    else
    {
        const_cast<ITstream&>(defaultInterpolationScheme_).rewind();
        return const_cast<ITstream&>(defaultInterpolationScheme_);
    }
}


Foam::ITstream& Foam::fvSchemes::divScheme(const word& name) const
{
    if (debug)
    {
        Info<< "Lookup divScheme for " << name << endl;
    }

    if (divSchemes_.found(name) || defaultDivScheme_.empty())
    {
        return divSchemes_.lookup(name);
    }
    else
    {
        const_cast<ITstream&>(defaultDivScheme_).rewind();
        return const_cast<ITstream&>(defaultDivScheme_);
    }
}


Foam::ITstream& Foam::fvSchemes::gradScheme(const word& name) const
{
    if (debug)
    {
        Info<< "Lookup gradScheme for " << name << endl;
    }

    if (gradSchemes_.found(name) || defaultGradScheme_.empty())
    {
        return gradSchemes_.lookup(name);
    }
    else
    {
        const_cast<ITstream&>(defaultGradScheme_).rewind();
        return const_cast<ITstream&>(defaultGradScheme_);
    }
}


Foam::ITstream& Foam::fvSchemes::snGradScheme(const word& name) const
{
    if (debug)
    {
        Info<< "Lookup snGradScheme for " << name << endl;
    }

    if (snGradSchemes_.found(name) || defaultSnGradScheme_.empty())
    {
        return snGradSchemes_.lookup(name);
    }
    else
    {
        const_cast<ITstream&>(defaultSnGradScheme_).rewind();
        return const_cast<ITstream&>(defaultSnGradScheme_);
    }
}


Foam::ITstream& Foam::fvSchemes::laplacianScheme(const word& name) const
{
    if (debug)
    {
        Info<< "Lookup laplacianScheme for " << name << endl;
    }

    if (laplacianSchemes_.found(name) || defaultLaplacianScheme_.empty())
    {
        return laplacianSchemes_.lookup(name);
    }
    else
    {
        const_cast<ITstream&>(defaultLaplacianScheme_).rewind();
        return const_cast<ITstream&>(defaultLaplacianScheme_);
    }
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
