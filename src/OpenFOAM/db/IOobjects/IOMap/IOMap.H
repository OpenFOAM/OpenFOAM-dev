/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

Class
    Foam::IOMap

Description
    A Map of objects of type \<T\> with automated input and output.
    Is a global object; i.e. can be read from undecomposed case.

SourceFiles
    IOMap.C

\*---------------------------------------------------------------------------*/

#ifndef IOMap_H
#define IOMap_H

#include "Map.H"
#include "regIOobject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                            Class IOMap Declaration
\*---------------------------------------------------------------------------*/

template<class T>
class IOMap
:
    public regIOobject,
    public Map<T>
{

public:

    //- Runtime type information
    TypeName("Map");


    // Constructors

        //- Construct from IOobject
        IOMap(const IOobject&);

        //- Construct from IOobject and size of Map
        IOMap(const IOobject&, const label);

        //- Construct from IOobject and a Map
        IOMap(const IOobject&, const Map<T>&);

        //- Move construct by transferring the Map contents
        IOMap(const IOobject&, Map<T>&&);

        //- Move constructor
        IOMap(IOMap<T>&&);


    //- Destructor
    virtual ~IOMap();


    // Member Functions

        bool writeData(Ostream&) const;

        //- Is object global
        virtual bool global() const
        {
            return true;
        }

        //- Return complete path + object name if the file exists
        //  either in the case/processor or case otherwise null
        virtual fileName filePath() const
        {
            return globalFilePath(type());
        }


    // Member Operators

        //- Assignment of other IOMap's entries to this IOMap
        void operator=(const IOMap<T>&);

        //- Move assignment of other IOMap's entries to this IOMap
        void operator=(IOMap<T>&&);

        //- Assignment of other Map's entries to this IOMap
        void operator=(const Map<T>&);

        //- Move assignment of other Map's entries to this IOMap
        void operator=(Map<T>&&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "IOMap.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
