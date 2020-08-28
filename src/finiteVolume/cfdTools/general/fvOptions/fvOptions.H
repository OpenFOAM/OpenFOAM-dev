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
    Foam::fv::options

Description
    Finite-volume options

SourceFiles
    options.C

\*---------------------------------------------------------------------------*/

#ifndef fvOptions_H
#define fvOptions_H

#include "fvOptionList.H"
#include "MeshObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                           Class options Declaration
\*---------------------------------------------------------------------------*/

class options
:
    public MeshObject<fvMesh, UpdateableMeshObject, options>,
    public dictionary,
    public optionList
{
    // Private Member Functions

        //- Create IO object if dictionary is present
        IOobject createIOobject(const fvMesh& mesh) const;


public:

    //- Runtime type information
    TypeName("fvOptions");


    // Constructors

        //- Construct from components with list of field names
        options(const fvMesh& mesh);

        //- Disallow default bitwise copy construction
        options(const options&) = delete;

        using MeshObject<fvMesh, UpdateableMeshObject, options>::New;


    //- Destructor
    virtual ~options()
    {}


    // Member Functions

        //- Declare fvOptions to be a global dictionary
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

        //- ReadData function required for regIOobject read operation
        virtual bool readData(Istream&);

        using optionList::read;

        //- Read dictionary
        virtual bool read();

        //- WriteData function required for regIOobject write operation
        virtual bool writeData(Ostream&) const;


        // Mesh changes

            //- Update for mesh changes
            virtual void updateMesh(const mapPolyMesh&);

            //- Update for mesh motion
            virtual bool movePoints();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const options&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Template function for obtaining global status
template<>
inline bool typeGlobal<fv::options>()
{
    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
