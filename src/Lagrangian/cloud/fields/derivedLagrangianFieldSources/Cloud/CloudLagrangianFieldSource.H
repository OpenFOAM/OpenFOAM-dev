/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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
    Foam::CloudLagrangianFieldSource

Description
    Templated base class for source conditions that refer to a cloud. Not a
    source condition in itself. Must be mixed in.

SourceFiles
    CloudLagrangianFieldSource.C

\*---------------------------------------------------------------------------*/

#ifndef CloudLagrangianFieldSource_H
#define CloudLagrangianFieldSource_H

#include "LagrangianFieldsFwd.H"
#include "LagrangianSubFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class LagrangianModel;
class cloud;

/*---------------------------------------------------------------------------*\
                 Class CloudLagrangianFieldSource Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class CloudLagrangianFieldSource
{
    // Private Member Data

        //- Reference to the derived field source
        const LagrangianFieldSource<Type>& field_;

        //- Reference to the cloud
        const Foam::cloud& cloud_;


protected:

    // Protected Member Functions

        //- Determine whether the cloud is one of the given types
        template<class Cloud, class ... Clouds>
        bool isCloud() const;

        //- Generate an error if the cloud is not one of the given types
        template<class Cloud, class ... Clouds>
        void assertCloud
        (
            const LagrangianModel& model,
            const LagrangianSubMesh& subMesh
        ) const;

        //- Get a reference to the cloud
        template<class Cloud = Foam::cloud>
        const Cloud& cloud
        (
            const LagrangianModel& model,
            const LagrangianSubMesh& subMesh
        ) const;


public:

    // Constructors

        //- Construct with reference to the derived field source
        template<class Derived>
        CloudLagrangianFieldSource(const Derived&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "CloudLagrangianFieldSource.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
