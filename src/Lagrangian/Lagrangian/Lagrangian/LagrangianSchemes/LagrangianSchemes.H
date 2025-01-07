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
    Foam::LagrangianSchemes

Description
    Selector class for Lagrangian schemes

SourceFiles
    LagrangianSchemes.C

\*---------------------------------------------------------------------------*/

#ifndef LagrangianSchemes_H
#define LagrangianSchemes_H

#include "IOdictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class LagrangianSchemes Declaration
\*---------------------------------------------------------------------------*/

class LagrangianSchemes
:
    public IOdictionary
{
    // Private Data

        //- Time schemes sub-dictionary
        dictionary ddtSchemes_;

        //- Default time scheme stream
        ITstream defaultDdtScheme_;

        //- Source schemes sub-dictionary
        dictionary SpSchemes_;

        //- Default source scheme stream
        ITstream defaultSpScheme_;

        //- Averaging schemes sub-dictionary
        dictionary averagingSchemes_;

        //- Default averaging scheme stream
        ITstream defaultAveragingScheme_;

        //- Interpolation schemes sub-dictionary
        dictionary interpolationSchemes_;

        //- Default interpolation scheme stream
        ITstream defaultInterpolationScheme_;

        //- Accumulation schemes sub-dictionary
        dictionary accumulationSchemes_;

        //- Default accumulation scheme stream
        ITstream defaultAccumulationScheme_;


    // Private Member Functions

        //- Clear the dictionaries and streams before reading
        void clear();

        //- Read a scheme from the dictionary
        static void read
        (
            const dictionary& dict,
            const word& type,
            dictionary& typeSchemes,
            ITstream& defaultTypeScheme
        );

        //- Read schemes from the dictionary
        void read(const dictionary&);

        //- Look up a scheme
        static ITstream& lookup
        (
            const word& name,
            const dictionary& typeSchemes,
            const ITstream& defaultTypeScheme
        );


public:

    //- Type information
    ClassName("LagrangianSchemes");


    // Constructors

        //- Construct for objectRegistry
        LagrangianSchemes(const objectRegistry& db);

        //- Disallow default bitwise copy construction
        LagrangianSchemes(const LagrangianSchemes&) = delete;


    //- Destructor
    ~LagrangianSchemes();


    // Member Functions

        // Access

            //- Access the dictionary
            const dictionary& schemesDict() const;

            //- Get the time scheme for the given field name
            ITstream& ddt(const word& name) const;

            //- Get the source scheme for the given field name
            ITstream& Sp(const word& name) const;

            //- Get the averaging scheme for the given field name
            ITstream& averaging(const word& name) const;

            //- Get the interpolation scheme for the given field name
            ITstream& interpolation(const word& name) const;

            //- Get the accumulation scheme for the given field name
            ITstream& accumulation(const word& name) const;


        // Read

            //- Read the LagrangianSchemes
            bool read();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const LagrangianSchemes&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
