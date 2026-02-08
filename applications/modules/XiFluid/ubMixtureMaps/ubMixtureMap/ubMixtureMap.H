/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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
    Foam::ubMixtureMap

Description
    Base class for unburnt/burnt gas composition mapping

SourceFiles
    ubMixtureMap.C

\*---------------------------------------------------------------------------*/

#ifndef ubMixtureMap_H
#define ubMixtureMap_H

#include "uRhoMulticomponentThermo.H"
#include "bRhoMulticomponentThermo.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class ubMixtureMap Declaration
\*---------------------------------------------------------------------------*/

class ubMixtureMap
{
protected:

    // Protected data

        const uRhoMulticomponentThermo& uThermo_;
        const bRhoMulticomponentThermo& bThermo_;


public:

    //- Runtime type information
    TypeName("ubMixtureMap");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            ubMixtureMap,
            thermo,
            (
                const uRhoMulticomponentThermo& uThermo,
                const bRhoMulticomponentThermo& bThermo
            ),
            (
                uThermo,
                bThermo
            )
        );


    // Constructors

        //- Construct from components
        ubMixtureMap
        (
            const uRhoMulticomponentThermo& uThermo,
            const bRhoMulticomponentThermo& bThermo
        );

        //- Disallow default bitwise copy construction
        ubMixtureMap(const ubMixtureMap&) = delete;


    // Selectors

        //- Return a reference to the selected XiEq model
        static autoPtr<ubMixtureMap> New
        (
            const uRhoMulticomponentThermo& uThermo,
            const bRhoMulticomponentThermo& bThermo
        );


    //- Destructor
    virtual ~ubMixtureMap();


    // Member Functions

        template<class To>
        inline const To& uMixtureCast() const
        {
            try
            {
                return dynamic_cast<const To&>(uThermo_);
            }
            catch (const std::bad_cast&)
            {
                FatalErrorInFunction
                    << "Attempt to cast type " << uThermo_.type()
                    << " to type " << To::typeName
                    << abort(FatalError);

                return dynamic_cast<const To&>(uThermo_);
            }
        }

        template<class To>
        inline const To& bMixtureCast() const
        {
            try
            {
                return dynamic_cast<const To&>(bThermo_);
            }
            catch (const std::bad_cast&)
            {
                FatalErrorInFunction
                    << "Attempt to cast type " << bThermo_.type()
                    << " to type " << To::typeName
                    << abort(FatalError);

                return dynamic_cast<const To&>(bThermo_);
            }
        }

        virtual PtrList<volScalarField::Internal> prompt
        (
            const PtrList<volScalarField>& Yu
        ) const = 0;

        virtual void reset
        (
            const volScalarField& b,
            UPtrList<volScalarField>& Yu,
            const volScalarField& c,
            const UPtrList<const volScalarField>& Yb
        ) const = 0;

        //- Update properties from the given XiProperties dictionary
        bool read(const dictionary& XiProperties);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const ubMixtureMap&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
