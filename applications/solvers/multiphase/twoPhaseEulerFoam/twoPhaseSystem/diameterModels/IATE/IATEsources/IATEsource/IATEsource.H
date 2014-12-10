/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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
    Foam::diameterModels::IATEsource

Description
    IATE (Interfacial Area Transport Equation) bubble diameter model
    run-time selectable sources.

SourceFiles
    IATEsource.C

\*---------------------------------------------------------------------------*/

#ifndef IATEsource_H
#define IATEsource_H

#include "IATE.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{

/*---------------------------------------------------------------------------*\
                           Class IATEsource Declaration
\*---------------------------------------------------------------------------*/

class IATEsource
{

protected:

    // Protected data

        //- Reference to the IATE this source applies to
        const IATE& iate_;


public:

    //- Runtime type information
    TypeName("IATEsource");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            IATEsource,
            dictionary,
            (
                const IATE& iate,
                const dictionary& dict
            ),
            (iate, dict)
        );


    //- Class used for the read-construction of
    //  PtrLists of IATE sources
    class iNew
    {
        const IATE& iate_;

    public:

        iNew(const IATE& iate)
        :
            iate_(iate)
        {}

        autoPtr<IATEsource> operator()(Istream& is) const
        {
            word type(is);
            dictionary dict(is);
            return IATEsource::New(type, iate_, dict);
        }
    };


    // Constructors

        IATEsource(const IATE& iate)
        :
            iate_(iate)
        {}

        autoPtr<IATEsource> clone() const
        {
            notImplemented("autoPtr<IATEsource> clone() const");
            return autoPtr<IATEsource>(NULL);
        }


    // Selectors

        static autoPtr<IATEsource> New
        (
            const word& type,
            const IATE& iate,
            const dictionary& dict
        );


    //- Destructor
    virtual ~IATEsource()
    {}


    // Member Functions

        const phaseModel& phase() const
        {
            return iate_.phase();
        }

        const twoPhaseSystem& fluid() const
        {
            return iate_.phase().fluid();
        }

        const phaseModel& otherPhase() const
        {
            return phase().otherPhase();
        }

        scalar phi() const
        {
            return 1.0/(36*constant::mathematical::pi);
        }

        //- Return the bubble relative velocity
        tmp<volScalarField> Ur() const;

        //- Return the bubble turbulent velocity
        tmp<volScalarField> Ut() const;

        //- Return the bubble Reynolds number
        tmp<volScalarField> Re() const;

        //- Return the bubble drag coefficient
        tmp<volScalarField> CD() const;

        //- Return the bubble Morton number
        tmp<volScalarField> Mo() const;

        //- Return the bubble Eotvos number
        tmp<volScalarField> Eo() const;

        //- Return the bubble Webber number
        tmp<volScalarField> We() const;

        virtual tmp<volScalarField> R() const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace diameterModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
