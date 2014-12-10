/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2014 OpenFOAM Foundation
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
    Foam::fv::MRFSource

Description
    Multiple Reference Frame (MRF) source

    \heading Source usage

    Example usage:
    \verbatim
    MRFSourceCoeffs
    {
        origin      (0 0 0);
        axis        (0 0 1);
        omega       104.72;
    }
    \endverbatim


SourceFiles
    MRFSource.C

\*---------------------------------------------------------------------------*/

#ifndef MRFSource_H
#define MRFSource_H

#include "fvOption.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class MRFZone;

namespace fv
{

/*---------------------------------------------------------------------------*\
                          Class MRFSource Declaration
\*---------------------------------------------------------------------------*/

class MRFSource
:
    public option
{

protected:

    // Protected data

        //- Run-time selectable porosity model
        autoPtr<MRFZone> mrfPtr_;

        //- Velocity field name, default = U
        word UName_;


    // Protected Member Functions

        //- Initialise
        void initialise();


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        MRFSource(const MRFSource&);

        //- Disallow default bitwise assignment
        void operator=(const MRFSource&);


public:

    //- Runtime type information
    TypeName("MRFSource");


    // Constructors

        //- Construct from components
        MRFSource
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~MRFSource()
    {}


    // Member Functions

        // Add explicit and implicit contributions

            //- Vector
            virtual void addSup
            (
                fvMatrix<vector>& eqn,
                const label fieldI
            );


        // Add explicit and implicit contributions to compressible equation

            //- Vector
            virtual void addSup
            (
                const volScalarField& rho,
                fvMatrix<vector>& eqn,
                const label fieldI
            );


        // Flux manipulations

            //- Make the given absolute flux relative
            virtual void makeRelative(surfaceScalarField& phi) const;

            //- Make the given absolute boundary flux relative
            virtual void makeRelative
            (
                FieldField<fvsPatchField, scalar>& phi
            ) const;

            //- Make the given absolute mass-flux relative
            virtual void makeRelative
            (
                const surfaceScalarField& rho,
                surfaceScalarField& phi
            ) const;

            //- Make the given relative flux absolute
            virtual void makeAbsolute(surfaceScalarField& phi) const;

            //- Make the given relative mass-flux absolute
            virtual void makeAbsolute
            (
                const surfaceScalarField& rho,
                surfaceScalarField& phi
            ) const;


        // I-O

            //- Write data
            virtual void writeData(Ostream&) const;

            //- Read dictionary
            virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
