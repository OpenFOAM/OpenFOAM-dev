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

#include "molecule.H"
#include "IOstreams.H"
#include "moleculeCloud.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const std::size_t Foam::molecule::sizeofFields_
(
    offsetof(molecule, siteForces_) - offsetof(molecule, Q_)
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::molecule::molecule
(
    const polyMesh& mesh,
    Istream& is,
    bool readFields
)
:
    particle(mesh, is, readFields),
    Q_(Zero),
    v_(Zero),
    a_(Zero),
    pi_(Zero),
    tau_(Zero),
    specialPosition_(Zero),
    potentialEnergy_(0.0),
    rf_(Zero),
    special_(0),
    id_(0),
    siteForces_(0),
    sitePositions_(0)
{
    if (readFields)
    {
        if (is.format() == IOstream::ASCII)
        {
            is  >> Q_;
            is  >> v_;
            is  >> a_;
            is  >> pi_;
            is  >> tau_;
            is  >> specialPosition_;
            potentialEnergy_ = readScalar(is);
            is  >> rf_;
            special_ = readLabel(is);
            id_ = readLabel(is);
            is  >> siteForces_;
            is  >> sitePositions_;
        }
        else
        {
            is.read(reinterpret_cast<char*>(&Q_), sizeofFields_);
            is  >> siteForces_ >> sitePositions_;
        }
    }

    // Check state of Istream
    is.check
    (
        "Foam::molecule::molecule"
        "(const Cloud<molecule>& cloud, Foam::Istream&), bool"
    );
}


void Foam::molecule::readFields(Cloud<molecule>& mC)
{
    bool valid = mC.size();

    particle::readFields(mC);

    IOField<tensor> Q(mC.fieldIOobject("Q", IOobject::MUST_READ), valid);
    mC.checkFieldIOobject(mC, Q);

    IOField<vector> v(mC.fieldIOobject("v", IOobject::MUST_READ), valid);
    mC.checkFieldIOobject(mC, v);

    IOField<vector> a(mC.fieldIOobject("a", IOobject::MUST_READ), valid);
    mC.checkFieldIOobject(mC, a);

    IOField<vector> pi(mC.fieldIOobject("pi", IOobject::MUST_READ), valid);
    mC.checkFieldIOobject(mC, pi);

    IOField<vector> tau(mC.fieldIOobject("tau", IOobject::MUST_READ), valid);
    mC.checkFieldIOobject(mC, tau);

    IOField<vector> specialPosition
    (
        mC.fieldIOobject("specialPosition", IOobject::MUST_READ),
        valid
    );
    mC.checkFieldIOobject(mC, specialPosition);

    IOField<label> special
    (
        mC.fieldIOobject("special", IOobject::MUST_READ),
        valid
    );
    mC.checkFieldIOobject(mC, special);

    IOField<label> id(mC.fieldIOobject("id", IOobject::MUST_READ), valid);
    mC.checkFieldIOobject(mC, id);

    label i = 0;
    forAllIter(moleculeCloud, mC, iter)
    {
        molecule& mol = iter();

        mol.Q_ = Q[i];
        mol.v_ = v[i];
        mol.a_ = a[i];
        mol.pi_ = pi[i];
        mol.tau_ = tau[i];
        mol.specialPosition_ = specialPosition[i];
        mol.special_ = special[i];
        mol.id_ = id[i];
        i++;
    }
}


void Foam::molecule::writeFields(const Cloud<molecule>& mC)
{
    particle::writeFields(mC);

    label np = mC.size();

    IOField<tensor> Q(mC.fieldIOobject("Q", IOobject::NO_READ), np);
    IOField<vector> v(mC.fieldIOobject("v", IOobject::NO_READ), np);
    IOField<vector> a(mC.fieldIOobject("a", IOobject::NO_READ), np);
    IOField<vector> pi(mC.fieldIOobject("pi", IOobject::NO_READ), np);
    IOField<vector> tau(mC.fieldIOobject("tau", IOobject::NO_READ), np);
    IOField<vector> specialPosition
    (
        mC.fieldIOobject("specialPosition", IOobject::NO_READ),
        np
    );
    IOField<label> special(mC.fieldIOobject("special", IOobject::NO_READ), np);
    IOField<label> id(mC.fieldIOobject("id", IOobject::NO_READ), np);

    // Post processing fields

    IOField<vector> piGlobal
    (
        mC.fieldIOobject("piGlobal", IOobject::NO_READ),
        np
    );

    IOField<vector> tauGlobal
    (
        mC.fieldIOobject("tauGlobal", IOobject::NO_READ),
        np
    );

    IOField<vector> orientation1
    (
        mC.fieldIOobject("orientation1", IOobject::NO_READ),
        np
    );

    IOField<vector> orientation2
    (
        mC.fieldIOobject("orientation2", IOobject::NO_READ),
        np
    );

    IOField<vector> orientation3
    (
        mC.fieldIOobject("orientation3", IOobject::NO_READ),
        np
    );

    label i = 0;
    forAllConstIter(moleculeCloud, mC, iter)
    {
        const molecule& mol = iter();

        Q[i] = mol.Q_;
        v[i] = mol.v_;
        a[i] = mol.a_;
        pi[i] = mol.pi_;
        tau[i] = mol.tau_;
        specialPosition[i] = mol.specialPosition_;
        special[i] = mol.special_;
        id[i] = mol.id_;

        piGlobal[i] = mol.Q_ & mol.pi_;
        tauGlobal[i] = mol.Q_ & mol.tau_;

        orientation1[i] = mol.Q_ & vector(1,0,0);
        orientation2[i] = mol.Q_ & vector(0,1,0);
        orientation3[i] = mol.Q_ & vector(0,0,1);

        i++;
    }

    const bool valid = np > 0;

    Q.write(valid);
    v.write(valid);
    a.write(valid);
    pi.write(valid);
    tau.write(valid);
    specialPosition.write(valid);
    special.write(valid);
    id.write(valid);

    piGlobal.write(valid);
    tauGlobal.write(valid);

    orientation1.write(valid);
    orientation2.write(valid);
    orientation3.write(valid);

    Info<< "writeFields " << mC.name() << endl;

    if (isA<moleculeCloud>(mC))
    {
        const moleculeCloud& m = dynamic_cast<const moleculeCloud&>(mC);

        m.writeXYZ
        (
            m.mesh().time().timePath()/cloud::prefix/"moleculeCloud.xmol"
        );
    }
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const molecule& mol)
{
    if (os.format() == IOstream::ASCII)
    {
        os  << token::SPACE << static_cast<const particle&>(mol)
            << token::SPACE << mol.Q_
            << token::SPACE << mol.v_
            << token::SPACE << mol.a_
            << token::SPACE << mol.pi_
            << token::SPACE << mol.tau_
            << token::SPACE << mol.specialPosition_
            << token::SPACE << mol.potentialEnergy_
            << token::SPACE << mol.rf_
            << token::SPACE << mol.special_
            << token::SPACE << mol.id_
            << token::SPACE << mol.siteForces_
            << token::SPACE << mol.sitePositions_;
    }
    else
    {
        os  << static_cast<const particle&>(mol);
        os.write
        (
            reinterpret_cast<const char*>(&mol.Q_),
            molecule::sizeofFields_
        );
        os  << mol.siteForces_ << mol.sitePositions_;
    }

    // Check state of Ostream
    os.check
    (
        "Foam::Ostream& Foam::operator<<"
        "(Foam::Ostream&, const Foam::molecule&)"
    );

    return os;
}


// ************************************************************************* //
