/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2021 hyStrath
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of hyStrath, a derivative work of OpenFOAM.

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

#include "basic2Thermo.H"
#include "stringOps.H"
#include "wordIOList.H"
#include "zeroGradientFvPatchFields.H"
#include "fixedEnergyFvPatchScalarField.H"
#include "gradientEnergyFvPatchScalarField.H"
#include "mixedEnergyFvPatchScalarField.H"
#include "fixedJumpFvPatchFields.H"
#include "fixedJumpAMIFvPatchFields.H"
#include "energyJumpFvPatchScalarField.H"
#include "energyJumpAMIFvPatchScalarField.H"


/* * * * * * * * * * * * * * * private static data * * * * * * * * * * * * * */

namespace Foam
{
    defineTypeNameAndDebug(basic2Thermo, 0);
    defineRunTimeSelectionTable(basic2Thermo, fvMesh);

    defineRunTimeSelectionTable(basic2Thermo, fvMeshDictPhase);
}

const Foam::word Foam::basic2Thermo::dictName("thermophysicalProperties");

//hystrath only uses the 7 component option
const Foam::wordList Foam::basic2Thermo::componentHeader7
({
    "type",
    "mixture",
    "transport",
    "thermo",
    "equationOfState",
    "specie",
    "energy"
});


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::Ostream& Foam::basic2Thermo::printThermoNames
(
    Ostream& os,
    const wordList& cmptNames,
    const wordList& thermoNames
)
{
    const int nCmpt = cmptNames.size();

    // Build a table of constituent parts by split name into constituent parts
    // - remove incompatible entries from the list
    // - note: row-0 contains the names of constituent parts (ie, the header)

    DynamicList<wordList> outputTbl;
    outputTbl.resize(thermoNames.size()+1);

    label rowi = 0;

    // Header
    outputTbl[rowi] = cmptNames;
    if (!outputTbl[rowi].empty())
    {
        ++rowi;
    }

    for (const word& thermoName : thermoNames)
    {
        outputTbl[rowi] = basic2Thermo::splitThermoName(thermoName, nCmpt);
        if (!outputTbl[rowi].empty())
        {
            ++rowi;
        }
    }

    if (rowi > 1)
    {
        outputTbl.resize(rowi);
        Foam::printTable(outputTbl, os);
    }

    return os;
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::word Foam::basic2Thermo::makeThermoName
(
    const dictionary& thermoTypeDict,
    const wordList*& cmptHeaderPtr
)
{
    if (cmptHeaderPtr)
    {
        cmptHeaderPtr = &(componentHeader7);
    }

    return word
    (
          thermoTypeDict.get<word>("type") + '<'
        + thermoTypeDict.get<word>("mixture") + '<'
        + thermoTypeDict.get<word>("transport") + '<'
        + thermoTypeDict.get<word>("thermo") + '<'
        + thermoTypeDict.get<word>("equationOfState") + '<'
        + thermoTypeDict.get<word>("specie") + ">>,"
        + thermoTypeDict.get<word>("energy") + ">>>"
    );
    
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::volScalarField& Foam::basic2Thermo::lookupOrConstruct
(
    const fvMesh& mesh,
    const word& fieldName,
    bool& isOwner
)
{
    auto* ptr = mesh.objectRegistry::getObjectPtr<volScalarField>(fieldName);

    isOwner = !ptr;

    if (!ptr)
    {
        ptr = new volScalarField
        (
            IOobject
            (
                fieldName,
                mesh.time().timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE,
                IOobject::REGISTER
            ),
            mesh
        );

        // Transfer ownership of this object to the objectRegistry
        ptr->store();
    }

    return *ptr;
}



Foam::basic2Thermo::basic2Thermo
(
    const fvMesh& mesh,
    const word& phaseName
)
:
    IOdictionary
    (
        IOobject
        (
            phasePropertyName(dictName, phaseName),
            mesh.time().constant(),
            mesh,
            IOobject::READ_MODIFIED,
            IOobject::NO_WRITE,
            IOobject::REGISTER
        )
    ),

    phaseName_(phaseName),

    pOwner_(false),
    TOwner_(false),

    p_(lookupOrConstruct(mesh, "p", pOwner_)),

    pe_
    (
        IOobject
        (
            "pe",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimPressure
    ),
    
    T_
    (
        IOobject
        (
            phasePropertyName("Tt"),
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimTemperature
    )
{
    this->readIfPresent("updateT", TOwner_);  // Manual override
}


Foam::basic2Thermo::basic2Thermo
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& phaseName
)
:
    IOdictionary
    (
        IOobject
        (
            phasePropertyName(dictName, phaseName),
            mesh.time().constant(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::REGISTER
        ),
        dict
    ),

    phaseName_(phaseName),

    pOwner_(false),
    TOwner_(false),

    p_(lookupOrConstruct(mesh, "p", pOwner_)),
    
    pe_
    (
        IOobject
        (
            phasePropertyName("pe"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("pe", dimPressure, 0.0)
    ),

    T_
    (
        IOobject
        (
            phasePropertyName("Tt"),
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimTemperature
    )
{
    this->readIfPresent("updateT", TOwner_);  // Manual override
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::basic2Thermo> Foam::basic2Thermo::New
(
    const fvMesh& mesh,
    const word& phaseName
)
{
    return New<basic2Thermo>(mesh, phaseName);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::basic2Thermo::~basic2Thermo()
{
    if (pOwner_)
    {
        db().checkOut(p_.name());
    }

    if (TOwner_)
    {
        db().checkOut(T_.name());
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


const Foam::basic2Thermo& Foam::basic2Thermo::lookupThermo
(
    const fvPatchScalarField& pf
)
{
    return pf.db().lookupObject<basic2Thermo>(dictName);   
}


void Foam::basic2Thermo::validate
(
    const string& app,
    const word& a
) const
{
    if (!("e" == phasePropertyName(a)))
    {
        FatalErrorIn(app)
            << "Supported energy type is " << phasePropertyName(a)
            << ", thermodynamics package provides " << "e"
            << exit(FatalError);
    }
}

void Foam::basic2Thermo::validate
(
    const string& app,
    const word& a,
    const word& b
) const
{
    if
    (
       !(
            "e" == phasePropertyName(a)
         || "e" == phasePropertyName(b)
        )
    )
    {
        FatalErrorIn(app)
            << "Supported energy types are " << phasePropertyName(a)
            << " and " << phasePropertyName(b)
            << ", thermodynamics package provides " << "e"
            << exit(FatalError);
    }
}

void Foam::basic2Thermo::validate
(
    const string& app,
    const word& a,
    const word& b,
    const word& c
) const
{
    if
    (
       !(
            "e" == phasePropertyName(a)
         || "e" == phasePropertyName(b)
         || "e" == phasePropertyName(c)
        )
    )
    {
        FatalErrorIn(app)
            << "Supported energy types are " << phasePropertyName(a)
            << ", " << phasePropertyName(b)
            << " and " << phasePropertyName(c)
            << ", thermodynamics package provides " << "e"
            << exit(FatalError);
    }
}

void Foam::basic2Thermo::validate
(
    const string& app,
    const word& a,
    const word& b,
    const word& c,
    const word& d
) const
{
    if
    (
       !(
            "e" == phasePropertyName(a)
         || "e" == phasePropertyName(b)
         || "e" == phasePropertyName(c)
         || "e" == phasePropertyName(d)
        )
    )
    {
        FatalErrorIn(app)
            << "Supported energy types are " << phasePropertyName(a)
            << ", " << phasePropertyName(b)
            << ", " << phasePropertyName(c)
            << " and " << phasePropertyName(d)
            << ", thermodynamics package provides " << "e"
            << exit(FatalError);
    }
}


Foam::wordList Foam::basic2Thermo::splitThermoName
(
    const std::string& thermoName,
    const int nExpectedCmpts
)
{
    /*
    wordList cmpts(nCmpt);

    string::size_type beg=0, end=0, endb=0, endc=0;
    int i = 0;

    while
    (
        (endb = thermoName.find('<', beg)) != string::npos
     || (endc = thermoName.find(',', beg)) != string::npos
    )
    {
        if (endb == string::npos)
        {
            end = endc;
        }
        else if ((endc = thermoName.find(',', beg)) != string::npos)
        {
            end = min(endb, endc);
        }
        else
        {
            end = endb;
        }

        if (beg < end)
        {
            cmpts[i] = thermoName.substr(beg, end-beg);
            cmpts[i++].replaceAll(">","");
        }
        beg = end + 1;
    }

    if (beg < thermoName.size())
    {
        cmpts[i] = thermoName.substr(beg, string::npos);
        cmpts[i++].replaceAll(">","");
    }

    return cmpts;
    */
    // Split on ",<>" but include space for good measure.
    // Splits things like
    // "hePsiThermo<pureMixture<const<hConst<perfectGas<specie>>,enthalpy>>>"

    const auto parsed = stringOps::splitAny<std::string>(thermoName, " ,<>");
    const int nParsed(parsed.size());

    wordList cmpts;

    if (!nExpectedCmpts || nParsed == nExpectedCmpts)
    {
        cmpts.resize(nParsed);

        auto iter = cmpts.begin();
        for (const auto& sub : parsed)
        {
            *iter = word(sub.str());
            ++iter;
        }
    }

    return cmpts;
}


const Foam::volScalarField& Foam::basic2Thermo::p() const
{
    return p_;
}


Foam::volScalarField& Foam::basic2Thermo::p()
{
    return p_;
}


const Foam::volScalarField& Foam::basic2Thermo::pe() const
{
    return pe_;
}


Foam::volScalarField& Foam::basic2Thermo::pe()
{
    return pe_;
}


const Foam::volScalarField& Foam::basic2Thermo::T() const
{
    return T_;
}


Foam::volScalarField& Foam::basic2Thermo::T()
{
    return T_;
}


bool Foam::basic2Thermo::read()
{
    return regIOobject::read();
}


// ************************************************************************* //
