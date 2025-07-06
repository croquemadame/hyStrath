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

Class
    Foam::basic2Thermo

Description
    Abstract base-class for fluid and solid thermodynamic properties

SourceFiles
    basic2Thermo.C

\*---------------------------------------------------------------------------*/

#ifndef basic2Thermo_H
#define basic2Thermo_H

#include "volFields.H"
#include "typeInfo.H"
#include "IOdictionary.H"
#include "autoPtr.H"
#include "wordIOList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class basic2Thermo Declaration
\*---------------------------------------------------------------------------*/

class basic2Thermo
:
    public IOdictionary
{
    // Private Data

        //- Components names/order
        static const wordList componentHeader4;

        //- Components names/order
        static const wordList componentHeader7;


    // Private Member Functions

        //- Construct name of thermo package from dictionary components
        static word makeThermoName
        (
            const dictionary& dict,
            const wordList*& cmptHeaderPtr
        );

        //- Look up field from registry or construct and store
        static volScalarField& lookupOrConstruct
        (
            const fvMesh& mesh,
            const word& fieldName,
            bool& isOwner  //!< Stored to registry by this instance
        );

        //- Generic lookup for thermodynamics package thermoTypeName
        //  \return constructor pointer, or FatalError
        template<class Thermo, class ThermoConstructTable>
        static typename ThermoConstructTable::mapped_type
        getThermoOrDie
        (
            const dictionary& thermoTypeDict,
            ThermoConstructTable& thermoTable,
            const word& thermoTypeName,
            const wordList& cmptNames
        );

        //- Generic lookup for each of the related thermodynamics packages
        //  \return constructor pointer, or FatalError
        template<class Thermo, class ThermoConstructTable>
        static typename ThermoConstructTable::mapped_type
        getThermoOrDie
        (
            const dictionary& thermoDict,
            ThermoConstructTable& thermoTable
        );




protected:

    // Protected data

        //- Phase-name
        const word& phaseName_;

        //including these just in case
        //- Pressure created and stored by this instance
        bool pOwner_;

        //- Temperature created and stored by this instance
        bool TOwner_;

        //- Include dpdt term in the enthalpy equation?
        bool dpdt_;

        // Fields

            //- Pressure [Pa]
            volScalarField& p_;
            
            //- Electron pressure [Pa]
            volScalarField pe_;

            //- Trans-rotational temperature [K]
            volScalarField T_;


    // Protected Member Functions

        //- Return the enthalpy/internal energy field boundary types
        //- by interrogating the temperature field boundary types
        //wordList heBoundaryTypes();

        //- Return the enthalpy/internal energy field boundary base types
        //- by interrogating the temperature field boundary types
        //wordList heBoundaryBaseTypes();

        //- No copy construct
        //basicThermo(const basicThermo&) = delete;

        //- Construct as copy (not implemented)
        basic2Thermo(const basic2Thermo&);

        volScalarField& lookupOrConstruct
        (
            const fvMesh& mesh,
            const char* name
        ) const;


public:

    //- Runtime type information
    TypeName("basic2Thermo");


    //- Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        basic2Thermo,
        fvMesh,
        (const fvMesh& mesh, const word& phaseName),
        (mesh, phaseName)
    );

    //- Declare run-time constructor selection table for dictionary based
    declareRunTimeSelectionTable
    (
        autoPtr,
        basic2Thermo,
        fvMeshDictPhase,
        (const fvMesh& mesh, const word& phaseName, const word& dictName),
        (mesh, phaseName, dictName)
    );


    // Constructors

        //- Construct from mesh and phase name
        basic2Thermo
        (
            const fvMesh&,
            const word& phaseName
        );

        //- Construct from mesh, dictionary and phase name
        basic2Thermo
        (
            const fvMesh&,
            const dictionary&,
            const word& phaseName
        );

        //- Construct from mesh, phase name and explicit naming of the
        //- dictionary (so it can be shared amongst phases).
        //basicThermo
        //(
        //    const fvMesh&,
        //    const word& phaseName,
        //    const word& dictionaryName
        //);


    // Selectors
        /*
        //- Generic lookup for each of the related thermodynamics packages
        template<class Thermo, class Table>
        static typename Table::iterator lookupThermo
        (
            const dictionary& thermoDict,
            Table* tablePtr
        );

        //- Generic New for each of the related thermodynamics packages
        template<class Thermo>
        static autoPtr<Thermo> New
        (
            const fvMesh&,
            const word& phaseName=word::null
        );

        //- Generic New for each of the related thermodynamics packages
        template<class Thermo>
        static autoPtr<Thermo> New
        (
            const fvMesh&,
            const dictionary&,
            const word& phaseName=word::null
        );

        //- Specialisation of the Generic New for basic2Thermo
        static autoPtr<basic2Thermo> New
        (
            const fvMesh&,
            const word& phaseName=word::null
        );

        */

        //- Generic New for each of the related thermodynamics packages
        template<class Thermo>
        static autoPtr<Thermo> New
        (
            const fvMesh&,
            const word& phaseName = word::null
        );

        //- Generic New for each of the related thermodynamics packages
        template<class Thermo>
        static autoPtr<Thermo> New
        (
            const fvMesh&,
            const dictionary&,
            const word& phaseName = word::null
        );

        //- Generic New for each of the related thermodynamics packages
        template<class Thermo>
        static autoPtr<Thermo> New
        (
            const fvMesh&,
            const word& phaseName,
            const word& dictName
        );

        //- Specialisation of the Generic New for basicThermo
        static autoPtr<basic2Thermo> New
        (
            const fvMesh&,
            const word& phaseName = word::null
        );

    //- Destructor
    virtual ~basic2Thermo();


    // Member functions

        static const word dictName;

        static word phasePropertyName
        (
            const word& name,
            const word& phaseName
        )
        {
            return IOobject::groupName(name, phaseName);
        }

        word phasePropertyName(const word& name) const
        {
            return basic2Thermo::phasePropertyName(name, phaseName_);
        }

        //- The phase property name scoped as scope:property.phase
        //- \sa IOobject::groupName and \sa IOobject::scopedName
        word phaseScopedName(const std::string& scope, const word& name) const
        {
            return IOobject::groupName
            (
                IOobject::scopedName(scope, name),
                phaseName_
            );
        }


        static const basic2Thermo& lookupThermo(const fvPatchScalarField& pf);

        //- Print (filtered) table of thermo names, splits on \c " ,<>"
        static Ostream& printThermoNames
        (
            Ostream& os,
            const wordList& cmptNames,
            const wordList& thermoNames
        );


        //- Check that the thermodynamics package is consistent
        //  with energy forms supported by the application
        void validate
        (
            const string& app,
            const word&
        ) const;

        //- Check that the thermodynamics package is consistent
        //  with energy forms supported by the application
        void validate
        (
            const string& app,
            const word&,
            const word&
        ) const;

        //- Check that the thermodynamics package is consistent
        //  with energy forms supported by the application
        void validate
        (
            const string& app,
            const word&,
            const word&,
            const word&
        ) const;

        //- Check that the thermodynamics package is consistent
        //  with energy forms supported by the application
        void validate
        (
            const string& app,
            const word&,
            const word&,
            const word&,
            const word&
        ) const;

        //- Split name of thermo package into a list of the components names
        static wordList splitThermoName
        (
            const std::string& thermoName,
            const int nExpectedCmpts
        );

        //- Update properties
        virtual void correct() = 0;

        //- Return the name of the thermo physics
        //virtual word thermoName() const = 0;

        //- Return true if the equation of state is incompressible
        //  i.e. rho != f(p)
        //virtual bool incompressible() const = 0;

        //- Return true if the equation of state is isochoric
        //  i.e. rho = const
        //virtual bool isochoric() const = 0;

        //- True if dpdt term should be included in enthalpy equation
        bool dpdt() const noexcept
        {
            return dpdt_;
        }

        //- True if T should be updated
        bool updateT() const noexcept
        {
            return TOwner_;
        }

        // Access to thermodynamic state variables

            //- Pressure [Pa]
            virtual const volScalarField& p() const;
            
            //- Non-const access to pressure [Pa]
            virtual volScalarField& p();

            //- Electron pressure [Pa]
            virtual const volScalarField& pe() const;
            
            //- Non-const access to electron pressure [Pa]
            virtual volScalarField& pe();

            //- Density [kg/m^3]
            virtual tmp<volScalarField> rho() const = 0;

            //- Density for patch [kg/m^3]
            virtual tmp<scalarField> rho(const label patchi) const = 0;

            //- Chemical enthalpy [J/kg]
            virtual tmp<volScalarField> hc() const = 0;

        // Fields derived from thermodynamic state variables

            //- Trans-rotational temperature [K]
            virtual const volScalarField& T() const;

            //- Non-const access to trans-rotational temperature [K]
            virtual volScalarField& T();

            //- Heat capacity at constant pressure [J/kg/K]
            //  TODO
            virtual tmp<volScalarField> Cp() const = 0;

            //- Heat capacity at constant pressure [J/kg/K] for patch
            //  TODO
            virtual tmp<scalarField> Cp
            (
                const scalarField& p,
                const scalarField& T,
                const label patchi
            ) const = 0;

            //- Heat capacity at constant volume [J/kg/K]
            //  TODO
            virtual tmp<volScalarField> Cv() const = 0;

            //- Heat capacity ratio, gamma = Cp/Cv
            //  TODO
            virtual tmp<volScalarField> gamma() const = 0;

        // Fields derived from transport state variables

            //- Effective thermal diffusivity for temperature
            //  of mixture [J/m/s/K]
            virtual tmp<volScalarField> kappaEff
            (
                const volScalarField&
            ) const = 0;

            //- Effective thermal diffusivity for temperature
            //  of mixture for patch [J/m/s/K]
            virtual tmp<scalarField> kappaEff
            (
                const scalarField& alphat,
                const label patchi
            ) const = 0;

            //- Effective thermal diffusivity of mixture [kg/m/s]
            virtual tmp<volScalarField> alphaEff
            (
                const volScalarField& alphat
            ) const = 0;

            //- Effective thermal diffusivity of mixture for patch [kg/m/s]
            virtual tmp<scalarField> alphaEff
            (
                const scalarField& alphat,
                const label patchi
            ) const = 0;


        //- Read thermophysical properties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "basic2ThermoTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
