module QuantumSystemTemplates

"""
Quantum system templates for common physical systems.

Organized by physical platform:

## Ions
- **IonChainSystem**: Trapped ion chains with motional modes
- **MolmerSorensenCoupling**: MS gate coupling helpers
- **RadialMSGateSystem**: Radial-mode MS gates (IEEE TQE 2024)

## Atoms  
- **RydbergChainSystem**: Rydberg atom chains with blockade

## Transmons
- **TransmonSystem**: Single transmon with anharmonicity
- **MultiTransmonSystem**: Coupled transmon systems
- **TransmonCavitySystem**: Transmons coupled to cavities

## Cats
- **CatSystem**: Quantum cat states in cavities
"""

using LinearAlgebra
using TestItems

# Import from parent PiccoloQuantumObjects modules
using ..QuantumSystems:
    QuantumSystem,
    OpenQuantumSystem,
    CompositeQuantumSystem,
    AbstractQuantumSystem,
    lift_operator
using ..QuantumObjectUtils: annihilate, operator_from_string, PAULIS
using ..Gates: GATES

const ⊗ = kron

# Ion systems
include("ions/ion_chain.jl")
include("ions/radial_ms.jl")

# Atom systems
include("atoms/rydberg_chain.jl")

# Transmon systems
include("transmons/transmon_system.jl")

# Cat systems
include("cats/cat_system.jl")

end  # module QuantumSystemTemplates
