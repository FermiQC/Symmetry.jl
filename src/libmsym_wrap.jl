@enum msym_geometry begin
        MSYM_GEOMETRY_UNKNOWN = 0
        MSYM_GEOMETRY_SPHERICAL = 1
        MSYM_GEOMETRY_LINEAR = 2
        MSYM_GEOMETRY_PLANAR_REGULAR = 3
        MSYM_GEOMETRY_PLANAR_IRREGULAR = 4
        MSYM_GEOMETRY_POLYHEDRAL_PROLATE = 5
        MSYM_GEOMETRY_POLYHEDRAL_OBLATE = 6
        MSYM_GEOMETRY_ASSYMETRIC = 7
end
    
struct msym_symmetry_operation
    @enum msym_symmetry_operation_type begin
        MSYM_SYMMETRY_OPERATION_TYPE_IDENTITY = 0
        MSYM_SYMMETRY_OPERATION_TYPE_PROPER_ROTATION = 1
        MSYM_SYMMETRY_OPERATION_TYPE_IMPROPER_ROTATION = 2
        MSYM_SYMMETRY_OPERATION_TYPE_REFLECTION = 3
        MSYM_SYMMETRY_OPERATION_TYPE_INVERSION = 4
    end
    order::Cint                              # Order of proper/improper rotation
    power::Cint                              # Power (e.g. C2^2 = I)
    @enum msym_symmetry_operation_orientation begin
        MSYM_SYMMETRY_OPERATION_ORIENTATION_NONE = 0
        MSYM_SYMMETRY_OPERATION_ORIENTATION_HORIZONTAL = 1
        MSYM_SYMMETRY_OPERATION_ORIENTATION_VERTICAL = 2
        MSYM_SYMMETRY_OPERATION_ORIENTATION_DIHEDRAL = 3
    end
    v::NTuple{3,Cdouble}                            # Proper/improper rotation vector or reflection plane normal
    cla::Cint                                # Class of symmetry operation (point group dependant)
end
    
@enum msym_point_group_type begin
    MSYM_POINT_GROUP_TYPE_Kh  = 0
    MSYM_POINT_GROUP_TYPE_K   = 1
    MSYM_POINT_GROUP_TYPE_Ci  = 2
    MSYM_POINT_GROUP_TYPE_Cs  = 3
    MSYM_POINT_GROUP_TYPE_Cn  = 4
    MSYM_POINT_GROUP_TYPE_Cnh = 5
    MSYM_POINT_GROUP_TYPE_Cnv = 6
    MSYM_POINT_GROUP_TYPE_Dn  = 7
    MSYM_POINT_GROUP_TYPE_Dnh = 8
    MSYM_POINT_GROUP_TYPE_Dnd = 9
    MSYM_POINT_GROUP_TYPE_Sn  = 10
    MSYM_POINT_GROUP_TYPE_T   = 11
    MSYM_POINT_GROUP_TYPE_Td  = 12
    MSYM_POINT_GROUP_TYPE_Th  = 13
    MSYM_POINT_GROUP_TYPE_O   = 14
    MSYM_POINT_GROUP_TYPE_Oh  = 15
    MSYM_POINT_GROUP_TYPE_I   = 16
    MSYM_POINT_GROUP_TYPE_Ih  = 17
end

struct msym_subgroup
    type::msym_point_group_type
    n::Cint
    order::Cint
    primary::Ref{msym_symmetry_operation}
    #msym_symmetry_operation_t *primary;
    sops::Ref{Ref{msym_symmetry_operation}}
    #msym_symmetry_operation_t **sops;
    generators::NTuple{2,Ref{msym_subgroup}}
    name::NTuple{8,Cchar}
end

struct msym_thresholds
    zero::Cdouble                            # For determining if something is zero (e.g. vectors close to center of mass)
    geometry::Cdouble                        # For translating inertial tensor eigenvalues to geometric structures
    angle::Cdouble                           # For determining angles, (e.g. if vectors are parallel)
    equivalence::Cdouble                     # Equivalence test threshold
    eigfact::Cdouble                         # Jacobi eigenvalue algorithm threshold
    permutation::Cdouble                     # Equality test when determining permutation for symmetry operation
    orthogonalization::Cdouble               # For orthogonalizing orbital subspaces
end

struct msym_element 
    id::Ref{Cvoid}                               # custom identifier
    m::Cdouble                              # Mass
    v::NTuple{3,Cdouble}                  # Position
    n::Cint                                 # Nuclear charge
    name::NTuple{4,Cchar}                 # Name
end

struct msym_equivalence_set
    elements::Ref{Ref{msym_element}}             # Pointers to elements
    err::Cdouble                            # Maximum error when detecting this equivalence set
    length::Cint                            # Number of elements
end

struct msym_real_spherical_harmonic
    n::Cint                                # Principal
    l::Cint                                # Azimuthal
    m::Cint                                # Liniear combination of magnetic quantum number (e.g. 2pz = 0, 2px = 1, 2py = -1)
end

struct msym_basis_function
    id::Ref{Cvoid}                               # custom identifier
    @enum msym_basis_type begin
        MSYM_BASIS_TYPE_REAL_SPHERICAL_HARMONIC = 0
        MSYM_BASIS_TYPE_CARTESIAN = 1
    end
    element::Ref{msym_element}
    rsh::msym_real_spherical_harmonic       # Atomic orbital basis
    name::NTuple{8,Cchar}
end

struct msym_partner_function
    i::Cint          # index of partner 0
    d::Cint          # component (dimension)
end

struct msym_salc
    d::Cint             # dimension of space (same as msym_character_table_t.s[msym_subrepresentation_space_t.s].d)
    fl::Cint            # number of basis functions
    pf::Ref{Cvoid}           # partner functions double[d][fl]
    f::Ref{Ref{msym_basis_function}}
end

struct msym_subrepresentation_space
    s::Cint      # symmetry species
    salcl::Cint  # nr of SALCs
    salc::Ref{msym_salc}
end

struct msym_symmetry_species
    d::Cint # dimensionality of a (ir)reducible representation of this species
    r::Cint # sum over all x ct->classc[x]*ct->table[i][x]^2/pg->order (can be decomposed into r irreducible representations)
    name::NTuple{8,Cchar}
end

struct msym_character_table
    d::Cint
    classc::Ref{Cint}
    sops::Ref{Ref{msym_symmetry_operation}}
    s::Ref{msym_symmetry_species}
    table::Ref{Cvoid}  #double[d][d]
end

struct msym_permutation_cycle_t
    l::Cint
    s::Cint
end

struct msym_permutation
    p::Ref{Cint}
    p_length::Cint
    c::Ref{msym_permutation_cycle_t}
    c_length::Cint
end

struct msym_context
    thresholds::Ref{msym_thresholds}
    elements::Ref{msym_element}
    pelements::Ref{Ref{msym_element}}
    basis::Ref{msym_basis_function}
    es::Ref{msym_equivalence_set}
    es_perm::Ref{Ref{msym_permutation}}
    srs::Ref{msym_subrepresentation_space}
    srsbf::Ref{Ref{msym_basis_function}}
    srs_span::Ref{Cint}
    flags::Culong
    elementsl::Cint
    basisl::Cint
    esl::Cint
    srsl::Cint
    es_perml::Cint
    sgl::Cint
    pg::Ref{msym_point_group_type}
    sg::Ref{msym_subgroup}
    cm::NTuple{3,Cdouble}
    geometry::msym_geometry
    eigval::NTuple{3,Cdouble}
    eigvec::NTuple{3,NTuple{Cdouble}}
    struct _external_data
        eesmap::Ref{Ref{msym_equivalence_set}}
        set_elements_ptr::Ref{msym_element}
        elements::Ref{msym_element}
        es::Ref{msym_equivalence_set}
    end
end

ctx = ccall((:msymCreateContext, "libmsym.so"), msym_context, ())