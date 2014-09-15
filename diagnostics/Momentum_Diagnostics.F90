!    Copyright (C) 2006 Imperial College London and others.
!    
!    Please see the AUTHORS file in the main source directory for a full list
!    of copyright holders.
!
!    Prof. C Pain
!    Applied Modelling and Computation Group
!    Department of Earth Science and Engineering
!    Imperial College London
!
!    amcgsoftware@imperial.ac.uk
!    
!    This library is free software; you can redistribute it and/or
!    modify it under the terms of the GNU Lesser General Public
!    License as published by the Free Software Foundation,
!    version 2.1 of the License.
!
!    This library is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!    Lesser General Public License for more details.
!
!    You should have received a copy of the GNU Lesser General Public
!    License along with this library; if not, write to the Free Software
!    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
!    USA

#include "fdebug.h"

module momentum_diagnostics

  use boundary_conditions
  use coriolis_module, only : two_omega => coriolis
  use diagnostic_source_fields
  use field_derivatives
  use field_options
  use fields
  use fldebug
  use geostrophic_pressure
  use global_parameters, only : OPTION_PATH_LEN
  use multimaterial_module
  use solvers
  use sparsity_patterns_meshes
  use spud
  use state_fields_module
  use state_module
  use fefields, only : compute_cv_mass
  use sediment, only : get_n_sediment_fields, get_sediment_item
  
  implicit none
  
  private
  
  public :: calculate_strain_rate, calculate_bulk_viscosity, calculate_strain_rate_second_invariant, &
            calculate_sediment_concentration_dependent_viscosity, &
            calculate_buoyancy, calculate_coriolis, calculate_tensor_second_invariant, &
            calculate_imposed_material_velocity_source, &
            calculate_imposed_material_velocity_absorption, &
            calculate_scalar_potential, calculate_projection_scalar_potential, &
            calculate_geostrophic_velocity, calculate_viscous_dissipation, & 
            calculate_adiabatic_heating_coefficient, calculate_adiabatic_heating_absorption, &
            calculate_viscous_dissipation_plus_surface_adiabat           
  
contains

  subroutine calculate_strain_rate(state, t_field)
    type(state_type), intent(inout) :: state
    type(tensor_field), intent(inout) :: t_field
    
    type(vector_field), pointer :: source_field
    type(vector_field), pointer :: positions

    ! Extract positions field:
    positions => extract_vector_field(state, "Coordinate")
    ! Extract source field:
    source_field => vector_source_field(state, t_field)
    ! Check that source field is not on a discontinuous mesh:
    call check_source_mesh_derivative(source_field, "strain_rate")
    ! Calculate strain_rate tensor:
    call strain_rate(source_field, positions, t_field)

  end subroutine calculate_strain_rate

  subroutine calculate_strain_rate_second_invariant(state, s_field)
    type(state_type), intent(inout) :: state
    type(scalar_field), intent(inout) :: s_field

    type(vector_field), pointer :: positions
    type(vector_field), pointer :: velocity

    type(tensor_field) :: strain_rate_tensor

    ewrite(1,*) 'In calculate_strain_rate_second_invariant'
    positions => extract_vector_field(state, "Coordinate")
    velocity  => extract_vector_field(state, "IteratedVelocity")

    ! Allocate strain_rate tensor:
    call allocate(strain_rate_tensor, s_field%mesh, name="strain_rate_II")

    call check_source_mesh_derivative(velocity, "strain_rate_second_invariant")

    ! Calculate strain_rate and second invariant:
    call strain_rate(velocity, positions, strain_rate_tensor)
    call tensor_second_invariant(strain_rate_tensor, s_field)

    ! Clean-up:
    call deallocate(strain_rate_tensor)

    ! Prin min and max:
    ewrite_minmax(s_field) 

  end subroutine calculate_strain_rate_second_invariant

  subroutine calculate_sediment_concentration_dependent_viscosity(state, t_field)
    ! calculates viscosity based upon total sediment concentration
    type(state_type), intent(inout) :: state
    type(tensor_field), intent(inout) :: t_field
    
    type(scalar_field_pointer), dimension(:), allocatable :: sediment_concs
    type(tensor_field), pointer :: zero_conc_viscosity
    type(scalar_field) :: rhs
    integer :: sediment_classes, i
    character(len = FIELD_NAME_LEN) :: field_name
    
    ewrite(1,*) 'In calculate_sediment_concentration_dependent_viscosity'

    sediment_classes = get_n_sediment_fields()

    if (sediment_classes > 0) then
        allocate(sediment_concs(sediment_classes))
        
        call get_sediment_item(state, 1, sediment_concs(1)%ptr)
        
        call allocate(rhs, sediment_concs(1)%ptr%mesh, name="Rhs")
        call set(rhs, 1.0)
        
        ! get sediment concentrations and remove c/0.65 from rhs
        do i=1, sediment_classes
           call get_sediment_item(state, i, sediment_concs(i)%ptr)
           call addto(rhs, sediment_concs(i)%ptr, scale=-(1.0/0.65))
        end do
        
        ! raise rhs to power of -1.625
        do i = 1, node_count(rhs)
           call set(rhs, i, node_val(rhs, i)**(-1.625))
        end do
        
        ! check for presence of ZeroSedimentConcentrationViscosity field
        if (.not. has_tensor_field(state, "ZeroSedimentConcentrationViscosity")) then
           FLExit("You must specify an zero sediment concentration viscosity to be able &
                &to calculate sediment concentration dependent viscosity field values")
        endif
        zero_conc_viscosity => extract_tensor_field(state, 'ZeroSedimentConcentrationViscosity')
        
        call set(t_field, zero_conc_viscosity)
        call scale(t_field, rhs)
        ewrite_minmax(t_field) 

        deallocate(sediment_concs)
        call deallocate(rhs)
    else
        ewrite(1,*) 'No sediment in problem definition'
    end if  
  end subroutine calculate_sediment_concentration_dependent_viscosity
  
  subroutine calculate_tensor_second_invariant(state, s_field)
    type(state_type), intent(inout) :: state
    type(scalar_field), intent(inout) :: s_field

    type(tensor_field), pointer :: source_field

    source_field => tensor_source_field(state, s_field)

    call tensor_second_invariant(source_field, s_field)

  end subroutine calculate_tensor_second_invariant

  subroutine calculate_viscous_dissipation(state, s_field)
    ! A routine to calculate the viscous dissipation. Currently
    ! assumes a constant viscosity tensor:
    type(state_type), intent(inout) :: state
    type(scalar_field), intent(inout) :: s_field

    type(vector_field), pointer :: positions
    type(vector_field), pointer :: velocity
    type(tensor_field), pointer :: viscosity

    type(scalar_field) :: velocity_divergence
    type(scalar_field) :: viscosity_component, viscosity_component_remap
    type(tensor_field) :: strain_rate_tensor

    integer :: dim1, dim2, node
    real :: val

    ewrite(1,*) 'In calculate_viscous_dissipation'

    ! Extract velocity field from state - will be used to calculate strain-
    ! rate tensor:
    velocity => extract_vector_field(state, "NonlinearVelocity")
    ! Check velocity field is not on a discontinous mesh:
    call check_source_mesh_derivative(velocity, "Viscous_Dissipation")

    ! Extract positions field from state:
    positions => extract_vector_field(state, "Coordinate")

    ! Allocate and initialize strain rate tensor:
    call allocate(strain_rate_tensor, s_field%mesh, "Strain_Rate_VD")
    call zero(strain_rate_tensor)

    ! Calculate strain rate tensor:
    call strain_rate(velocity, positions, strain_rate_tensor)

    ! Calculate velocity divergence for correct definition of stress:
    call allocate(velocity_divergence, s_field%mesh, 'Velocity_divergence')
    call div(velocity, positions, velocity_divergence)
    ewrite_minmax(velocity_divergence)

    ! Extract viscosity from state and remap to s_field mesh:
    viscosity => extract_tensor_field(state, "Viscosity")
    ! Extract first component of viscosity tensor from full tensor:
    !*** This is not ideal - only valid for constant viscosity tensors
    !*** though they can still vary spatially and temporally.
    viscosity_component = extract_scalar_field(viscosity,1,1)  
    call allocate(viscosity_component_remap, s_field%mesh, "RemappedViscosityComponent")
    call remap_field(viscosity_component, viscosity_component_remap)

    ! Calculate viscous dissipation (scalar s_field):
    do node=1,node_count(s_field)
       val = 0.
       do dim1 = 1, velocity%dim
          do dim2 = 1, velocity%dim
             if(dim1==dim2) then
                ! Add divergence of velocity term to diagonal only: 
                val = val + 2.*node_val(viscosity_component_remap, node) * & 
                     & (node_val(strain_rate_tensor,dim1,dim2,node)      - &
                     & 1./3. * node_val(velocity_divergence, node))**2
             else
                val = val + 2.*node_val(viscosity_component_remap, node) * & 
                     & node_val(strain_rate_tensor,dim1,dim2,node)**2   
             end if
          end do
       end do
       call set(s_field, node, val)
    end do

    ewrite_minmax(s_field)

    ! Deallocate:
    call deallocate(strain_rate_tensor)
    call deallocate(viscosity_component_remap)
    call deallocate(velocity_divergence)

  end subroutine calculate_viscous_dissipation

  subroutine calculate_viscous_dissipation_plus_surface_adiabat(state, s_field)
    ! A routine to calculate the viscous dissipation plus the surface adiabatic term. 
    ! Currently assumes a constant viscosity:
    type(state_type), intent(inout) :: state
    type(scalar_field), intent(inout) :: s_field

    type(scalar_field) :: surface_adiabat_field
    logical :: include_surface_adiabat
    
    ewrite(1,*) 'In calculate_viscous_dissipation_plus_surface_adiabat'

    ! This must be done first since it sets the scalar field to the viscous dissipation
    call calculate_viscous_dissipation(state, s_field)

    ! Allocate surface adiabat field, which will later be added to viscous dissipation:
    call allocate(surface_adiabat_field, s_field%mesh, "SurfaceAdiabatField")
    surface_adiabat_field%option_path = s_field%option_path

    ! Calculate surface adiabat term:
    call adiabatic_heating_coefficient(state, surface_adiabat_field, include_surface_adiabat=.true.)

    ! Add surface adiabat to viscous dissipation and clean up:
    call addto(s_field, surface_adiabat_field)
    call deallocate(surface_adiabat_field)

    ewrite_minmax(s_field)

  end subroutine calculate_viscous_dissipation_plus_surface_adiabat

  subroutine calculate_adiabatic_heating_coefficient(state, s_field)
    ! Calculates adiabatic heating coefficient, which is later multiplied
    ! by Temperature to form an absorption term :
    type(state_type), intent(inout) :: state
    type(scalar_field), intent(inout) :: s_field

    logical :: include_surface_adiabat

    ewrite(1,*) 'In calculate_adiabatic_heating_coefficient'

    call adiabatic_heating_coefficient(state, s_field, include_surface_adiabat=.false.)

    ewrite_minmax(s_field)

  end subroutine calculate_adiabatic_heating_coefficient

  subroutine adiabatic_heating_coefficient(state, s_field, include_surface_adiabat)
    ! Calculates adiabatic heating coefficient, which is later multiplied
    ! by Temperature to form an absorption term:
    type(state_type), intent(inout) :: state
    type(scalar_field), intent(inout) :: s_field
    logical :: include_surface_adiabat

    type(scalar_field), pointer :: thermal_expansion_local, density_local
    type(scalar_field), pointer :: reference_temperature_local
    type(vector_field), pointer :: velocity, gravity_direction

    type(vector_field) :: velocity_remap
    type(scalar_field) :: velocity_component, thermal_expansion_remap
    type(scalar_field) :: density_remap
    type(scalar_field) :: reference_temperature_remap
    real :: gravity_magnitude, gamma, rho0, T0
    integer :: node, stat_vel, stat_rho, stat_gamma, stat_reft

    character(len=OPTION_PATH_LEN) :: eos_option_path
    character(len=FIELD_NAME_LEN) :: density_name
    logical :: have_linear_eos, have_linearised_mantle_compressible_eos

    ewrite(1,*) 'In adiabatic_heating_coefficient'

    ! Get velocity field from state and verify that it can be remapped on to s_field%mesh
    ! in order to calculate inner product with gravity direction vector:
    velocity => extract_vector_field(state, "NonlinearVelocity")
    call allocate(velocity_remap, velocity%dim, s_field%mesh, "RemappedVelocity")
    call test_remap_validity(velocity,velocity_remap,stat_vel)
    call deallocate(velocity_remap)

    ! Extract gravitational info from state:
    gravity_direction => extract_vector_field(state, "GravityDirection")
    call get_option("/physical_parameters/gravity/magnitude", gravity_magnitude)

    if(stat_vel == 0) then

       call allocate(velocity_component, s_field%mesh, "VerticalVelocityComponent")
       ! Take inner product of velocity and gravity to determine vertical component of velocity:
       call inner_product(velocity_component, velocity, gravity_direction)
       
       ! Determine which EOS is relevant:
       eos_option_path='/material_phase::'//trim(state%name)//'/equation_of_state'
       have_linear_eos = (have_option(trim(eos_option_path)//'/fluids/linear'))
       have_linearised_mantle_compressible_eos = (have_option(trim(eos_option_path)//'/compressible/linearised_mantle'))

       if(have_linear_eos) then

          ! Get value for thermal expansion coefficient (constant)
          call get_option(trim(eos_option_path)//'/fluids/linear/temperature_dependency/thermal_expansion_coefficient', gamma)
          ! Get value for density:
          call get_option(trim(eos_option_path)//'/fluids/linear/reference_density', rho0)
          ! Calculate and set adiabatic heating coefficient:
          do node = 1, node_count(s_field)
             call set(s_field, node, -gamma*rho0*gravity_magnitude*node_val(velocity_component,node))
          end do

       elseif(have_linearised_mantle_compressible_eos) then

          ! Get spatially varying thermal expansion field and remap to s_field%mesh if possible and required:
          thermal_expansion_local=>extract_scalar_field(state,'IsobaricThermalExpansivity')       
          call allocate(thermal_expansion_remap, s_field%mesh, 'RemappedIsobaricThermalExpansivity')
          call test_remap_validity(thermal_expansion_local,thermal_expansion_remap,stat_gamma)
          if(stat_gamma == 0) then
             ! Remap to s_field%mesh:
             call remap_field(thermal_expansion_local, thermal_expansion_remap)
          else
             ! Remap is not possible:
             ewrite(-1,*) "IsobaricThermalExpansivity cannot be remapped to the adiabatic_heating_coefficient mesh,"
             ewrite(-1,*) "in the adiabatic heating coefficient algorithm."
             FLExit("Please put the IsobaricThermalExpansivity Field on a mesh that can be remapped to the adiabatic_heating_coefficient mesh.")
          end if

          ! Get CompressibleReferenceDensity and remap to s_field%mesh if possible and required:
          density_local => extract_scalar_field(state,'CompressibleReferenceDensity')
          call allocate(density_remap, s_field%mesh, 'RemappedCompressibleReferenceDensity')
          call test_remap_validity(density_local,density_remap,stat_rho)
          if(stat_rho == 0) then
             ! Remap to s_field%mesh:
             call remap_field(density_local, density_remap)
          else
             ! Remap is not possible:
             ewrite(-1,*) "CompressibleReferenceDensity cannot be remapped to the adiabatic_heating_coefficient mesh,"
             ewrite(-1,*) "in the adiabatic heating coefficient algorithm."
             FLExit("Please put the CompressibleReferenceDensity field on a mesh that can be remapped to the adiabatic_heating_coefficient mesh.")
          end if

          ! Calculate and set adiabatic heating coefficient:
          call get_option(trim(state%option_path)//'/scalar_field::Temperature/prognostic/equation[0]/density[0]/name', density_name)

          if (trim(density_name)=="CompressibleReferenceDensity") then

             do node = 1, node_count(s_field)
                gamma = node_val(thermal_expansion_remap,node)
                call set(s_field, node, -gamma * node_val(density_remap, node) &
                     * gravity_magnitude * node_val(velocity_component,node))
             end do

          elseif (trim(density_name)=="Density") then

             reference_temperature_local=>extract_scalar_field(state,'CompressibleReferenceTemperature')       
             call allocate(reference_temperature_remap, s_field%mesh, 'RemappedCompressibleReferenceTemperature')
             call test_remap_validity(reference_temperature_local,reference_temperature_remap,stat_reft)
             if(stat_reft == 0) then
                ! Remap to s_field%mesh:
                call remap_field(reference_temperature_local, reference_temperature_remap)
             else
                ! Remap is not possible:
                ewrite(-1,*) "Reference Temperature cannot be remapped to the adiabatic_heating_coefficient mesh,"
                ewrite(-1,*) "in the adiabatic heating coefficient algorithm."
                FLExit("Please put the CompressiblReferenceTemperature Field on a mesh that can be remapped to the adiabatic_heating_coefficient mesh.")
             end if

             do node = 1, node_count(s_field)
                gamma = node_val(thermal_expansion_remap,node)
                call set(s_field, node, ( (-gamma * node_val(density_remap, node) &
                     * gravity_magnitude * node_val(velocity_component,node) ) &
                     * (1. + gamma*node_val(reference_temperature_remap,node)) ) )
             end do

             call deallocate(reference_temperature_remap)

          else

             ewrite(-1,*) "Something has gone wrong in the adiabatic heating coefficient algorithm."
             FLExit("Density has an unknown name!")

          end if

          call deallocate(thermal_expansion_remap)
          call deallocate(density_remap)

       else

          FLExit("Selected EOS not yet configured for adiabatic_heating_algorithm")

       end if

       call deallocate(velocity_component)

    else

       ! Velocity cannot be remapped to the s_field%mesh. We therefore calculate the adiabatic heating coefficient
       ! using a 'projection':
       call adiabatic_heating_coefficient_projection(state, s_field)

    end if

    if(include_surface_adiabat) then
       ! Scale solutions by T0 and reverse sign:
       call get_option(trim(complete_field_path(trim(s_field%option_path))) // &
            "/algorithm[0]/surface_temperature", T0)
       call scale(s_field, -T0)
    end if

  end subroutine adiabatic_heating_coefficient

  subroutine adiabatic_heating_coefficient_projection(state, s_field)
    ! Calculates adiabatic heating coefficient, which is later multiplied
    ! by Temperature to form an absorption term. 
    type(state_type), intent(inout) :: state
    type(scalar_field), intent(inout) :: s_field

    type(vector_field), pointer :: positions, gravity_direction, velocity
    type(scalar_field), pointer :: thermal_expansion, density_local
    type(scalar_field), pointer :: dummyscalar, reference_temperature_local

    type(scalar_field) :: lumped_mass

    real :: gravity_magnitude, rho0, gamma
    integer :: ele

    character(len=OPTION_PATH_LEN) eos_option_path
    character(len=FIELD_NAME_LEN) :: density_name
    logical :: have_linear_eos, have_linearised_mantle_compressible_eos

    ewrite(1,*) 'In adiabatic_heating_coefficient_projection'

    ! Allocate dummy scalar for instances where a spatially varying reference density/thermal_expansion 
    ! field is not required:
    allocate(dummyscalar)
    call allocate(dummyscalar, s_field%mesh, name="DummyScalar", field_type=FIELD_TYPE_CONSTANT)
    call zero(dummyscalar)

    ! Extract gravitational unit vector and velocity from state:
    call get_option("/physical_parameters/gravity/magnitude", gravity_magnitude)
    gravity_direction => extract_vector_field(state, "GravityDirection")
    velocity => extract_vector_field(state, "NonlinearVelocity")
    
    ! Extract coordinates from state (required for detwei in integration):
    positions=>extract_vector_field(state,'Coordinate')

    ! Determine which EOS is relevant:
    eos_option_path='/material_phase::'//trim(state%name)//'/equation_of_state'
    have_linear_eos = (have_option(trim(eos_option_path)//'/fluids/linear'))
    have_linearised_mantle_compressible_eos = (have_option(trim(eos_option_path)//'/compressible/linearised_mantle'))

    ! Extract relevant parameters from options/state:
    if(have_linear_eos) then
       ! Get value for thermal expansion coefficient (constant)
       call get_option(trim(eos_option_path)//'/fluids/linear/temperature_dependency/thermal_expansion_coefficient', gamma)
       ! Get value for reference density (constant)
       call get_option(trim(eos_option_path)//'/fluids/linear/reference_density', rho0)
       ! As these values are constant for this eos, we need to point reference density and 
	   ! thermal_expansion to dummy scalars:
       density_local => dummyscalar
       thermal_expansion => dummyscalar
    elseif(have_linearised_mantle_compressible_eos) then
       ! Get spatially varying thermal expansion field:
       thermal_expansion=>extract_scalar_field(state,'IsobaricThermalExpansivity')
       ! Get CompressibleReferenceDensity field:
       density_local => extract_scalar_field(state,'CompressibleReferenceDensity')
    else
       FLExit("Selected EOS not yet configured for adiabatic_heating_coefficient_CV algorithm")
    endif

    ! Integrate to determine RHS:
    call zero(s_field)

    if(have_linear_eos) then

       do ele = 1, element_count(s_field)
          call integrate_RHS_ele(s_field, positions, velocity, gravity_direction, thermal_expansion, density_local, gravity_magnitude, ele, rho0=rho0, gamma=gamma)
       end do
       
    else if(have_linearised_mantle_compressible_eos) then             
       
       call get_option(trim(state%option_path)//'/scalar_field::Temperature/prognostic/equation[0]/density[0]/name', density_name)
       
       if (trim(density_name)=="CompressibleReferenceDensity") then
          
          do ele = 1, element_count(s_field)
             call integrate_RHS_ele(s_field, positions, velocity, gravity_direction, thermal_expansion, density_local, gravity_magnitude, ele)
          end do
          
       elseif (trim(density_name)=="Density") then
          
          reference_temperature_local=>extract_scalar_field(state,'CompressibleReferenceTemperature')       
          
          do ele = 1, element_count(s_field)
             call integrate_RHS_ele(s_field, positions, velocity, gravity_direction, thermal_expansion, density_local, gravity_magnitude, ele, reference_temperature=reference_temperature_local)
          end do
          
       else

          ewrite(-1,*) "Something has gone wrong in the adiabatic heating coefficient algorithm."
          FLExit("Density has an unknown name!")
          
       end if

    end if

    ! Compute inverse lumped mass matrix:
    call allocate(lumped_mass, s_field%mesh, name="Lumped_mass")        
    call compute_cv_mass(positions, lumped_mass)
    call invert(lumped_mass)
    
    ! Multiply RHS by inverse lumped mass matrix to determine final adiabatic_heating_coefficient:
    call scale(s_field, lumped_mass)
    
    ! Clean up:
    call deallocate(lumped_mass)
    call deallocate(dummyscalar)

  contains

    subroutine integrate_RHS_ele(s_field, positions, velocity, gravity_direction, thermal_expansion, density_local, gravity_magnitude, ele, reference_temperature, rho0, gamma)
      type(scalar_field), intent(inout) :: s_field
      type(vector_field), intent(in), pointer :: positions
      type(scalar_field), intent(in), pointer :: thermal_expansion, density_local
      type(vector_field), intent(in), pointer :: velocity, gravity_direction
      real, intent(in)    :: gravity_magnitude
      type(scalar_field), optional, intent(in), pointer :: reference_temperature
      real, optional      :: rho0, gamma
      integer, intent(in) :: ele
      
      ! For integration:
      integer :: dim, gi
      real, dimension(ele_loc(s_field,ele)) :: ele_val
      real, dimension(positions%dim, ele_ngi(positions, ele)) :: positions_quad, velocity_quad, gravity_direction_quad
      real, dimension(ele_ngi(positions, ele)) :: detwei, density_quad, gamma_quad, inner_prod, reference_temperature_quad

      ! Evaluate key parameters at gauss points:
      call transform_to_physical(positions, ele, detwei = detwei)
      positions_quad              = ele_val_at_quad(positions, ele)
      velocity_quad               = ele_val_at_quad(velocity, ele)
      gravity_direction_quad      = ele_val_at_quad(gravity_direction, ele)
      density_quad                = ele_val_at_quad(density_local, ele)
      gamma_quad                  = ele_val_at_quad(thermal_expansion, ele)

      ! Calculate inner product of velocity and gravity unit vector at gauss points:
      inner_prod = 0.
      do dim = 1, velocity%dim
         do gi = 1, ele_ngi(velocity, ele)
            inner_prod(gi) = inner_prod(gi) + velocity_quad(dim,gi)*gravity_direction_quad(dim,gi)
         end do
      end do

      ! Evaluate nodal values of integral:
      if(present(rho0) .and. present(gamma)) then ! Linear EOS
         ele_val = shape_rhs(ele_shape(s_field,ele), -inner_prod*gravity_magnitude*gamma*rho0*detwei)
      elseif(present(reference_temperature)) then ! Density Coefficient
         reference_temperature_quad  = ele_val_at_quad(reference_temperature, ele)
         ele_val = shape_rhs(ele_shape(s_field,ele), (-inner_prod*gravity_magnitude*gamma_quad*density_quad*detwei) * (1.+gamma_quad*reference_temperature_quad))
      else ! CompressibleReferenceDensity Coefficient
         ele_val = shape_rhs(ele_shape(s_field,ele), -inner_prod*gravity_magnitude*gamma_quad*density_quad*detwei)
      end if

      ! Add this to the global RHS vector:
      call addto(s_field, ele_nodes(s_field, ele), ele_val)
           
    end subroutine integrate_RHS_ele

  end subroutine adiabatic_heating_coefficient_projection

  subroutine calculate_adiabatic_heating_absorption(state, s_field)
    ! Calculates adiabatic heating absorption term - i.e. the adiabatic heating
    ! coefficient, multiplied by T:
    type(state_type), intent(inout) :: state
    type(scalar_field), intent(inout) :: s_field

    type(scalar_field), pointer :: temperature
    type(scalar_field) :: temperature_remap

    ewrite(1,*) 'In calculate_adiabatic_heating_absorption'

    call calculate_adiabatic_heating_coefficient(state,s_field)

    ! Extract temperature from state:
    temperature => extract_scalar_field(state, "Temperature")
    call allocate(temperature_remap, s_field%mesh, "TemperatureRemap")
    call remap_field(temperature, temperature_remap)

    ! Multiply adiabatic coefficient by Temperature to attain full absorption term:
    call scale(s_field, temperature_remap)

    call deallocate(temperature_remap)

  end subroutine calculate_adiabatic_heating_absorption

  subroutine calculate_bulk_viscosity(states, t_field)
    type(state_type), dimension(:), intent(inout) :: states
    type(tensor_field), intent(inout) :: t_field

    character(len = OPTION_PATH_LEN) :: mean_type
    
    call get_option(trim(complete_field_path(trim(t_field%option_path))) // &
                    "/algorithm[0]/mean/name", mean_type, default="arithmetic")

    call calculate_bulk_property(states, t_field, "MaterialViscosity", &
      & mean_type = mean_type, momentum_diagnostic = .true.)
  
  end subroutine calculate_bulk_viscosity
  
  subroutine calculate_imposed_material_velocity_source(states, state_index, v_field)
    type(state_type), dimension(:), intent(inout) :: states
    integer, intent(in) :: state_index
    type(vector_field), intent(inout) :: v_field
    
    logical :: prescribed
    integer :: i, stat
    type(vector_field), pointer :: absorption, mat_vel
  
    call zero(v_field)
    
    do i = 1, size(states)
  
      mat_vel => extract_vector_field(states(i), "MaterialVelocity", stat)
      
      if(stat==0) then
      
        call add_scaled_material_property(states(i), v_field, mat_vel, &
                                          momentum_diagnostic=.true.)
                                          
      else
        ! alternatively use the Velocity field from the state
        
        mat_vel => extract_vector_field(states(i), "Velocity", stat)
        
        if(stat==0) then
          prescribed = have_option(trim(mat_vel%option_path)//"/prescribed")
          
          if(prescribed.and.(.not.aliased(mat_vel))) then
            ! but make sure it's prescribed and not aliased
          
            call add_scaled_material_property(states(i), v_field, mat_vel, &
                                              momentum_diagnostic=.true.)
                                              
          end if
        
        end if
        
      end if

    end do

    absorption => extract_vector_field(states(state_index), "VelocityAbsorption")
    call scale(v_field, absorption)
  
  end subroutine calculate_imposed_material_velocity_source

  subroutine calculate_imposed_material_velocity_absorption(states, v_field)
    type(state_type), dimension(:), intent(inout) :: states
    type(vector_field), intent(inout) :: v_field
    
    logical :: prescribed
    integer :: i, stat
    real :: dt
    real, dimension(v_field%dim) :: factor
    type(vector_field) :: temp_abs
    type(vector_field), pointer :: mat_vel
    
    call get_option("/timestepping/timestep", dt)
    call get_option(trim(complete_field_path(trim(v_field%option_path))) // &
                    "/algorithm[0]/relaxation_factor", factor, default=spread(1.0, 1, v_field%dim))
    
    call allocate(temp_abs, v_field%dim, v_field%mesh, "TemporaryAbsorption", &
                  field_type=FIELD_TYPE_CONSTANT)
    call set(temp_abs, factor/dt)
        
    call zero(v_field)
    
    do i = 1, size(states)
  
      mat_vel => extract_vector_field(states(i), "MaterialVelocity", stat)
      
      if(stat==0) then
      
        call add_scaled_material_property(states(i), v_field, temp_abs, &
                                          momentum_diagnostic=.true.)

      else
        ! alternatively use the Velocity field from the state
        
        mat_vel => extract_vector_field(states(i), "Velocity", stat)
        
        if(stat==0) then
          prescribed = have_option(trim(mat_vel%option_path)//"/prescribed")
          
          if(prescribed.and.(.not.aliased(mat_vel))) then
            ! but make sure it's prescribed and not aliased
          
            call add_scaled_material_property(states(i), v_field, temp_abs, &
                                              momentum_diagnostic=.true.)
            
          end if
        
        end if
        
      end if

    end do
    
    call deallocate(temp_abs)
    
  end subroutine calculate_imposed_material_velocity_absorption

  subroutine calculate_buoyancy(state, v_field)
    type(state_type), intent(in) :: state
    type(vector_field), intent(inout) :: v_field
    
    integer :: i, stat
    real :: gravity_magnitude
    type(scalar_field), pointer :: buoyancy_density
    type(vector_field), pointer :: gravity
  
    ewrite(1, *) "In calculate_buoyancy"
    
    buoyancy_density => extract_scalar_field(state, "VelocityBuoyancyDensity", stat = stat)
    if(stat /= 0) then
      ewrite(0, *) "Warning: Cannot calculate Buoyancy without VelocityBuoyancyDensity field"
      call zero(v_field)
      ewrite(1, *) "Exiting calculate_buoyancy"
      return
    end if    
    ewrite_minmax(buoyancy_density)
    
    gravity => extract_vector_field(state, "GravityDirection", stat = stat)
    if(stat /= 0) then
      ewrite(0, *) "Warning: Cannot calculate Buoyancy without GravityDirection field"
      call zero(v_field)
      ewrite(1, *) "Exiting calculate_buoyancy"
      return
    end if    
    ewrite_minmax(gravity)
    
    call get_option("/physical_parameters/gravity/magnitude", gravity_magnitude)
    ewrite(2, *) "Gravity magnitude = ", gravity_magnitude
    
    if(.not. v_field%mesh == buoyancy_density%mesh) then
      ewrite(-1, *) "VelocityBuoyancyDensity mesh: " // trim(buoyancy_density%mesh%name)
      FLExit("Buoyancy must be on the VelocityBuoyancyDensity mesh")
    end if
    
    do i = 1, node_count(v_field)
      call set(v_field, i, node_val(gravity, i) * node_val(buoyancy_density, i) * gravity_magnitude)
    end do
    
    ewrite(1, *) "Exiting calculate_buoyancy"
  
  end subroutine calculate_buoyancy
  
  subroutine calculate_coriolis(state, v_field)
    type(state_type), intent(inout) :: state
    type(vector_field), intent(inout) :: v_field
   
    character(len = OPTION_PATH_LEN) :: base_path
    
    base_path = trim(complete_field_path(v_field%option_path)) // "/algorithm"

    if(have_option(trim(base_path) // "/consistent_interpolation")) then
      call compute_coriolis_ci(state, v_field)
    else if(have_option(trim(base_path) // "/galerkin_projection")) then
      if(have_option(trim(base_path) // "/galerkin_projection/lump_mass")) then
        call compute_coriolis_gp_lumped(state, v_field)
      else
        call compute_coriolis_gp(state, v_field, option_path = trim(base_path) // "/galerkin_projection")
      end if
    else
      FLAbort("Failed to determine interpolation method")
    end if  
      
  end subroutine calculate_coriolis
  
  subroutine compute_coriolis_ci(state, coriolis)
    type(state_type), intent(in) :: state
    type(vector_field), intent(inout) :: coriolis
  
    integer :: i
    type(vector_field) :: positions, velocity_remap
    type(vector_field), pointer :: velocity
    
    positions = get_nodal_coordinate_field(state, coriolis%mesh)
    velocity => extract_vector_field(state, "Velocity")
    
    if(velocity%mesh == coriolis%mesh) then
      velocity_remap = velocity
      call incref(velocity_remap)
    else
      call allocate(velocity_remap, velocity%dim, coriolis%mesh, "VelocityRemap")
      call remap_field(velocity, velocity_remap)
    end if
    
    do i = 1, node_count(coriolis)
      call set(coriolis, i, coriolis_val(node_val(positions, i), node_val(velocity_remap, i)))
    end do
    
    call deallocate(positions)
    call deallocate(velocity_remap)
    
  end subroutine compute_coriolis_ci
  
  subroutine compute_coriolis_gp(state, coriolis, option_path)
    type(state_type), intent(inout) :: state
    type(vector_field), intent(inout) :: coriolis
    character(len = *), optional, intent(in) :: option_path
  
    integer :: i
    type(csr_matrix), pointer :: mass
    type(vector_field) :: rhs
    type(vector_field), pointer :: positions, velocity
    
    positions => extract_vector_field(state, "Coordinate")
    velocity => extract_vector_field(state, "Velocity")
    
    mass => get_mass_matrix(state, coriolis%mesh)
    call allocate(rhs, coriolis%dim, coriolis%mesh, "CoriolisRhs")

    call zero(rhs)
    do i = 1, ele_count(rhs)
      call assemble_coriolis_ele(i, positions, velocity, rhs)
    end do

    call petsc_solve(coriolis, mass, rhs, option_path = option_path)

    call deallocate(rhs)
  
  end subroutine compute_coriolis_gp
  
  subroutine compute_coriolis_gp_lumped(state, coriolis)
    type(state_type), intent(inout) :: state
    type(vector_field), intent(inout) :: coriolis
    
    integer :: i
    type(scalar_field), pointer :: masslump
    type(vector_field), pointer :: positions, velocity
    
    positions => extract_vector_field(state, "Coordinate")
    velocity => extract_vector_field(state, "Velocity")
    
    masslump => get_lumped_mass(state, coriolis%mesh)

    call zero(coriolis)
    do i = 1, ele_count(coriolis)
      call assemble_coriolis_ele(i, positions, velocity, coriolis)
    end do
    
    do i = 1, coriolis%dim
      coriolis%val(i,:) = coriolis%val(i,:) / masslump%val
    end do
    
  end subroutine compute_coriolis_gp_lumped
    
  subroutine assemble_coriolis_ele(ele, positions, velocity, rhs)
    integer, intent(in) :: ele
    type(vector_field), intent(in) :: positions
    type(vector_field), intent(in) :: velocity
    type(vector_field), intent(inout) :: rhs

    real, dimension(ele_ngi(rhs, ele)) :: detwei

    call transform_to_physical(positions, ele, detwei = detwei)

    call addto(rhs, ele_nodes(rhs, ele), &
      & shape_vector_rhs(ele_shape(rhs, ele), &
        & coriolis_val(ele_val_at_quad(positions, ele), ele_val_at_quad(velocity, ele)), &
      & detwei))

  end subroutine assemble_coriolis_ele
  
  subroutine calculate_scalar_potential(state, s_field)
    type(state_type), intent(inout) :: state
    type(scalar_field), intent(inout) :: s_field

    type(vector_field), pointer :: source_field

    source_field => vector_source_field(state, s_field)
    call geopressure_decomposition(state, source_field, s_field, &
      & option_path = trim(complete_field_path(s_field%option_path)) // "/algorithm")

  end subroutine calculate_scalar_potential
  
  subroutine calculate_projection_scalar_potential(state, s_field)
    type(state_type), intent(inout) :: state
    type(scalar_field), intent(inout) :: s_field

    character(len = OPTION_PATH_LEN) :: bcfield_name, path
    type(scalar_field), pointer :: gp
    type(vector_field), pointer :: bcfield, source_field

    source_field => vector_source_field(state, s_field, index = 1)
    
    path = trim(complete_field_path(s_field%option_path)) // "/algorithm"

    if(have_option(trim(path) // "/bc_field")) then
      call get_option(trim(path) // "/bc_field/name", bcfield_name)
      bcfield => extract_vector_field(state, bcfield_name)
    else
      bcfield => source_field
    end if
    
    if(have_option(trim(path) // "/source_field_2_name")) then
      gp => scalar_source_field(state, s_field, index = 2)
      call projection_decomposition(state, source_field, s_field, &
        & bcfield = bcfield, gp = gp, option_path = path)
    else
      call projection_decomposition(state, source_field, s_field, &
        & bcfield = bcfield, option_path = path)
    end if

  end subroutine calculate_projection_scalar_potential

  subroutine calculate_geostrophic_velocity(state, v_field)
    type(state_type), intent(inout) :: state
    type(vector_field), intent(inout) :: v_field
    
    character(len = OPTION_PATH_LEN) :: path
    integer :: stat
    real :: scale_factor
    type(scalar_field), pointer :: source_field
    type(cmc_matrices) :: matrices
    type(vector_field), pointer :: velocity
    
    source_field => scalar_source_field(state, v_field)
    velocity => extract_vector_field(state, "Velocity")
    path = trim(complete_field_path(v_field%option_path)) // "/algorithm"
    call allocate(matrices, state, velocity, source_field, option_path = path, add_cmc = .false.)
    
    call geostrophic_velocity(matrices, state, v_field, source_field) 
    
    call deallocate(matrices)
    
    call get_option(trim(path) // "/scale_factor", scale_factor, stat = stat)
    if(stat == SPUD_NO_ERROR) call scale(v_field, scale_factor)
  
  end subroutine calculate_geostrophic_velocity

end module momentum_diagnostics
