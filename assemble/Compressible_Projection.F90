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

module compressible_projection
  use fldebug
  use state_module
  use sparse_tools
  use spud
  use fields
  use sparse_matrices_fields
  use field_options
  use equation_of_state, only: compressible_eos, compressible_material_eos
  use global_parameters, only: OPTION_PATH_LEN
  use fefields, only: compute_lumped_mass
  use state_fields_module
  use upwind_stabilisation
  implicit none 

  ! Buffer for output messages.
  character(len=255), private :: message

  private
  public :: assemble_compressible_projection_cv, assemble_compressible_projection_cg, &
            update_compressible_density, include_implicit_pressure_buoyancy, compressible_projection_check_options

  ! Stabilisation schemes
  integer, parameter :: STABILISATION_NONE = 0, &
    & STABILISATION_STREAMLINE_UPWIND = 1, STABILISATION_SUPG = 2
  ! Stabilisation scheme
  integer :: stabilisation_scheme
  integer :: nu_bar_scheme
  real :: nu_bar_scale

contains

  subroutine assemble_compressible_projection_cv(state, cmc, dt, theta_pg, theta_divergence, cmcget, rhs)

    ! inputs:
    ! bucket full of fields
    type(state_type), dimension(:), intent(inout) :: state

    type(csr_matrix), intent(inout) :: cmc

    real, intent(in) :: dt
    real, intent(in) :: theta_pg, theta_divergence
    logical, intent(in) :: cmcget
    
    ! not assembled if just assembling the auxilliary schur matrix
    type(scalar_field), intent(inout), optional :: rhs

    if((size(state)==1).and.(.not.has_scalar_field(state(1), "MaterialVolumeFraction"))) then
    
      call assemble_1mat_compressible_projection_cv(state(1), cmc, dt, &
                                                    theta_pg, theta_divergence, cmcget, rhs=rhs)
      
    else
    
      call assemble_mmat_compressible_projection_cv(state, cmc, dt, cmcget, rhs=rhs)
      
    end if
    

  end subroutine assemble_compressible_projection_cv

  subroutine assemble_1mat_compressible_projection_cv(state, cmc, dt, &
                                                      theta_pg, theta_divergence, cmcget, rhs)

    ! inputs:
    ! bucket full of fields
    type(state_type), intent(inout) :: state

    type(csr_matrix), intent(inout) :: cmc

    real, intent(in) :: dt
    real, intent(in) :: theta_pg, theta_divergence
    logical, intent(in) :: cmcget

    type(scalar_field), intent(inout), optional :: rhs
    
    ! local:
    type(scalar_field) :: eospressure, drhodp
    type(scalar_field), pointer :: density, olddensity
    type(scalar_field), pointer :: pressure
    type(scalar_field), pointer :: p_lumpedmass
    type(scalar_field) :: lhsfield, absrhs
    
    type(scalar_field), pointer :: source, absorption
    integer :: stat

    real :: atmospheric_pressure, theta
    
    logical :: assemble_rhs, exclude_mass

    ewrite(1,*) 'Entering assemble_1mat_compressible_projection_cv'
    
    assemble_rhs = present(rhs)
    if(assemble_rhs) call zero(rhs)

    pressure=>extract_scalar_field(state, "Pressure")
    call get_option(trim(pressure%option_path)//'/prognostic/atmospheric_pressure', &
                    atmospheric_pressure, default=0.0)
    
    if(pressure%mesh%shape%degree>1) then
      ! try lumping on the submesh
      p_lumpedmass => get_lumped_mass_on_submesh(state, pressure%mesh)
    else
      ! find the lumped mass
      p_lumpedmass => get_lumped_mass(state, pressure%mesh)
    end if
    ewrite_minmax(p_lumpedmass)
          
    if(cmcget) call allocate(lhsfield, pressure%mesh, "LHSField")

    call allocate(eospressure, pressure%mesh, 'EOSPressure')
    call allocate(drhodp, pressure%mesh, 'DerivativeDensityWRTBulkPressure')

    call zero(eospressure)
    call zero(drhodp)

    call compressible_eos(state, pressure=eospressure, drhodp=drhodp)

    density=>extract_scalar_field(state,'Density')
    ewrite_minmax(density)
    olddensity=>extract_scalar_field(state,'OldDensity')
    ewrite_minmax(olddensity)

    exclude_mass = have_option(trim(density%option_path)//&
                                "/prognostic/spatial_discretisation/control_volumes/mass_terms/exclude_mass_terms")

    call get_option(trim(density%option_path)//"/prognostic/temporal_discretisation/theta", theta)

    if(cmcget) then
      call set(lhsfield, p_lumpedmass)
      call scale(lhsfield, drhodp)
      if(.not.exclude_mass) then
        call addto_diag(cmc, lhsfield, scale=1./(dt*dt*theta_divergence*theta_pg))
      end if
    end if
    
    if(assemble_rhs) then
    
      !     rhs = p_lumpedmass* &
      !      ( (1./dt)*(olddensity - density + drhodp*(eospressure - (pressure + atmospheric_pressure)))
      !       +(absorption)*(drhodp*theta_pg*(eospressure - (pressure + atmospheric_pressure)) - theta_pg*density - (1-theta_pg)*olddensity)
      !       +source)
      if(exclude_mass) then
        call zero(rhs)
      else
        call set(rhs, pressure)
        call addto(rhs, atmospheric_pressure)
        call scale(rhs, -1.0)
        call addto(rhs, eospressure)
        call scale(rhs, drhodp)
        call addto(rhs, density, -1.0)
        call addto(rhs, olddensity)
        call scale(rhs, (1./dt))
      end if
      
      source => extract_scalar_field(state, "DensitySource", stat=stat)
      if(stat==0) then
        call addto(rhs, source)
      end if
      
    end if
    
    absorption => extract_scalar_field(state, "DensityAbsorption", stat=stat)
    if(stat==0) then
      if(assemble_rhs) then
        call allocate(absrhs, absorption%mesh, "AbsorptionRHS")
        
        call set(absrhs, pressure)
        call addto(absrhs, atmospheric_pressure)
        call scale(absrhs, -1.0)
        call addto(absrhs, eospressure)
        call scale(absrhs, drhodp)
        call scale(absrhs, theta)
        call addto(absrhs, density, -theta)
        call addto(absrhs, olddensity, -(1-theta))
        call scale(absrhs, absorption)
        
        call addto(rhs, absrhs)
        
        call deallocate(absrhs)
      end if
      
      if(cmcget) then
        call scale(lhsfield, absorption)
        call addto_diag(cmc, lhsfield, scale=(theta/(dt*theta_divergence*theta_pg)))
      end if
    end if
    
    if(assemble_rhs) call scale(rhs, p_lumpedmass)
    
    call deallocate(eospressure)
    call deallocate(drhodp)

    if(cmcget) call deallocate(lhsfield)

  end subroutine assemble_1mat_compressible_projection_cv

  subroutine assemble_mmat_compressible_projection_cv(state, cmc, dt, cmcget, rhs)

    ! inputs:
    ! bucket full of fields
    type(state_type), dimension(:), intent(inout) :: state

    type(csr_matrix), intent(inout) :: cmc

    real, intent(in) :: dt
    logical, intent(in) :: cmcget

    type(scalar_field), intent(inout), optional :: rhs
    
    ! local:
    integer :: i, stat
    character(len=OPTION_PATH_LEN) :: pressure_option_path

    type(scalar_field) :: materialpressure, materialdrhodp, density, &
                          olddensity, matdrhodpp, drhodp
    type(scalar_field), pointer :: volumefraction, oldvolumefraction, materialdensity, oldmaterialdensity
    type(scalar_field), pointer :: dummy_ones

    type(scalar_field), pointer :: pressure
    type(vector_field), pointer :: positions
    type(scalar_field) :: lumped_mass, tempfield
    
    logical :: compressible_eos, assemble_rhs

    real :: atmospheric_pressure

    ewrite(1,*) 'Entering assemble_mmat_compressible_projection_cv'

    pressure=>extract_prognostic_pressure(state, stat=stat)
    if(stat/=0) then
       ! how did we end up here?
       FLAbort("In assemble_mmat_compressible_projection_cv without a pressure")
    end if
    pressure_option_path=trim(pressure%option_path)
    
    compressible_eos = .false.
    state_loop: do i = 1, size(state)
      compressible_eos = have_option("/material_phase::"//trim(state(i)%name)//"/equation_of_state/compressible")
      if(compressible_eos) then
        exit state_loop
      end if
    end do state_loop
    
    assemble_rhs = present(rhs)
    if(assemble_rhs) call zero(rhs)
   
    if (compressible_eos) then
    
      positions=>extract_vector_field(state(1), "Coordinate")
      call allocate(lumped_mass, pressure%mesh, "LumpedMassField")
      if(cmcget) call allocate(tempfield, pressure%mesh, "TemporaryAssemblyField")
      call compute_lumped_mass(positions, lumped_mass)

      allocate(dummy_ones)
      call allocate(dummy_ones, pressure%mesh, "DummyOnesField")
      call set(dummy_ones, 1.0)

      call get_option(trim(pressure_option_path)//'/prognostic/atmospheric_pressure', &
                      atmospheric_pressure, default=0.0)

      call allocate(materialpressure, pressure%mesh, 'MaterialEOSPressure')
      call allocate(materialdrhodp, pressure%mesh, 'DerivativeMaterialdensityWRTBulkPressure')

      call allocate(density, pressure%mesh, 'MaterialDensity')
      call allocate(olddensity, pressure%mesh, 'OldMaterialDensity')
      call allocate(matdrhodpp, pressure%mesh, 'MaterialPressure')
      call allocate(drhodp, pressure%mesh, 'Drhodp')

      density%val = 0.0
      olddensity%val = 0.0
      matdrhodpp%val = 0.0
      drhodp%val=0.0

      do i = 1,size(state)

        materialpressure%val=0.0
        materialdrhodp%val=0.0

        call compressible_material_eos(state(i), materialpressure=materialpressure, materialdrhodp=materialdrhodp)

        volumefraction=>extract_scalar_field(state(i),'MaterialVolumeFraction', stat=stat)
        if(stat==0) then
          oldvolumefraction=>extract_scalar_field(state(i),'OldMaterialVolumeFraction')
          materialdensity=>extract_scalar_field(state(i),'MaterialDensity')
          oldmaterialdensity=>extract_scalar_field(state(i),'OldMaterialDensity')

          density%val = density%val &
                            + materialdensity%val*volumefraction%val
          olddensity%val = olddensity%val &
                              + oldmaterialdensity%val*oldvolumefraction%val
          matdrhodpp%val = matdrhodpp%val &
                                + materialpressure%val*materialdrhodp%val*volumefraction%val
          drhodp%val = drhodp%val &
                            + materialdrhodp%val*volumefraction%val
        endif

      end do

      if(cmcget) then
        call zero(tempfield)
        tempfield%val = (1./(dt*dt))*lumped_mass%val*drhodp%val
        call addto_diag(cmc, tempfield)
      end if

      if(assemble_rhs) then
        rhs%val = (1./dt)*lumped_mass%val* &
                          ( &
                            olddensity%val &
                          - density%val &
                          ) &
              +(1./dt)*lumped_mass%val* &
                          ( &
                            matdrhodpp%val &
                          - drhodp%val*(pressure%val+atmospheric_pressure) &
                          )
      end if

      call deallocate(density)
      call deallocate(olddensity)
      call deallocate(matdrhodpp)
      call deallocate(drhodp)

      call deallocate(materialpressure)
      call deallocate(materialdrhodp)

      call deallocate(lumped_mass)
      if(cmcget) call deallocate(tempfield)
      call deallocate(dummy_ones)
      deallocate(dummy_ones)

    end if

  end subroutine assemble_mmat_compressible_projection_cv
  
  subroutine assemble_compressible_projection_cg(state, cmc, dt, theta_pg, theta_divergence, cmcget, rhs)

    ! inputs:
    ! bucket full of fields
    type(state_type), dimension(:), intent(inout) :: state

    type(csr_matrix), intent(inout) :: cmc

    real, intent(in) :: dt
    real, intent(in) :: theta_pg, theta_divergence
    logical, intent(in) :: cmcget
    
    ! this isn't required to be assembled when we're just assembling the auxilliary schur complement matrix
    type(scalar_field), intent(inout), optional :: rhs

    if((size(state)==1).and.(.not.has_scalar_field(state(1), "MaterialVolumeFraction"))) then
    
      call assemble_1mat_compressible_projection_cg(state(1), cmc, dt, &
                                                    theta_pg, theta_divergence, cmcget, rhs=rhs)
      
    else
      
        FLExit("Multimaterial compressible continuous_galerkin pressure not possible.")
      
    end if
    

  end subroutine assemble_compressible_projection_cg

  subroutine assemble_1mat_compressible_projection_cg(state, cmc, dt, &
                                                      theta_pg, theta_divergence, cmcget, rhs)

    ! inputs:
    ! bucket full of fields
    type(state_type), intent(inout) :: state

    type(csr_matrix), intent(inout) :: cmc

    real, intent(in) :: dt
    real, intent(in) :: theta_pg, theta_divergence
    logical, intent(in) :: cmcget

    type(scalar_field), intent(inout), optional :: rhs
    
    ! local
    type(mesh_type), pointer :: test_mesh

    type(vector_field), pointer :: field

    integer, dimension(:), pointer :: test_nodes

    real, dimension(:), allocatable :: ele_rhs
    type(element_type), pointer :: test_shape_ptr
    type(element_type) :: test_shape
    real, dimension(:,:,:), allocatable :: dtest_t
    real, dimension(:), allocatable :: detwei
    real, dimension(:,:,:), allocatable :: j_mat
    
    real, dimension(:), allocatable :: density_at_quad, olddensity_at_quad, p_at_quad, &
                                      drhodp_at_quad, eosp_at_quad, abs_at_quad
    real, dimension(:,:), allocatable :: nlvelocity_at_quad

    ! loop integers
    integer :: ele

    ! pointer to coordinates
    type(vector_field), pointer :: coordinate, nonlinearvelocity, velocity
    type(scalar_field), pointer :: pressure, density, olddensity
    type(scalar_field), pointer :: source, absorption
    type(scalar_field) :: eospressure, drhodp
    real :: theta, atmospheric_pressure

    real, dimension(:,:), allocatable :: ele_mat
    
    logical :: have_absorption, have_source, exclude_mass, assemble_rhs
    integer :: stat

    ! =============================================================
    ! Subroutine to construct the matrix CT_m (a.k.a. C1/2/3T).
    ! =============================================================

    ewrite(1,*) 'Entering assemble_1mat_compressible_projection_cg'
    
    assemble_rhs = present(rhs)
    if(assemble_rhs) call zero(rhs)

    ! only do all this if we need to make cmc (otherwise we'd be adding repeatedly)
    coordinate=> extract_vector_field(state, "Coordinate")
    
    density => extract_scalar_field(state, "Density")
    olddensity => extract_scalar_field(state, "OldDensity")
    
    absorption => extract_scalar_field(state, "DensityAbsorption", stat=stat)
    have_absorption = (stat==0)
    if(have_absorption) then
      ewrite(2,*) 'Have DensityAbsorption'
    end if
    
    source => extract_scalar_field(state, "DensitySource", stat=stat)
    have_source = (stat==0)
    if(have_source) then
      ewrite(2,*) 'Have DensitySource'
    end if
    
    exclude_mass = have_option(trim(density%option_path)//&
                                "/prognostic/spatial_discretisation/continuous_galerkin/mass_terms/exclude_mass_terms")

    velocity=>extract_vector_field(state, "Velocity")
    nonlinearvelocity=>extract_vector_field(state, "NonlinearVelocity") ! maybe this should be updated after the velocity solve?
    
    pressure => extract_scalar_field(state, "Pressure")

    call get_option(trim(pressure%option_path)//'/prognostic/atmospheric_pressure', &
                    atmospheric_pressure, default=0.0)

    ! these are put on the density mesh, which should be of sufficient order to represent
    ! the multiplication of the eos (of course that may not be possible in which case
    ! something should be done at the gauss points instead)
    call allocate(eospressure, density%mesh, 'EOSPressure')
    call allocate(drhodp, density%mesh, 'DerivativeDensityWRTBulkPressure')

    call zero(eospressure)
    call zero(drhodp)

    ! this needs to be changed to be evaluated at the quadrature points!
    call compressible_eos(state, pressure=eospressure, drhodp=drhodp)

    ewrite_minmax(density)
    ewrite_minmax(olddensity)

    if(have_option(trim(density%option_path) // &
                        "/prognostic/spatial_discretisation/continuous_galerkin/&
                        &stabilisation/streamline_upwind_petrov_galerkin")) then
      ewrite(2, *) "SUPG stabilisation"
      stabilisation_scheme = STABILISATION_SUPG
      call get_upwind_options(trim(density%option_path) // & 
                              "/prognostic/spatial_discretisation/continuous_galerkin/&
                              &stabilisation/streamline_upwind_petrov_galerkin", &
                              & nu_bar_scheme, nu_bar_scale)
    else
      ewrite(2, *) "No stabilisation"
      stabilisation_scheme = STABILISATION_NONE
    end if

    call get_option(trim(density%option_path)//"/prognostic/temporal_discretisation/theta", theta)
    
    test_mesh => pressure%mesh
    field => velocity
    
    allocate(dtest_t(ele_loc(test_mesh, 1), ele_ngi(test_mesh, 1), field%dim), &
            detwei(ele_ngi(field, 1)), &
            density_at_quad(ele_ngi(density, 1)), &
            olddensity_at_quad(ele_ngi(density, 1)), &
            nlvelocity_at_quad(field%dim, ele_ngi(field, 1)), &
            j_mat(field%dim, field%dim, ele_ngi(density, 1)), &
            drhodp_at_quad(ele_ngi(drhodp, 1)), &
            eosp_at_quad(ele_ngi(eospressure, 1)), &
            abs_at_quad(ele_ngi(density, 1)), &
            p_at_quad(ele_ngi(pressure, 1)))
    
    if(cmcget) then
      allocate(ele_mat(ele_loc(test_mesh, 1), ele_loc(test_mesh, 1)))
    end if
    if(assemble_rhs) then
      allocate(ele_rhs(ele_loc(test_mesh, 1)))
    end if
    
    do ele=1, element_count(test_mesh)
    
      test_nodes=>ele_nodes(test_mesh, ele)

      test_shape_ptr => ele_shape(test_mesh, ele)
      
      density_at_quad = ele_val_at_quad(density, ele)
      olddensity_at_quad = ele_val_at_quad(olddensity, ele)
      
      p_at_quad = ele_val_at_quad(pressure, ele) + atmospheric_pressure
                        
      nlvelocity_at_quad = ele_val_at_quad(nonlinearvelocity, ele)
      
      drhodp_at_quad = ele_val_at_quad(drhodp, ele)
      eosp_at_quad = ele_val_at_quad(eospressure, ele)
      
      select case(stabilisation_scheme)
        case(STABILISATION_SUPG)
          call transform_to_physical(coordinate, ele, test_shape_ptr, dshape = dtest_t, detwei=detwei, j=j_mat)
          test_shape = make_supg_shape(test_shape_ptr, dtest_t, nlvelocity_at_quad, j_mat, &
            & nu_bar_scheme = nu_bar_scheme, nu_bar_scale = nu_bar_scale)
        case default
          call transform_to_physical(coordinate, ele, detwei=detwei)
          test_shape = test_shape_ptr
          call incref(test_shape)
      end select
      ! Important note: with SUPG the test function derivatives have not been
      ! modified.

      if(cmcget) then
        if(exclude_mass) then
          ele_mat = 0.0
        else
          ele_mat = (1./(dt*dt*theta_divergence*theta_pg))*shape_shape(test_shape, test_shape_ptr, detwei*drhodp_at_quad)
        end if
      end if
      
      if(assemble_rhs) then
        !       /
        ! rhs = |test_shape* &
        !       /
        !      ((1./dt)*(drhodp*(eospressure - (pressure + atmospheric_pressure)) + olddensity - density)
        ! +(absorption)*(drhodp*theta*(eospressure - (pressure + atmospheric_pressure)) 
        !                - theta*density - (1-theta)*olddensity)
        ! +source)dV
        if(exclude_mass) then
          ele_rhs = 0.0
        else
          ele_rhs = (1./dt)*shape_rhs(test_shape, detwei*((drhodp_at_quad*(eosp_at_quad - p_at_quad)) &
                                                      +(olddensity_at_quad - density_at_quad)))
        end if
        
        if(have_source) then
          ele_rhs = ele_rhs + shape_rhs(test_shape, detwei*ele_val_at_quad(source, ele))
        end if
      end if
        
      if(have_absorption) then
        abs_at_quad = ele_val_at_quad(absorption, ele)
        if(cmcget) then
          ele_mat = ele_mat + &
                    (theta/(dt*theta_divergence*theta_pg))*shape_shape(test_shape, test_shape_ptr, &
                                                                       detwei*drhodp_at_quad*abs_at_quad)
        end if
        if(assemble_rhs) then
          ele_rhs = ele_rhs + &
                    shape_rhs(test_shape, detwei*abs_at_quad*(theta*(drhodp_at_quad*(eosp_at_quad - p_at_quad)-density_at_quad) &
                                                            -(1-theta)*olddensity_at_quad))
        end if
      end if
        
      if(assemble_rhs) call addto(rhs, test_nodes, ele_rhs)
      
      if((.not.exclude_mass).and.(cmcget)) then
        call addto(cmc, test_nodes, test_nodes, ele_mat)
      end if
      
      call deallocate(test_shape)
      
    end do

    call deallocate(drhodp)
    call deallocate(eospressure)

  end subroutine assemble_1mat_compressible_projection_cg

  subroutine update_compressible_density(state)
  
    type(state_type), dimension(:), intent(inout) :: state
    
    type(scalar_field), pointer :: density
    
    if((size(state)==1).and.(.not.has_scalar_field(state(1), "MaterialVolumeFraction"))) then
    
      density=>extract_scalar_field(state(1),'Density')
      
      if(have_option(trim(density%option_path)//"/prognostic")) then
        
        call compressible_eos(state(1), density=density)
      
      end if
    
    end if
  
  end subroutine update_compressible_density

  subroutine include_implicit_pressure_buoyancy(ct_m, state, pressure, velocity, get_ct)
    !!< Add the dependency of the buoyancy density on pressure directly into the divergence matrix
    !!< (As this is only used with compressible this is only actually used as a gradient matrix)
    type(block_csr_matrix), intent(inout) :: ct_m
    type(state_type), intent(in) :: state
    type(scalar_field), intent(in) :: pressure
    type(vector_field), intent(in) :: velocity
    logical, intent(in) :: get_ct

    type(scalar_field), pointer :: buoyancy
    type(vector_field), pointer :: gravity, coordinate
    type(scalar_field) :: drhodp

    integer, dimension(:), pointer :: pressure_nodes, velocity_nodes
    real, dimension(:), allocatable :: detwei
    real, dimension(:,:,:), allocatable :: ele_mat

    real :: gravity_magnitude
    integer :: ele, dim, stat
    logical :: have_gravity, have_compressible_eos, implicit_pressure_buoyancy

    ! nothing to do if we're not assembly ct_m
    if(.not.get_ct) return

    call get_option("/physical_parameters/gravity/magnitude", gravity_magnitude, &
         stat=stat)
    have_gravity = stat == 0

    ! nothing to do if we don't have gravity turned on
    if(.not.have_gravity) return

    have_compressible_eos = have_option(trim(state%option_path)//"/equation_of_state/compressible")

    ! nothing to do if we're not compressible
    if(.not.have_compressible_eos) return

    implicit_pressure_buoyancy = have_option(trim(pressure%option_path)// &
                                 "/prognostic/spatial_discretisation/compressible/implicit_pressure_buoyancy")

    ! nothing to do if we don't have the option turned on
    if(.not.implicit_pressure_buoyancy) return
  
    ! made it in!
    ewrite(1,*) 'Including implicit pressure buoyancy term in ct_m'

    buoyancy=>extract_scalar_field(state, "VelocityBuoyancyDensity")
    gravity=>extract_vector_field(state, "GravityDirection", stat)
    coordinate=>extract_vector_field(state, "Coordinate")

    call allocate(drhodp, buoyancy%mesh, "Implicitdrhodp")
    call compressible_eos(state, drhodp=drhodp)

    allocate(detwei(ele_ngi(velocity, 1)), &
             ele_mat(velocity%dim, ele_loc(pressure, 1), ele_loc(velocity, 1)))

    element_loop: do ele = 1, ele_count(pressure)
      pressure_nodes=>ele_nodes(pressure, ele)
      velocity_nodes=>ele_nodes(velocity, ele)

      call transform_to_physical(coordinate, ele, detwei=detwei)

      ele_mat = ele_mat + &
                shape_shape_vector(ele_shape(pressure, ele), ele_shape(velocity, ele), &
                                   detwei*ele_val_at_quad(drhodp, ele)*gravity_magnitude, &
                                   ele_val_at_quad(gravity, ele))

      do dim = 1, velocity%dim
        call addto(ct_m, 1, dim, pressure_nodes, velocity_nodes, ele_mat(dim,:,:))
      end do

    end do element_loop
    
  end subroutine include_implicit_pressure_buoyancy
  
  subroutine compressible_projection_check_options

    integer :: iphase
    character(len=OPTION_PATH_LEN) :: pressure_option_path, density_option_path

    do iphase = 1, option_count("/material_phase")
      pressure_option_path = "/material_phase["//int2str(iphase-1)//"]/scalar_field::Pressure"
      if(have_option(trim(pressure_option_path)//"/prognostic/spatial_discretisation/compressible/implicit_pressure_buoyancy")) then
        if(have_option(trim(pressure_option_path)//"/prognostic/spatial_discretisation/control_volumes")) then
          FLExit("Compressible option implicit_pressure_buoyancy does not work with control volume Pressure discretisations.")
        end if
      end if
      density_option_path = "/material_phase["//int2str(iphase-1)//"]/scalar_field::Density"
      if(have_option(trim(density_option_path)//"/prognostic/spatial_discretisation/use_reference_density")) then
        if(.not.have_option(trim(density_option_path)// &
            "/prognostic/spatial_discretisation/continuous_galerkin/mass_terms/exclude_mass_terms").and. &
           .not.have_option(trim(density_option_path)// &
            "/prognostic/spatial_discretisation/control_volumes/mass_terms/exclude_mass_terms")) then
          FLExit("Using reference density in the continuity equation only really makes sense with exclude_mass_terms.")
        end if
        if(have_option(trim(density_option_path)//"/prognostic/scalar_field::Absorption")) then
          FLExit("Using reference density in the continuity equation doesn't make sense with absorption.")
        end if
      end if
    end do

  end subroutine compressible_projection_check_options

end module compressible_projection

