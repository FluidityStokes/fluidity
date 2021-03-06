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
!    Found

#include "fdebug.h"

module mesh_diagnostics

  use fldebug
  use spud
  use global_parameters, only: FIELD_NAME_LEN
  use halos_numbering
  use fields
  use state_module
  use field_options
  use mesh_quality
  use diagnostic_source_fields

  implicit none
  
  private
  
  public :: calculate_column_ids, calculate_universal_column_ids
  public :: calculate_mesh_quality, calculate_region_ids

contains

  subroutine calculate_column_ids(state, s_field)
    type(state_type), intent(in) :: state
    type(scalar_field), intent(inout) :: s_field
    
    ewrite(1, *) "In calculate_column_ids"

    if(.not.associated(s_field%mesh%columns)) then
      if(have_option(trim(s_field%mesh%option_path)//"/from_mesh/extrude")) then
        FLAbort("No columns associated with an extruded mesh.")
      else
        FLExit("Requested column_id output on non-extruded mesh.")
      end if
    end if
    
    call set_all(s_field, float(s_field%mesh%columns))
    
    ewrite(1, *) "Exiting calculate_column_ids"
  
  end subroutine calculate_column_ids

  subroutine calculate_universal_column_ids(state, s_field)
    type(state_type), intent(in) :: state
    type(scalar_field), intent(inout) :: s_field

    type(mesh_type), pointer :: from_mesh
    character(len=FIELD_NAME_LEN) :: from_mesh_name
    integer :: nhalos
    
    ewrite(1, *) "In calculate_universal_column_ids"

    if(.not.associated(s_field%mesh%columns)) then
      if(have_option(trim(s_field%mesh%option_path)//"/from_mesh/extrude")) then
        FLAbort("No columns associated with an extruded mesh.")
      else
        FLExit("Requested column_id output on non-extruded mesh.")
      end if
    end if

    call get_option(trim(s_field%mesh%option_path)//"/from_mesh/mesh/name", from_mesh_name)
    from_mesh => extract_mesh(state, trim(from_mesh_name))
    nhalos = halo_count(s_field)
    if(nhalos>0) then
      call set_all(s_field, float(halo_universal_numbers(from_mesh%halos(nhalos), s_field%mesh%columns)))
    else
      call set_all(s_field, float(s_field%mesh%columns))
    end if
    
    ewrite(1, *) "Exiting calculate_universal_column_ids"
  
  end subroutine calculate_universal_column_ids

  subroutine calculate_mesh_quality(state, s_field)
    type(state_type), intent(in) :: state
    type(scalar_field), intent(inout) :: s_field

    type(vector_field), pointer :: positions

    integer :: measure
    character(len=FIELD_NAME_LEN) :: measure_name

    positions => extract_vector_field(state,"Coordinate")

    if (element_count(s_field) /= node_count(s_field)) then
       FLAbort("Mesh quality measures must be on a P0 mesh")
    end if

    call get_option(trim(complete_field_path(trim(s_field%option_path))) // &
                    "/algorithm[0]/quality_function/name", measure_name)

    select case(trim(measure_name))
    case("radius_ratio")
       measure=VTK_QUALITY_RADIUS_RATIO
    case("aspect_ratio")
       measure=VTK_QUALITY_ASPECT_RATIO
    case("aspect_frobenius")
       measure=VTK_QUALITY_ASPECT_FROBENIUS
    case("edge_ratio")
       measure=VTK_QUALITY_EDGE_RATIO
    case("condition")
       measure=VTK_QUALITY_CONDITION
    case("min_angle")
       measure=VTK_QUALITY_MIN_ANGLE
    case("max_angle")
       measure=VTK_QUALITY_MAX_ANGLE
    case("shape")
       measure=VTK_QUALITY_SHAPE
    case("min_angl")
       measure=VTK_QUALITY_SHAPE_AND_SIZE
    case("area_or_volume")
       measure=VTK_QUALITY_AREA
    case default
        FLAbort("Unknown quality function for "//trim(s_field%name) )
    end select

    call get_mesh_quality(positions, s_field, measure)

  end subroutine calculate_mesh_quality

  subroutine calculate_region_ids(state, s_field)

    type(state_type), intent(in) :: state
    type(scalar_field), intent(inout) :: s_field

    ewrite(1, *) "In calculate_region_ids"

    if(.not.associated(s_field%mesh%region_ids)) then
      FLAbort("No region ids stored on mesh.")
    end if

    if(.not.(s_field%mesh%shape%degree==0 .and. s_field%mesh%continuity<0)) then
      FLExit("Diagnostic region_ids field should be on a P0(DG) mesh")
    end if

    call set_all(s_field, float(s_field%mesh%region_ids))

    ewrite(1, *) "Exiting calculate_region_ids"

  end subroutine calculate_region_ids

end module mesh_diagnostics
