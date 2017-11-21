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
module shape_functions
  !!< Generate shape functions for elements of arbitrary polynomial degree.
  use FLDebug
  use futils
  use polynomials
  use element_numbering
  use quadrature
  use elements
  use Superconvergence
  use ieee_arithmetic, only: ieee_quiet_nan, ieee_value
  
  implicit none

  private :: lagrange_polynomial, nonconforming_polynomial

  interface make_element_shape
     module procedure make_element_shape_from_element, make_element_shape
  end interface
  
contains

  function make_element_shape_from_element(model, vertices, dim, degree,&
       & quad, type, quad_s, constraint_type_choice, stat)  result (shape)
    !!< This function enables element shapes to be derived from other
    !!< element shapes by specifying which attributes to change.
    type(element_type) :: shape
    type(element_type), intent(in) :: model
    !! Vertices is the number of vertices of the element, not the number of nodes!
    !! dim may be 1, 2, or 3.
    !! Degree is the degree of the Lagrange polynomials.
    integer, intent(in), optional :: vertices, dim, degree
    type(quadrature_type), intent(in), target, optional :: quad
    integer, intent(in), optional :: type
    type(quadrature_type), intent(in), optional, target :: quad_s
    !! Element constraints
    integer, intent(in), optional :: constraint_type_choice
    integer, intent(out), optional :: stat

    integer :: lvertices, ldim, ldegree, lconstraint_type_choice
    type(quadrature_type) :: lquad
    type(quadrature_type), pointer :: lquad_s
    integer :: ltype

    if (present(vertices)) then
       lvertices=vertices
    else
       lvertices=model%numbering%vertices
    end if

    if (present(dim)) then
       ldim=dim
    else
       ldim=model%dim
    end if

    if(present(degree)) then
       ldegree=degree
    else
       ldegree=model%degree
    end if

    if(present(quad)) then
       lquad=quad
    else
       lquad=model%quadrature
    end if

    if(present(type)) then
       ltype=type
    else
       ltype=model%numbering%type
    end if

    if(present(quad_s)) then
       lquad_s=>quad_s
    else if (associated(model%surface_quadrature)) then
       lquad_s=>model%surface_quadrature
    else
       lquad_s=>null()
    end if

    if(present(constraint_type_choice)) then
       lconstraint_type_choice=constraint_type_choice
    else if (associated(model%constraints)) then
       lconstraint_type_choice=model%constraints%type
    else
       lconstraint_type_choice=CONSTRAINT_NONE
    end if

    
    if (associated(lquad_s)) then
       shape = make_element_shape(lvertices, ldim, ldegree, lquad, ltype,&
            lquad_s, constraint_type_choice=lconstraint_type_choice, stat=stat)
    else
       shape = make_element_shape(lvertices, ldim, ldegree, lquad, ltype,&
            constraint_type_choice=lconstraint_type_choice, stat=stat)
    end if

  end function make_element_shape_from_element

  function make_element_shape(vertices, dim, degree, quad, type,&
       quad_s, constraint_type_choice, stat)  result (shape)
    !!< Generate the shape functions for an element. The result is a suitable
    !!< element_type.
    !!
    !!< At this stage only Lagrange family polynomial elements are supported.
    type(element_type) :: shape
    !! Vertices is the number of vertices of the element, not the number of nodes!
    !! dim \in [1,2,3] is currently supported.
    !! Degree is the degree of the Lagrange polynomials.
    integer, intent(in) :: vertices, dim, degree
    type(quadrature_type), intent(in), target :: quad
    integer, intent(in), optional :: type
    type(quadrature_type), intent(in), optional, target :: quad_s
    integer, intent(in), optional :: constraint_type_choice
    integer, intent(out), optional :: stat

    real, pointer :: g(:)=> null()

    type(ele_numbering_type), pointer :: ele_num
    ! Count coordinates of each point 
    integer, dimension(dim+1) :: counts
    integer :: i,j,k
    integer :: ltype, coords
    real :: dx
    type(constraints_type), pointer :: constraint

    ! Check that the quadrature and the element shapes match.
    assert(quad%vertices==vertices)
    assert(quad%dim==dim)

    if (present(type)) then
       ltype=type
    else
       ltype=ELEMENT_LAGRANGIAN
    end if

    if (present(stat)) stat=0

    ! Get the local numbering of our element
    ele_num=>find_element_numbering(vertices, dim, degree, type)

    if (.not.associated(ele_num)) then
       if (present(stat)) then
          stat=1
          return
       else
          FLAbort('Element numbering unavailable.')
       end if
    end if

    shape%numbering=>ele_num
    shape%quadrature=quad
    call incref(quad)

    ! The number of local coordinates depends on the element family.
    select case(ele_num%family)
    case (FAMILY_SIMPLEX)
       coords=dim+1
    case (FAMILY_CUBE)
       if(ele_num%type==ELEMENT_TRACE .and. dim==2) then
          !For trace elements the local coordinate is face number
          !then the local coordinates on the face
          !For quads, the face is an interval element which has
          !two local coordinates.
          coords=3
       else
          coords=dim
       end if
    case default
       FLAbort('Illegal element family.')
    end select

    if (present(quad_s) .and. ele_num%type/=ELEMENT_TRACE .and. ele_num%family==FAMILY_SIMPLEX) then
       allocate(shape%surface_quadrature)
       shape%surface_quadrature=quad_s
       call incref(quad_s)
       call allocate(shape, ele_num, quad%ngi, ngi_s=quad_s%ngi)
       shape%n_s=0.0
       shape%dn_s=0.0
    else
       call allocate(shape, ele_num, quad%ngi)
    end if
    shape%degree=degree
    shape%n=0.0
    shape%dn=0.0

    ! Construct shape for each node
    do i=1,shape%loc

       counts(1:coords)=ele_num%number2count(:,i)

       ! Construct appropriate polynomials.
       do j=1,coords

          select case(ltype)
          case(ELEMENT_LAGRANGIAN)
             if (degree == 0) then
                dx = 0.0
             else
                dx = 1.0/degree
             end if 
             select case(ele_num%family)
             case (FAMILY_SIMPLEX)
                ! Raw polynomial.
                shape%spoly(j,i)&
                     =lagrange_polynomial(counts(j), counts(j), dx)
             case(FAMILY_CUBE)
                ! note that local coordinates run from -1.0 to 1.0
                shape%spoly(j,i)&
                     =lagrange_polynomial(counts(j), degree, 2.0*dx, &
                     origin=-1.0)
             end select

          case(ELEMENT_TRACE)
             shape%spoly(j,i) = (/ieee_value(0.0,ieee_quiet_nan)/)

          case(ELEMENT_BUBBLE)
             if(i==shape%loc) then

                ! the last node is the bubble shape function
                shape%spoly(j,i) = (/1.0, 0.0/)

             else

                select case(ele_num%family)
                case (FAMILY_SIMPLEX)
                   ! Raw polynomial.
                   shape%spoly(j,i)&
                        =lagrange_polynomial(counts(j)/coords, counts(j)/coords, 1.0/degree)

                end select

             end if

          case(ELEMENT_NONCONFORMING)

             shape%spoly(j,i)=nonconforming_polynomial(counts(j))

          case default

             FLAbort('An unsupported element type has been selected.')

          end select

          ! Derivative
          if(ele_num%type==ELEMENT_TRACE) then
             shape%dspoly(j,i) = (/ieee_value(0.0,ieee_quiet_nan)/)
          else
             shape%dspoly(j,i)=ddx(shape%spoly(j,i))
          end if
       end do

       if(ele_num%type==ELEMENT_TRACE) then
          !No interior functions, hence NaNs
          shape%n = ieee_value(0.0,ieee_quiet_nan)
          shape%dn = ieee_value(0.0,ieee_quiet_nan)
       else
          ! Loop over all the quadrature points.
          do j=1,quad%ngi

             ! Raw shape function
             shape%n(i,j)=eval_shape(shape, i, quad%l(j,:))

             ! Directional derivatives.
             shape%dn(i,j,:)=eval_dshape(shape, i, quad%l(j,:))
          end do

          if (present(quad_s)) then
             select case(ele_num%family)
             case(FAMILY_SIMPLEX)
                allocate(g(dim+1))
                do j=1,quad_s%ngi
                   g(1) = 0.0
                   do k=1,dim
                      g(k+1)=quad_s%l(j,k)
                   end do
                   ! In order to match the arbitrary face node ordering
                   ! these must get reoriented before use so we don't care
                   ! about which local facet they're with respect to.
                   shape%n_s(i,j)=eval_shape(shape,i,g)
                   shape%dn_s(i,j,:)=eval_dshape(shape,i,g)
                end do
                deallocate(g)
             end select
          end if
       end if
    end do

    if(ele_num%type.ne.ELEMENT_TRACE) then
       shape%superconvergence => get_superconvergence(shape)
    end if

    if(present(constraint_type_choice)) then
       if(constraint_type_choice/=CONSTRAINT_NONE) then
          allocate(constraint)
          shape%constraints=>constraint
          call allocate(shape%constraints,shape,constraint_type_choice)
       end if
    end if

  end function make_element_shape

  function lagrange_polynomial(n,degree,dx, origin) result (poly)
    ! nth equispaced lagrange polynomial of specified degree and point
    ! spacing dx.
    integer, intent(in) :: n, degree
    real, intent(in) :: dx
    type(polynomial) :: poly
    ! fixes location of n=0 location (0.0 if not specified)
    real, intent(in), optional :: origin

    real lorigin
    integer :: i

    ! This shouldn't be necessary but there appears to be a bug in initial
    ! component values in gfortran:
    poly%coefs=>null()
    poly%degree=-1
    
    if (present(origin)) then
       lorigin=origin
    else
       lorigin=0.0
    end if

    poly=(/1.0/)
    
    degreeloop: do i=0,degree
       if (i==n) cycle degreeloop

       poly=poly*(/1.0, -(lorigin+i*dx) /)

    end do degreeloop

    ! normalize to 1.0 in the n-th location
    poly=poly/eval(poly, lorigin+n*dx)
    
  end function lagrange_polynomial

  function nonconforming_polynomial(n) result (poly)
    ! nth P1 nonconforming polynomial.
    integer, intent(in) :: n
    type(polynomial) :: poly

    ! This shouldn't be necessary but there appears to be a bug in initial
    ! component values in gfortran:
    poly%coefs=>null()
    poly%degree=-1

    poly=(/1.0/)    

    if (n==0) then

       ! polynomial is -2x+1
       poly=(/-2.0, 1.0/)

    end if
       
  end function nonconforming_polynomial

  function transformation_lagrangian_to_monotonic_p2(ele_num) result (matrix)
    !!< Returns a matrix that transforms coefficients with respect to the standard Lagrangian P2 basis
    !!< to coefficients wrt the monotonic P2 basis. The transpose of this provides the linear combination
    !!< of Lagrangian P2 basis functions that form the monotonic P2 basis
    type(ele_numbering_type), intent(in) :: ele_num
    real, dimension(ele_num%nodes, ele_num%nodes) :: matrix

    integer, dimension(ele_num%vertices) :: vertices
    integer, dimension(1) :: edge_node
    integer :: i, j

    matrix = 0.
    vertices = local_vertices(ele_num)
    do i=1, size(vertices)
      ! the vertex shape function is the sum of the corr. nodal lagrangian shape function
      ! (which isn't monotonic)
      matrix(vertices(i), vertices(i)) = 1.0

      ! plus a quarter of the adjacent edge shape functions to fix monotonicity
      do j=1, size(vertices)
        if (i==j) exit
        edge_node = edge_local_num((/ i, j /), ele_num, interior=.true.)

        matrix(edge_node(1), vertices(i)) = 0.25

        ! the edge basis function itself is scaled such that the sum of all basis functions is still 1
        matrix(edge_node(1), edge_node(1)) = 0.5
      end do
    end do

  end function transformation_lagrangian_to_monotonic_p2

  function monotonic_p2_shape(lagrangian_p2_shape) result (shape)
    !!< Given the standard lagrangian P2 shape, creates a P2 shape 
    !!< that uses a non-nodal basis which is monotonic. The basis 
    !!< functions are still associated with the standard P2 nodes,
    !!< in the sense that each basis function is 1 in the 
    !!< associated P2 node - however the vertex basis functions is
    !!< not zero in the adjacent edge P2 nodes. The basis functions are
    !!< however still continuous between elements.
    type(element_type), intent(in):: lagrangian_p2_shape
    type(element_type) :: shape

    type(ele_numbering_type), pointer :: ele_num
    real, dimension(lagrangian_p2_shape%loc, lagrangian_p2_shape%loc) :: matrix
    integer :: i, j, k

    ele_num => lagrangian_p2_shape%numbering
    assert(ele_num%type == FAMILY_SIMPLEX)
    assert(ele_num%family == FAMILY_SIMPLEX)
    assert(lagrangian_p2_shape%degree == 2)
    assert(.not. associated(lagrangian_p2_shape%superconvergence))
    assert(.not. associated(lagrangian_p2_shape%constraints))

    ! This follows the same order of setting the shape attributes
    ! as in make_element_shape: we first set numbering and quadrature
    ! /before/ calling allocate (which is intent(inout)!) - yuck!
    shape%numbering => ele_num
    shape%quadrature = lagrangian_p2_shape%quadrature
    call incref(shape%quadrature)

    if (associated(lagrangian_p2_shape%surface_quadrature)) then
      allocate(shape%surface_quadrature)
      shape%surface_quadrature = lagrangian_p2_shape%surface_quadrature
      call incref(shape%surface_quadrature)
      call allocate(shape, ele_num, shape%quadrature%ngi, ngi_s=shape%surface_quadrature%ngi)
    else
      call allocate(shape, ele_num, shape%quadrature%ngi)
    end if
    shape%degree = 2

    matrix = transpose(transformation_lagrangian_to_monotonic_p2(ele_num))

    shape%n = matmul(matrix, lagrangian_p2_shape%n)
    do i=1, size(shape%dn, 3)
      shape%dn(:,:,i) = matmul(matrix, lagrangian_p2_shape%dn(:,:,i))
    end do
    do i=1, size(matrix, 1)
      do j=1, size(shape%spoly, 1)
        shape%spoly(j,i) = (/ 0. /)
        do k=1, size(matrix, 2)
          shape%spoly(j,i) = shape%spoly(j,i) + matrix(i, k) * lagrangian_p2_shape%spoly(j,k)
          shape%dspoly(j,i) = shape%dspoly(j,i) + matrix(i, k) * lagrangian_p2_shape%dspoly(j,k)
        end do
      end do
    end do

    if (associated(shape%surface_quadrature)) then
      shape%n_s = matmul(matrix, lagrangian_p2_shape%n_s)
      do i=1, size(shape%dn_s, 3)
        shape%dn_s(:,:,i) = matmul(matrix, lagrangian_p2_shape%dn_s(:,:,i))
      end do
    end if

    ! This indicates that the basis functions do *not* have the property: 1 in associated node
    ! 0 in other nodes. Thus we cannot equate "value in node" with coeffient. We do however assume
    ! that any linear combination of basis functions is uniquely determined by the values in all
    ! nodes of the element numbering
    shape%nodal = .false.

  end function monotonic_p2_shape

end module shape_functions
