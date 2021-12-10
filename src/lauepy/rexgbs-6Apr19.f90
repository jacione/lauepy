module var

  !  program rexgbs is at line number 845, approx.
  
  integer :: numsamp , numsymm
  integer, parameter :: numsigmas = 300
  integer :: numcsl , nsig(numsigmas)
  double precision :: quatsymm(4,48) , quatsamp(4,4)
  double precision :: quatcsl(4,50)
  integer, parameter :: norients = 1000

end module var

!+{+{+{+{+{+{+{+{+{+{+{+{+{+{+{+{+{+{+{+{+{+{+{+{+{+{+{+{+{+{+{+{+{+{+{+
module precision

  ! Real kinds

  integer, parameter :: kr4 = selected_real_kind(6,37)       ! single precision real
  integer, parameter :: kr8 = selected_real_kind(15,307)     ! double precision real

  ! Integer kinds

  integer, parameter :: ki4 = selected_int_kind(9)           ! single precision integer
  integer, parameter :: ki8 = selected_int_kind(18)          ! double precision integer

  !Complex kinds

  integer, parameter :: kc4 = kr4                            ! single precision complex
  integer, parameter :: kc8 = kr8                            ! double precision complex

end module precision

!..............................................
module strings

  use precision

  private :: value_dr,value_sr,value_di,value_si
  private :: write_dr,write_sr,write_di,write_si
  private :: writeq_dr,writeq_sr,writeq_di,writeq_si

  interface value  ! Generic operator for converting a number string to a 
     ! number. Calling syntax is 'call value(numstring,number,ios)' 
     ! where 'numstring' is a number string and 'number' is a 
     ! real number or an integer (single or double precision).         
     module procedure value_dr
     module procedure value_sr
     module procedure value_di
     module procedure value_si
  end interface value

  interface writenum  ! Generic  interface for writing a number to a string. The 
     ! number is left justified in the string. The calling syntax
     ! is 'call writenum(number,string,format)' where 'number' is
     ! a real number or an integer, 'string' is a character string
     ! containing the result, and 'format' is the format desired, 
     ! e.g., 'e15.6' or 'i5'.
     module procedure write_dr
     module procedure write_sr
     module procedure write_di
     module procedure write_si
  end interface writenum

  interface writeq  ! Generic interface equating a name to a numerical value. The
     ! calling syntax is 'call writeq(unit,name,value,format)' where
     ! unit is the integer output unit number, 'name' is the variable
     ! name, 'value' is the real or integer value of the variable, 
     ! and 'format' is the format of the value. The result written to
     ! the output unit has the form <name> = <value>.
     module procedure writeq_dr
     module procedure writeq_sr
     module procedure writeq_di
     module procedure writeq_si
  end interface writeq


  !**********************************************************************

contains

  !**********************************************************************

  subroutine parse(str,delims,args,nargs)

    ! Parses the string 'str' into arguments args(1), ..., args(nargs) based on
    ! the delimiters contained in the string 'delims'. Preceding a delimiter in
    ! 'str' by a backslash (\) makes this particular instance not a delimiter.
    ! The integer output variable nargs contains the number of arguments found.

    character(len=*) :: str,delims
    character(len=len_trim(str)) :: strsav
    character(len=*),dimension(:) :: args

    strsav=str
    call compact(str)
    na=size(args)
    do i=1,na
       args(i)=' '
    end do
    nargs=0
    lenstr=len_trim(str)
    if(lenstr==0) return
    k=0

    do
       if(len_trim(str) == 0) exit
       nargs=nargs+1
       call split(str,delims,args(nargs))
       call removebksl(args(nargs))
    end do
    str=strsav

  end subroutine parse

  !**********************************************************************

  subroutine compact(str)

    ! Converts multiple spaces and tabs to single spaces; deletes control characters;
    ! removes initial spaces.

    character(len=*):: str
    character(len=1):: ch
    character(len=len_trim(str)):: outstr

    str=adjustl(str)
    lenstr=len_trim(str)
    outstr=' '
    isp=0
    k=0

    do i=1,lenstr
       ch=str(i:i)
       ich=iachar(ch)

       select case(ich)

       case(9,32)     ! space or tab character
          if(isp==0) then
             k=k+1
             outstr(k:k)=' '
          end if
          isp=1

       case(33:)      ! not a space, quote, or control character
          k=k+1
          outstr(k:k)=ch
          isp=0

       end select

    end do

    str=adjustl(outstr)

  end subroutine compact

  !**********************************************************************

  subroutine removesp(str)

    ! Removes spaces, tabs, and control characters in string str

    character(len=*):: str
    character(len=1):: ch
    character(len=len_trim(str))::outstr

    str=adjustl(str)
    lenstr=len_trim(str)
    outstr=' '
    k=0

    do i=1,lenstr
       ch=str(i:i)
       ich=iachar(ch)
       select case(ich)    
       case(0:32)  ! space, tab, or control character
          cycle       
       case(33:)  
          k=k+1
          outstr(k:k)=ch
       end select
    end do

    str=adjustl(outstr)

  end subroutine removesp

  !**********************************************************************

  subroutine value_dr(str,rnum,ios)

    ! Converts number string to a double precision real number

    character(len=*)::str
    real(kr8)::rnum
    integer :: ios

    ilen=len_trim(str)
    ipos=scan(str,'Ee')
    if(.not.is_digit(str(ilen:ilen)) .and. ipos/=0) then
       ios=3
       return
    end if
    read(str,*,iostat=ios) rnum

  end subroutine value_dr

  !**********************************************************************

  subroutine value_sr(str,rnum,ios)

    ! Converts number string to a single precision real number

    character(len=*)::str
    real(kr4) :: rnum
    real(kr8) :: rnumd 

    call value_dr(str,rnumd,ios)
    if( abs(rnumd) > huge(rnum) ) then
       ios=15
       return
    end if
    if( abs(rnumd) < tiny(rnum) ) rnum=0.0_kr4
    rnum=rnumd

  end subroutine value_sr

  !**********************************************************************

  subroutine value_di(str,inum,ios)

    ! Converts number string to a double precision integer value

    character(len=*)::str
    integer(ki8) :: inum
    real(kr8) :: rnum

    call value_dr(str,rnum,ios)
    if(abs(rnum)>huge(inum)) then
       ios=15
       return
    end if
    inum=nint(rnum,ki8)

  end subroutine value_di

  !**********************************************************************

  subroutine value_si(str,inum,ios)

    ! Converts number string to a single precision integer value

    character(len=*)::str
    integer(ki4) :: inum
    real(kr8) :: rnum

    call value_dr(str,rnum,ios)
    if(abs(rnum)>huge(inum)) then
       ios=15
       return
    end if
    inum=nint(rnum,ki4)

  end subroutine value_si

  !**********************************************************************

  subroutine shiftstr(str,n)

    ! Shifts characters in in the string 'str' n positions (positive values
    ! denote a right shift and negative values denote a left shift). Characters
    ! that are shifted off the end are lost. Positions opened up by the shift 
    ! are replaced by spaces.

    character(len=*):: str

    lenstr=len(str)
    nabs=iabs(n)
    if(nabs>=lenstr) then
       str=repeat(' ',lenstr)
       return
    end if
    if(n<0) str=str(nabs+1:)//repeat(' ',nabs)  ! shift left
    if(n>0) str=repeat(' ',nabs)//str(:lenstr-nabs)  ! shift right 
    return

  end subroutine shiftstr

  !**********************************************************************

  subroutine insertstr(str,strins,loc)

    ! Inserts the string 'strins' into the string 'str' at position 'loc'. 
    ! Characters in 'str' starting at position 'loc' are shifted right to
    ! make room for the inserted string. Trailing spaces of 'strins' are 
    ! removed prior to insertion

    character(len=*):: str,strins
    character(len=len(str))::tempstr

    lenstrins=len_trim(strins)
    tempstr=str(loc:)
    call shiftstr(tempstr,lenstrins)
    tempstr(1:lenstrins)=strins(1:lenstrins)
    str(loc:)=tempstr
    return

  end subroutine insertstr

  !**********************************************************************

  subroutine delsubstr(str,substr)

    ! Deletes first occurrence of substring 'substr' from string 'str' and
    ! shifts characters left to fill hole. Trailing spaces or blanks are
    ! not considered part of 'substr'.

    character(len=*):: str,substr

    lensubstr=len_trim(substr)
    ipos=index(str,substr)
    if(ipos==0) return
    if(ipos == 1) then
       str=str(lensubstr+1:)
    else
       str=str(:ipos-1)//str(ipos+lensubstr:)
    end if
    return

  end subroutine delsubstr

  !**********************************************************************

  subroutine delall(str,substr)

    ! Deletes all occurrences of substring 'substr' from string 'str' and
    ! shifts characters left to fill holes.

    character(len=*):: str,substr

    lensubstr=len_trim(substr)
    do
       ipos=index(str,substr)
       if(ipos == 0) exit
       if(ipos == 1) then
          str=str(lensubstr+1:)
       else
          str=str(:ipos-1)//str(ipos+lensubstr:)
       end if
    end do
    return

  end subroutine delall

  !**********************************************************************

  function uppercase(str) result(ucstr)

    ! convert string to upper case

    character (len=*):: str
    character (len=len_trim(str)):: ucstr

    ilen=len_trim(str)
    ioffset=iachar('A')-iachar('a')     
    iquote=0
    ucstr=str
    do i=1,ilen
       iav=iachar(str(i:i))
       if(iquote==0 .and. (iav==34 .or.iav==39)) then
          iquote=1
          iqc=iav
          cycle
       end if
       if(iquote==1 .and. iav==iqc) then
          iquote=0
          cycle
       end if
       if (iquote==1) cycle
       if(iav >= iachar('a') .and. iav <= iachar('z')) then
          ucstr(i:i)=achar(iav+ioffset)
       else
          ucstr(i:i)=str(i:i)
       end if
    end do
    return

  end function uppercase

  !**********************************************************************

  function lowercase(str) result(lcstr)

    ! convert string to lower case

    character (len=*):: str
    character (len=len_trim(str)):: lcstr

    ilen=len_trim(str)
    ioffset=iachar('A')-iachar('a')
    iquote=0
    lcstr=str
    do i=1,ilen
       iav=iachar(str(i:i))
       if(iquote==0 .and. (iav==34 .or.iav==39)) then
          iquote=1
          iqc=iav
          cycle
       end if
       if(iquote==1 .and. iav==iqc) then
          iquote=0
          cycle
       end if
       if (iquote==1) cycle
       if(iav >= iachar('A') .and. iav <= iachar('Z')) then
          lcstr(i:i)=achar(iav-ioffset)
       else
          lcstr(i:i)=str(i:i)
       end if
    end do
    return

  end function lowercase

  !**********************************************************************

  subroutine readline(nunitr,line,ios)

    ! Reads line from unit=nunitr, ignoring blank lines
    ! and deleting comments beginning with an exclamation point(!)

    character (len=*):: line

    do  
       read(nunitr,'(a)', iostat=ios) line      ! read input line
       if(ios /= 0) return
       line=adjustl(line)
       ipos=index(line,'!')
       if(ipos == 1) cycle
       if(ipos /= 0) line=line(:ipos-1)
       if(len_trim(line) /= 0) exit
    end do
    return

  end subroutine readline

  !**********************************************************************

  subroutine match(str,ipos,imatch)

    ! Sets imatch to the position in string of the delimiter matching the delimiter
    ! in position ipos. Allowable delimiters are (), [], {}, <>.

    character(len=*) :: str
    character :: delim1,delim2,ch

    lenstr=len_trim(str)
    delim1=str(ipos:ipos)
    select case(delim1)
    case('(')
       idelim2=iachar(delim1)+1
       istart=ipos+1
       iend=lenstr
       inc=1
    case(')')
       idelim2=iachar(delim1)-1
       istart=ipos-1
       iend=1
       inc=-1
    case('[','{','<')
       idelim2=iachar(delim1)+2
       istart=ipos+1
       iend=lenstr
       inc=1
    case(']','}','>')
       idelim2=iachar(delim1)-2
       istart=ipos-1
       iend=1
       inc=-1
    case default
       write(*,*) delim1,' is not a valid delimiter'
       return
    end select
    if(istart < 1 .or. istart > lenstr) then
       write(*,*) delim1,' has no matching delimiter'
       return
    end if
    delim2=achar(idelim2) ! matching delimiter

    isum=1
    do i=istart,iend,inc
       ch=str(i:i)
       if(ch /= delim1 .and. ch /= delim2) cycle
       if(ch == delim1) isum=isum+1
       if(ch == delim2) isum=isum-1
       if(isum == 0) exit
    end do
    if(isum /= 0) then
       write(*,*) delim1,' has no matching delimiter'
       return
    end if
    imatch=i

    return

  end subroutine match

  !**********************************************************************

  subroutine write_dr(rnum,str,fmt)

    ! Writes double precision real number rnum to string str using format fmt

    real(kr8) :: rnum
    character(len=*) :: str,fmt
    character(len=80) :: formt

    formt='('//trim(fmt)//')'
    write(str,formt) rnum
    str=adjustl(str)

  end subroutine write_dr

  !***********************************************************************

  subroutine write_sr(rnum,str,fmt)

    ! Writes single precision real number rnum to string str using format fmt

    real(kr4) :: rnum
    character(len=*) :: str,fmt
    character(len=80) :: formt

    formt='('//trim(fmt)//')'
    write(str,formt) rnum
    str=adjustl(str)

  end subroutine write_sr

  !***********************************************************************

  subroutine write_di(inum,str,fmt)

    ! Writes double precision integer inum to string str using format fmt

    integer(ki8) :: inum
    character(len=*) :: str,fmt
    character(len=80) :: formt

    formt='('//trim(fmt)//')'
    write(str,formt) inum
    str=adjustl(str)

  end subroutine write_di

  !***********************************************************************

  subroutine write_si(inum,str,fmt)

    ! Writes single precision integer inum to string str using format fmt

    integer(ki4) :: inum
    character(len=*) :: str,fmt
    character(len=80) :: formt

    formt='('//trim(fmt)//')'
    write(str,formt) inum
    str=adjustl(str)

  end subroutine write_si

  !***********************************************************************

  subroutine trimzero(str)

    ! Deletes nonsignificant trailing zeroes from number string str. If number
    ! string ends in a decimal point, one trailing zero is added.

    character(len=*) :: str
    character :: ch
    character(len=10) :: exp

    ipos=scan(str,'eE')
    if(ipos>0) then
       exp=str(ipos:)
       str=str(1:ipos-1)
    endif
    lstr=len_trim(str)
    do i=lstr,1,-1
       ch=str(i:i)
       if(ch=='0') cycle          
       if(ch=='.') then
          str=str(1:i)//'0'
          if(ipos>0) str=trim(str)//trim(exp)
          exit
       endif
       str=str(1:i)
       exit
    end do
    if(ipos>0) str=trim(str)//trim(exp)

  end subroutine trimzero

  !**********************************************************************

  subroutine writeq_dr(unit,namestr,value,fmt)

    ! Writes a string of the form <name> = value to unit

    real(kr8) :: value
    integer :: unit
    character(len=*) :: namestr,fmt
    character(len=32) :: tempstr

    call writenum(value,tempstr,fmt)
    call trimzero(tempstr)
    write(unit,*) trim(namestr)//' = '//trim(tempstr)

  end subroutine writeq_dr

  !**********************************************************************

  subroutine writeq_sr(unit,namestr,value,fmt)

    ! Writes a string of the form <name> = value to unit

    real(kr4) :: value
    integer :: unit
    character(len=*) :: namestr,fmt
    character(len=32) :: tempstr

    call writenum(value,tempstr,fmt)
    call trimzero(tempstr)
    write(unit,*) trim(namestr)//' = '//trim(tempstr)

  end subroutine writeq_sr

  !**********************************************************************

  subroutine writeq_di(unit,namestr,ivalue,fmt)

    ! Writes a string of the form <name> = ivalue to unit

    integer(ki8) :: ivalue
    integer :: unit
    character(len=*) :: namestr,fmt
    character(len=32) :: tempstr
    call writenum(ivalue,tempstr,fmt)
    call trimzero(tempstr)
    write(unit,*) trim(namestr)//' = '//trim(tempstr)

  end subroutine writeq_di

  !**********************************************************************

  subroutine writeq_si(unit,namestr,ivalue,fmt)

    ! Writes a string of the form <name> = ivalue to unit

    integer(ki4) :: ivalue
    integer :: unit
    character(len=*) :: namestr,fmt
    character(len=32) :: tempstr
    call writenum(ivalue,tempstr,fmt)
    call trimzero(tempstr)
    write(unit,*) trim(namestr)//' = '//trim(tempstr)

  end subroutine writeq_si

  !**********************************************************************

  function is_letter(ch) result(res)

    ! Returns .true. if ch is a letter and .false. otherwise

    character :: ch
    logical :: res

    select case(ch)
    case('A':'Z','a':'z')
       res=.true.
    case default
       res=.false.
    end select
    return

  end function is_letter

  !**********************************************************************

  function is_digit(ch) result(res)

    ! Returns .true. if ch is a digit (0,1,...,9) and .false. otherwise

    character :: ch
    logical :: res

    select case(ch)
    case('0':'9')
       res=.true.
    case default
       res=.false.
    end select
    return

  end function is_digit

  !**********************************************************************

  subroutine split(str,delims,before,sep)

    ! Routine finds the first instance of a character from 'delims' in the
    ! the string 'str'. The characters before the found delimiter are
    ! output in 'before'. The characters after the found delimiter are
    ! output in 'str'. The optional output character 'sep' contains the 
    ! found delimiter. A delimiter in 'str' is treated like an ordinary 
    ! character if it is preceded by a backslash (\). If the backslash 
    ! character is desired in 'str', then precede it with another backslash.

    character(len=*) :: str,delims,before
    character,optional :: sep
    logical :: pres
    character :: ch,cha

    pres=present(sep)
    str=adjustl(str)
    call compact(str)
    lenstr=len_trim(str)
    if(lenstr == 0) return        ! string str is empty
    k=0
    ibsl=0                        ! backslash initially inactive
    before=' '
    do i=1,lenstr
       ch=str(i:i)
       if(ibsl == 1) then          ! backslash active
          k=k+1
          before(k:k)=ch
          ibsl=0
          cycle
       end if
       if(ch == '\') then          ! backslash with backslash inactive
          k=k+1
          before(k:k)=ch
          ibsl=1
          cycle
       end if
       ipos=index(delims,ch)         
       if(ipos == 0) then          ! character is not a delimiter
          k=k+1
          before(k:k)=ch
          cycle
       end if
       if(ch /= ' ') then          ! character is a delimiter that is not a space
          str=str(i+1:)
          if(pres) sep=ch
          exit
       end if
       cha=str(i+1:i+1)            ! character is a space delimiter
       iposa=index(delims,cha)
       if(iposa > 0) then          ! next character is a delimiter
          str=str(i+2:)
          if(pres) sep=cha
          exit
       else
          str=str(i+1:)
          if(pres) sep=ch
          exit
       end if
    end do
    if(i >= lenstr) str=''
    str=adjustl(str)              ! remove initial spaces
    return

  end subroutine split

  !**********************************************************************

  subroutine removebksl(str)

    ! Removes backslash (\) characters. Double backslashes (\\) are replaced
    ! by a single backslash.

    character(len=*):: str
    character(len=1):: ch
    character(len=len_trim(str))::outstr

    str=adjustl(str)
    lenstr=len_trim(str)
    outstr=' '
    k=0
    ibsl=0                        ! backslash initially inactive

    do i=1,lenstr
       ch=str(i:i)
       if(ibsl == 1) then          ! backslash active
          k=k+1
          outstr(k:k)=ch
          ibsl=0
          cycle
       end if
       if(ch == '\') then          ! backslash with backslash inactive
          ibsl=1
          cycle
       end if
       k=k+1
       outstr(k:k)=ch              ! non-backslash with backslash inactive
    end do

    str=adjustl(outstr)

    end subroutine removebksl

  !**********************************************************************

  subroutine ReplaceCommaWithDot(str)

    ! Replaces any comma with a dot
    ! intended for making European floats into American

    character(len=*):: str
    character(len=1):: ch
    character(len=len_trim(str))::outstr

    str=adjustl(str)
    lenstr=len_trim(str)

    do i=1,lenstr
       if( str(i:i) == ',' ) str(i:i) = '.'
    end do

  end subroutine ReplaceCommaWithDot

!**********************************************************************

end module strings  
!+{+{+{+{+{+{+{+{+{+{+{+{+{+{+{+{+{+{+{+{+{+{+{+{+{+{+{+{+{+{+{+{+{+{+{+


!+++++++++++++++++++++

!  Check the section where OutputForAnalysis.txt is written
!   the version on the old MacBook may have the "near100, near110 "etc included
!   as well as the fix for the textur2 input of pairs of Eulers

program rexgbs
  use var
  implicit none

  !  compile with
  !   gfortran -ffixed-line-length-none -O3 -fno-align-commons -o rexgbs rexgbs-[DATE].f90
  !  or
  !   gfortran -fbounds-check -Wall -gdwarf-2 -O0 -fno-align-commons -o rexgbs-debug rexgbs-[DATE].f90

  !      quatn(4) = rtmp * 0.25  !  CHECK THIS!  line 1465

  double precision :: scalar
  double precision, parameter :: pi  = 4.d0 * atan ( 1.0d0 ) 
  double precision, parameter :: pi2  = 8.d0 * atan ( 1.0d0 ) 
  double precision, parameter :: rad = 180. / pi , rad2 = 360. / pi
  double precision :: rfvector(3) , rfvector2(3)
  character :: afile*200 , bfile*200

  double precision :: gbnorma(3),gbdirna(3),gbnormb(3),gbdirnb(3),Bunge(3)
  double precision :: matrixA(3,3) , matrixB(3,3) , matrixTwin(3,3) ,matrix(3,3)

!!$  double precision :: quatsymm(4,48) ,quatsamp(4,4)
!!$  integer :: numsymm ,numsamp
!!$  common/a1/ quatsymm,numsymm,quatsamp,numsamp
!!$  double precision :: quatcsl(4,50)
!!$  integer :: numcsl ,nsig(50)
!!$  common/a2/ quatcsl ,numcsl ,nsig

  !     common/a3/aquat(4,100),ngrain,phi1(100),capphi(100),phi2(100)
  double precision :: aquat(4,norients) ,phi1(norients),capphi(norients),phi2(norients),grwt(norients)
  double precision :: p1(2) , pp(2) , p2(2)
  integer :: aspin(norients) , ngrain
  common/a3/ aquat ,phi1 ,capphi ,phi2 ,grwt , aspin , ngrain

  !  for the 2nd list
  double precision :: bquat(4,norients),bphi1(norients),bcapphi(norients),bphi2(norients),bgrwt(norients)
  integer :: bspin(norients) , bngrain
  common/a4/ bquat , bphi1 , bcapphi , bphi2 , bgrwt , bspin , bngrain

  integer :: iq_symm
  character :: inline*120

  integer :: iq_extra , iex , ulimit  !  for insertion of extra symm. oper.
  integer :: iq_source  !  0=terminal, 1=file
  double precision :: qtmp(4) , qtmp1(4) , qtmp2(4) , d1 , d2 , d3 , rnorm
  double precision :: qdisq1(4) , qdisq2(4)
  double precision :: samp1(3) , samp2(3)  ! for testing rotations of the misorient axis
  double precision :: qgrain1(4) , qgrain2(4) ,  qgrain3(4) , qgrain4(4)
  integer :: combine

  integer :: ior, i, ij, ijk, i2, iq, iq2, j, jk, k, k1, k2, jkl
  integer :: iqPairedList
  ! temporary control of getting a sequence of orientation pairs from a file
  !   rather than calculating all combinations in a list of orientations
  integer :: index1, index1a, index2, index2a
  double precision :: ang2 , shear !  needed for twin analysis
  double precision :: HPI(3,3)
  !      data HPI/-1.,0.,0.,0.,-1.,0.,0.,0.,1./
  double precision :: b3(3) , b4(4) , n3(3) , n4(4) , car , bar  !  lattice param ratios
  double precision, parameter :: sqr3 = sqrt(3.)
  integer :: iq_choice
  double precision :: omegaNorm
  double precision :: lat , long  !  for GMT plotting
  double precision :: MillerAxis(4) , MillerAxisSorted(4) , factor  !  misor axis in Miller indices
  integer :: nsize , iq_data_type , innerLoopStart , innerLoopLimit , outerLoopLimit , OuterLoopIncrement
  double precision :: minMisor  !  tracks min misor angle
  integer :: jindexMin
  double precision :: misaxis(3) , axis(3) , qr1(4) , qr2(4) , qtwo( 4 , 2 ) , qresult(4)
  double precision :: QtwoForSampleFrame( 4 , 2 ), AxisForSampleFrame(3)
  ! see last output section near line # 2347
  double precision :: tmpvec(3) , tmpout(3)
  !  for storing results between multiple pairs
  double precision :: qresulta( norients , norients , 4 ) , gbaxis( 3 , norients , norients ) , rnear( norients , norients )
  double precision :: QresultForSampleFrame( norients , norients , 4 )
  ! see last output section near line # 2347
  integer :: lsig( norients , norients )
  double precision :: angle( norients , norients )
  integer*4 ::  ac_min(3)
  double precision :: angle_input, brandon, crit, theta, thetmin, tmp
  double precision :: vec100(3) , vec110(3), vec111(3) ! test nearness to simple axes
  double precision :: rtmp100, rtmp110 , rtmp111
  double precision, parameter :: QuatReverseTwin(1:4) = (/ -0.28867513459, -0.28867513459, -0.28867513459, 0.86602540378 /) 
  double precision, parameter :: QuatTwin(1:4) = (/ 0.28867513459, 0.28867513459, 0.28867513459, 0.86602540378 /)
  double precision :: QuatAfterTwinApplied(4)
  double precision :: AngleAfterTwinApplied, AxisAfterTwinApplied(3)
  
  !  CODE:
  call load_CSL_quat( .false. )

  vec100(1) = 0.; vec100(2) = 0.; vec100(3) = 1.
  vec110(1) = 0.; vec110(2) = 1.; vec110(3) = 1.
  vec111(1) = 1.; vec111(2) = 1.; vec111(3) = 1.
  call vecnorm(vec100) ; call vecnorm(vec110) ; call vecnorm(vec111)
  
  iq_extra = 0  !  keep it turned off by default
  do i = 1 , 3
     do j = 1 , 3
        hpi(i,j) = 0.
        if( i.eq.j ) hpi(i,j) = -1.
        if( i.eq.j .and. i.eq.3 ) hpi(i,j) = 1.
     end do
  end do
  !     CALL GetArg(1, afile)
  print*,'Welcome to REXGBS'
  print*,'Last edited 16 Jan 16, ADR'
  write(*,*) 'Calculates a set of misorientations etc, from a list'
  write(*,*) ' of orientations in a file, or from the command line.'
  print*,'Be sure that you include '
  write(*,*) ' all the sample-symmetry related variants for each'
  write(*,*) ' component or you will not get the full picture!'

  !     read(*,*) iq
  !iq = 0  !  this governs how to get the input

  if ( iargc() <= 3 ) then
     print *,'REXGBS: No input file given as argument: input from terminal'
     print*,' Sample usage:'
     print*,'[/code/texture/]rexgbs symmetry[0|1|2|3|4] input-type[0-8] data_source[0|1] data_file_type[0|1|2] file_with_data'
     print*
     print*, 'Crystal symmetry: 0=cubic 1=Hex 2=Orth 3=Tetr 4=Tetr_4/m'
     print*
     print*,'input-type (ONLY for terminal input):'
     write(*,*) 'Pair of Euler angles in radians [ = 0 ]'
     write(*,*) ' or a pair of Euler angles in degrees [ = 1 ]'
     print*,' or a single axis-angle misorientation [ = 2 ]?'
     write(*,*) ' or as UVW-HKL-A + UVW-HKL-B  [ = 3 ]?'
     print*,' or as a single Rodrigues-Frank vector [ = 4 ]?'
     print*,' or as a PAIR of Rodrigues-Frank vectors [ = 5 ]?'
     print*,' or as a PAIR of axis-angle orientations [ = 6 ]?'
     print*,' or as a PAIR of quaternions [ = 7 ]?'
     print*,' or as a twin (reflection or rotation) [ = 8 ]?'
     print*
     print*, 'data_source choice: enter orientations directly from terminal [=0]?'
     print*, 'or, read orientations from a file [=1]'
     !print*, 'or, work on a pair of lists [=2]?'
     print*
     print*,' data_file_type can be text with 2 Euler sets per line, Bunge only [=0]'
     print*,' or a TEXIN/WTS file with 1 Euler set per line [=1]'
     print*,' or a Sukbin-style file with 1 Euler set per line as a quaternion [=2]'
  end if

  iq_symm = 0
  CALL GetArg( 1 , inline )
  if( inline == '' ) then
     print *,' No symmetry as argument: input from terminal'
     print*,'Enter xtal symmetry: 0=cubic 1=Hex 2=Orth 3=Tetr 4=Tetr_4/m'
     read* , iq_symm
  else
     !         print*,'second argument =  ',inline
     read(inline,"(i1)") iq_symm
     !         print*,'iq_symm = ',iq_symm
     if(iq_symm >= 0 .and. iq_symm < 5 ) then
        print*,'Crystal symmetry: 0=cubic 1=Hex 2=Orth 3=Tetr 4=Tetr_4/m'
        print*,'Symmetry choice = ' , iq_symm
     else
        stop 'invalid symmetry choice'
     end if
  end if
  call load_orthorhombic_sample_symmetry_quat( .false. )
  if( iq_symm == 0 ) call load_cubic_crystal_symmetry_quat( .false. )
  if( iq_symm == 1 ) call load_hcp_crystal_symmetry_quat( .false. )
  if( iq_symm == 2 ) call load_orthorhombic_symmetry_quat( .false. )
  if( iq_symm == 3 ) stop 'apologies, need to add tetr symm'
  if( iq_symm == 4 ) stop 'apologies, need to add tetr_4/m symm'
  
  CALL GetArg( 2 , inline )
  if( inline == '' ) then  !  no input
     write(*,*) 'Choose a parameterization: '
     write(*,*) 'A single entry defines the misorientation '
     write(*,*) 'A paired entry defines two orientations:'
     write(*,*) 'Pair of Euler angles in radians [ = 0 ]'
     write(*,*) ' or a pair of Euler angles in degrees [ = 1 ]'
     print*,' or a single axis-angle misorientation [ = 2 ]?'
     write(*,*) ' or as UVW-HKL-A + UVW-HKL-B  [ = 3 ]?'
     print*,' or as a single Rodrigues-Frank vector [ = 4 ]?'
     print*,' or as a PAIR of Rodrigues-Frank vectors [ = 5 ]?'
     print*,' or as a PAIR of axis-angle orientations [ = 6 ]?'
     print*,' or as a PAIR of quaternions [ = 7 ]?'
     print*,' or as a twin (reflection or rotation) [ = 8 ]?'
     print*,'BEWARE: remember that a pair with the same '
     print*,'orientation is treated as zero misorientation!'
     print*
     read(*,*) iq2
  else
     read(inline, * ) iq2
  end if
  if( iq2 .ge. 0 .and. iq2 .le. 8 ) then  !  7 -> 8, 31 xii 15
     print* , 'input type selected = ' , iq2
  else
     stop 'invalid input type choice'
  end if

  CALL GetArg( 3 , inline )
  if( inline == '' ) then
     print *,' Which source of data?  Command-line [0] or file [1] or pair of files [2]?'
     read*, iq
  else
     read(inline, * ) iq
  end if
  if( iq >= 0 .and. iq < 3 ) then  !  fixed to allow command line, 31 xii 15, ADR
     print* , 'data source selected = ' , iq
     if ( iq == 2 ) print*,' two filenames must be on same line with a space in between the names'     
  else
     stop 'invalid input source choice'
  end if

  iq_data_type = -99  !  assign a definite value (for terminal input)
  if ( iq > 0 ) then
     CALL GetArg( 4 , inline )
     if( inline == '' ) then
        print *,' Which type of data file?  text with Euler pairs [0] or TEXIN/WTS file [1] or Sukbin file [2]?'
        read*, iq_data_type
     else
        read( inline , * ) iq_data_type
        if( iq_data_type >= 0 .and. iq_data_type < 3 ) then
           print* , 'data file type selected = ' , iq_data_type
        else
           stop 'invalid data file type choice'
        end if
     end if
     if( iq_data_type == 0 ) print*,'Will compute misorientations for each pair (only), one per input line'
  
     if ( iq > 0 ) then
        CALL GetArg( 5 , afile )  ! if there are command line arguments, this will not work for 2 files
        if( afile == '' ) then
           if ( iq == 1 ) then
              print*,'Enter the filename of the dataset'
              read*, afile
              print*, 'input file name with data: ', afile
           end if
           if ( iq == 2 ) then
              print*,'Enter the filenames of the two datasets'
              read*, afile, bfile
              print*, 'input filenames with data: ', afile, bfile
           end if
!!$     else
!!$        read ( inline, '(a)' ) afile
        end if  !  if( afile == '' ) 
     end if    !  if ( iq > 0 )
  end if

!!$  CALL GetArg( 3 , bfile )
!!$  if( bfile .eq. '' ) then
!!$     print *,' No 2nd input file given as argument'
!!$     print*,' normal calculation'
!!$  else 
!!$     print *,'2nd input file given as argument; will search for matches'
!!$     iq = 2
!!$  end if

!!$  iq_data_type = 1
!!$  if ( iq == 2 ) then
!!$     CALL GetArg( 4 , inline )
!!$     if( inline.eq.'' .and. iq /= 1 ) then
!!$        print *,' No input type as argument'
!!$        print*,'Will assume TEXIN-like (LApp) input'
!!$        iq_data_type = 1
!!$     else
!!$        !         print*,'second argument =  ',inline
!!$        read(inline,"(i1)") iq_data_type
!!$        !         print*,'iq_data_type = ',iq_data_type
!!$        if( iq_data_type <= 0 .or. iq_data_type >= 3 ) then
!!$           stop 'invalid choice of input data type'
!!$        endif
!!$        if( iq_data_type == 1 ) print*,'texin/wts input'
!!$        if( iq_data_type == 2 ) print*,'Sukbin quat list input'
!!$     end if
!!$  end if

  !call textur4( iq_symm )     !  read in the symmetry operators

  print*,'CAUTION!!!'
  print*,'If you are combining 2 successive transformations'
  print*,'as in, e.g., twin chains, you must be v careful'
  print*,'about the order.'
  print*,'Remember that the second rotation is negated'
  print*,'whenever a pair of orientations is combined.'
  print*,'Accordingly, you should give the negative (inverse)'
  print*,'of the 2nd rotation to get the desired effect.'
  print*,'I.e. enter gA, then -gB.'
  print*,'To reverse the order, enter gB, then -gA'
  print*,'________________________________________'
  print*


  if ( iq == 1 ) call textur2( afile , iq_data_type )

  if ( iq == 0 ) then
     !  only if you go with command line input can you get all the variations in input type
     ngrain=2  !  limit to 2 grains
!!$     CALL GetArg(4, inline)
!!$     if(inline.eq.'') then
        Print*,'Insert additional symmetry operation [ yes=1 ]'
        read*,iq_extra
!!$     else  !  no argument present
!!$        read(inline, * ) iq_extra
!!$        if( iq_extra .eq. 1 ) then
!!$           if ( iq2 .eq. 2 .or. iq2 .eq. 4 ) then
!!$              print*,'Sorry, must enter TWO orientations'
!!$              stop
!!$           endif
!!$           print*,'will insert extra symmetry operator '
!!$        end if
!!$     end if

     if( iq2.eq.0 .or. iq2.eq.1 .or. (iq2.ge.5 .and. iq2.le.7 ) ) then
        print*,'Combine the two into one misor [ yes = 1 ]?'
        read*,combine
     end if

     if( iq2 .eq. 0 ) then
        write(*,*) 'Enter phi1A PhiA phi2A phi1B PhiB phi2B [rad]'
        read(*,*) phi1(1),capphi(1),phi2(1) ,phi1(2),capphi(2),phi2(2)
        call quatB(phi1(1),capphi(1),phi2(1),qr1)
        call quatB(phi1(2),capphi(2),phi2(2),qr2)
     endif
     if( iq2 .eq. 1 ) then
        write(*,*) 'Enter phi1A PhiA phi2A phi1B PhiB phi2B [degr]'
        read(*,*) phi1(1),capphi(1),phi2(1) ,phi1(2),capphi(2),phi2(2)
        call quatB(phi1(1)/rad,capphi(1)/rad,phi2(1)/rad,qr1)
        call quatB(phi1(2)/rad,capphi(2)/rad,phi2(2)/rad,qr2)
     endif
     if( iq2 .eq. 0 .or. iq2 .eq. 1 ) then
        do ijk=1,4
           aquat(ijk,1)=qr1(ijk)
           aquat(ijk,2)=qr2(ijk)
        enddo
        !         endif

     elseif( iq2 .eq. 2 ) then
        write(*,*) 'Enter Angle(degr) axis1, axis2, axis3'
        read(*,*) angle_input,misaxis
        call donorm(misaxis) ! normalize
        phi1(1)=0.
        capphi(1)=0.
        phi2(1)=0.
        aquat(1,1)=0.
        aquat(2,1)=0.
        aquat(3,1)=0.
        aquat(4,1)=1.
        aquat(1,2) = misaxis(1) * sin(angle_input*pi/360.)
        aquat(2,2)=misaxis(2)*sin(angle_input*pi/360.)
        aquat(3,2)=misaxis(3)*sin(angle_input*pi/360.)
        aquat(4,2)=cos(angle_input*pi/360.)
        do ijk=1,4
           qr2(ijk)=aquat(ijk,2)
        enddo
        write(*,"('quaternion = ',4f8.3,2i4)") qr2
        call q2eulB(phi1(2),capphi(2),phi2(2),qr2)
        phi1(2)=phi1(2)*180./pi
        capphi(2)=capphi(2)*180./pi
        phi2(2)=phi2(2)*180./pi
        if(phi1(2).gt.360.) phi1(2)=phi1(2)-360.
        if(phi1(2).lt.0.) phi1(2)=phi1(2)+360.
        if(phi2(2).gt.360.) phi2(2)=phi2(2)-360.
        if(phi2(2).lt.0.) phi2(2)=phi2(2)+360.

     elseif( iq2 .eq. 3 ) then
        print*,'Enter  H K L  U V W   for xtal A'
        read*,(gbnorma(i),i=1,3),(gbdirna(i),i=1,3)
        call uvwhkl2quat(gbnorma,gbdirna,qr1,Bunge,matrixA)
        do i = 1,4
           aquat(i,1) = qr1(i)
        enddo
        phi1(1) = Bunge(1)
        capphi(1) = Bunge(2)
        phi2(1) = Bunge(3)
        print*
        print*,'Orientation matrix for grain A'
        do i = 1,3
           print"('[ ',3(1x,f8.3),' ]')",(matrixA(i,j),j=1,3)
        enddo

        print*,'Enter  H K L    U V W   for xtal B'
        read*,(gbnormb(i),i=1,3),(gbdirnb(i),i=1,3)
        !  here we re-use code from rexgbs.f (aug 06)
        call uvwhkl2quat(gbnormb,gbdirnb,qr2,Bunge,matrixB)
        phi1(2) = Bunge(1)
        capphi(2) = Bunge(2)
        phi2(2) = Bunge(3)
        do i = 1,4
           aquat(i,2) = qr2(i)
        enddo
        print*,'Orientation matrix for grain B'
        do i = 1,3
           print"('[ ',3(1x,f8.3),' ]')",(matrixB(i,j),j=1,3)
        enddo
        print*
        print*,'Misorientation matrix from A to B'
        do i = 1,3
           do j = 1,3
              matrix(i,j) = 0.
              do k = 1,3
                 matrix(i,j) = matrix(i,j) + matrixB(i,k)*matrixA(j,k)
                 !  note inversion for matrix A
              enddo
           enddo
           print"('[ ',3(1x,f8.3),' ]')",(matrix(i,j),j=1,3)
        enddo
        tmp = (-1.d0 + matrix(1,1)+matrix(2,2)+matrix(3,3)) * 0.5d0
        if(tmp.lt.-1.) tmp = -1.
        if(tmp.gt.1.) tmp = 1.
        print*,'Angle = ',acos(tmp) * rad

        !  RF vector entry (single)
     elseif( iq2 .eq. 4 ) then
        write(*,*) 'Enter R1  R2  R3'
        read(*,*) rfvector
        call rod2q(rfvector,qr2)
        phi1(1) = 0.
        capphi(1) = 0.
        phi2(1) = 0.
        aquat(1,1) = 0.
        aquat(2,1) = 0.
        aquat(3,1) = 0.
        aquat(4,1) = 1.

        aquat(1,2) = qr2(1)
        aquat(2,2) = qr2(2)
        aquat(3,2) = qr2(3)
        aquat(4,2) = qr2(4)
        write(*,"('quaternion = ',4f8.3,2i4)") qr2
        call q2eulB(phi1(2),capphi(2),phi2(2),qr2)
        phi1(2)=phi1(2)*180./pi
        capphi(2)=capphi(2)*180./pi
        phi2(2)=phi2(2)*180./pi
        if(phi1(2).gt.360.) phi1(2)=phi1(2)-360.
        if(phi1(2).lt.0.) phi1(2)=phi1(2)+360.
        if(phi2(2).gt.360.) phi2(2)=phi2(2)-360.
        if(phi2(2).lt.0.) phi2(2)=phi2(2)+360.

        !  pair of RF vectors
     elseif( iq2.eq. 5 ) then
        write(*,*) 'Enter R1  R2  R3 (for A), R1  R2  R3 (for B)'
        read(*,*) rfvector,rfvector2

        call rod2q(rfvector,qr2)
        aquat(1,1) = qr2(1)
        aquat(2,1) = qr2(2)
        aquat(3,1) = qr2(3)
        aquat(4,1) = qr2(4)
        write(*,"('quaternion for A = ',4f8.3,2i4)") qr2
        call q2eulB(phi1(1),capphi(1),phi2(1),qr2)
        phi1(1)=phi1(1)*180./pi
        capphi(1)=capphi(1)*180./pi
        phi2(1)=phi2(1)*180./pi
        if(phi1(1).gt.360.) phi1(1)=phi1(1)-360.
        if(phi1(1).lt.0.) phi1(1)=phi1(1)+360.
        if(phi2(1).gt.360.) phi2(1)=phi2(1)-360.
        if(phi2(1).lt.0.) phi2(1)=phi2(1)+360.

        call rod2q(rfvector2,qr2)
        aquat(1,2) = qr2(1)
        aquat(2,2) = qr2(2)
        aquat(3,2) = qr2(3)
        aquat(4,2) = qr2(4)
        write(*,"('quaternion for B = ',4f8.3,2i4)") qr2
        call q2eulB(phi1(2),capphi(2),phi2(2),qr2)
        phi1(2) = phi1(2)*180./pi
        capphi(2) = capphi(2)*180./pi
        phi2(2) = phi2(2)*180./pi
        if(phi1(2).gt.360.) phi1(2) = phi1(2)-360.
        if(phi1(2).lt.0.) phi1(2) = phi1(2)+360.
        if(phi2(2).gt.360.) phi2(2) = phi2(2)-360.
        if(phi2(2).lt.0.) phi2(2) = phi2(2)+360.

        !  pair of axis-angle specs.
     elseif( iq2 .eq. 6 ) then
        write(*,*) 'Enter 1st Angle(degr) axis1, axis2, axis3'
        read(*,*) angle_input,misaxis
        call donorm(misaxis) ! normalize
        aquat(1,1) = misaxis(1)*sin(angle_input / rad2)
        aquat(2,1) = misaxis(2)*sin(angle_input / rad2)
        aquat(3,1) = misaxis(3)*sin(angle_input / rad2)
        aquat(4,1) = cos(angle_input / rad2 )
        do ijk = 1 , 4
           qr1(ijk) = aquat(ijk,1)
        enddo
        write(*,"('1st quaternion = ',4f8.3,2i4)") qr1
        call q2eulB(phi1(1),capphi(1),phi2(1),qr1)
        phi1(1) = phi1(1) * rad
        capphi(1) = capphi(1) * rad
        phi2(1) = phi2(1) * rad
        if(phi1(1).gt.360.) phi1(2)=phi1(2)-360.
        if(phi1(1).lt.0.) phi1(2)=phi1(2)+360.
        if(phi2(1).gt.360.) phi2(2)=phi2(2)-360.
        if(phi2(1).lt.0.) phi2(2)=phi2(2)+360.
        write(*,*) 'Enter 2nd Angle(degr) axis1, axis2, axis3'
        read(*,*) angle_input,misaxis
        call donorm(misaxis) ! normalize
        aquat(1,2) = misaxis(1)*sin( angle_input / rad2 )
        aquat(2,2) = misaxis(2)*sin( angle_input / rad2 )
        aquat(3,2) = misaxis(3)*sin( angle_input / rad2 )
        aquat(4,2) = cos( angle_input / rad2 )
        do ijk=1,4
           qr2(ijk) = aquat(ijk,2)
        enddo
        write(*,"('2nd quaternion = ',4f8.3,2i4)") qr2
        call q2eulB(phi1(2),capphi(2),phi2(2),qr2)
        phi1(2)=phi1(2)* rad
        capphi(2)=capphi(2)* rad
        phi2(2)=phi2(2)* rad
        if(phi1(2).gt.360.) phi1(2)=phi1(2)-360.
        if(phi1(2).lt.0.) phi1(2)=phi1(2)+360.
        if(phi2(2).gt.360.) phi2(2)=phi2(2)-360.
        if(phi2(2).lt.0.) phi2(2)=phi2(2)+360.

        !  pair of quaternions
     elseif( iq2.eq. 7 ) then
        print*,'Enter 2 quaternions, with cos(t/2) as q4;'
        print*,'Enter q1 q2 q3 q4 (for A), q1 q2 q3 q4 (for B):'
        read(*,*) qr1 , qr2
        call qnorm(qr1)
        call qnorm(qr2)
        print*,'The quaternions were normalized'
        aquat(1,1) = qr1(1)
        aquat(2,1) = qr1(2)
        aquat(3,1) = qr1(3)
        aquat(4,1) = qr1(4)
        write(*,"('check: quaternion for A = ',4f8.3,2i4)") qr1
        call q2eulB(phi1(1),capphi(1),phi2(1),qr2)
        phi1(1)=phi1(1)*180./pi
        capphi(1)=capphi(1)*180./pi
        phi2(1)=phi2(1)*180./pi
        if(phi1(1).gt.360.) phi1(1)=phi1(1)-360.
        if(phi1(1).lt.0.) phi1(1)=phi1(1)+360.
        if(phi2(1).gt.360.) phi2(1)=phi2(1)-360.
        if(phi2(1).lt.0.) phi2(1)=phi2(1)+360.

        aquat(1,2) = qr2(1)
        aquat(2,2) = qr2(2)
        aquat(3,2) = qr2(3)
        aquat(4,2) = qr2(4)
        write(*,"('check: quaternion for B = ',4f8.3,2i4)") qr2
        call q2eulB(phi1(2),capphi(2),phi2(2),qr2)
        phi1(2) = phi1(2)*180./pi
        capphi(2) = capphi(2)*180./pi
        phi2(2) = phi2(2)*180./pi
        if(phi1(2).gt.360.) phi1(2) = phi1(2)-360.
        if(phi1(2).lt.0.) phi1(2) = phi1(2)+360.
        if(phi2(2).gt.360.) phi2(2) = phi2(2)-360.
        if(phi2(2).lt.0.) phi2(2) = phi2(2)+360.

        !  pair of quaternions
     elseif( iq2.eq. 8 ) then
        open ( 13 , file = 'REXGBS-pstext1.txt',status = 'unknown' )
        print* , 'Use DrawStereograms.sh to plot REXGBS-pstext1.txt'
        write( 13 , * ) '.Twin     output for pstext'
        write( 13 ,"(' 7. 3. 10 0 0 CM X')" )
        write( 13 ,"(' 85. 4. 10 0 0 CM Y')" )
        write( 13 ,"(' 60. 80. 10 0 0 CM Z')" )
        if ( iq_symm >= 1 .or. iq_symm <= 4 ) then
           print*,'enter c over a ratio'
           read*, car
        end if
        if ( iq_symm .eq. 2 ) then
           print*,'enter b over a ratio'
           read*, bar
        end if
        phi1(1) = 0.
        capphi(1) = 0.
        phi2(1) = 0.
        aquat(1,1) = 0.
        aquat(2,1) = 0.
        aquat(3,1) = 0.
        aquat(4,1) = 1.
        print*
        print*,'Choose Tome method [1] (OK for hcp)'
        print* , ' or 180 about K1, for reflection twin [2]'
        print* , ' or 180 about eta1, for rotation twin [3]'
        print* , 'note: either method is OK for compound twins'
        print* , 'note: both direction and plane input for plotting'
        !            print*,', or b+n+shear method [x] (NOT advised at present)'
        !  everything related to this last method is commented out because
        !    it will not give the correct result; twin reorientation must be
        !    treated differently than slip
        read * , iq_choice
        if ( iq_choice < 1 .or. iq_choice > 3 ) then
           print * , ' invalid choice!'
           call exit(1)
        end if

        if ( iq_choice == 1 ) then
           print*,'The twin orientation is given by a 180 degrees '
           print*,'rotation around the Burgers vector of the '
           print* , 'active twinning system [from Carlos Tome]'
        else if ( iq_choice == 2 ) then
           print*,'The twin orientation is given by a 180 degrees '
           print*,'rotation around K1. '
        else if ( iq_choice == 3 ) then
           print*,'The twin orientation is given by a 180 degrees '
           print*,'rotation around eta1. '
        end if
        print*,'Enter b1 b2 b3 '
        print*,'or  b1 b2 b3 b4  if HCP :'
        if ( iq_symm .eq. 1 ) then
           read * , b4
        else
           read * , b3
        end if
        print*,'Enter n1 n2 n3 '
        print*,'or  n1 n2 n3 n4  if HCP :'
        if ( iq_symm .eq. 1 ) then
           read * , n4
        else
           read * , n3
        end if

        !  code from: subroutine twinor(ktw,ng,nomen,dbca)  in lapp68.large.f
        !  from examining VPSC7, it is clear that Carlos+Ricardo use the Bunge
        !  convention as their standard
        if ( iq_symm == 1 ) then   !  hex
           gbnormb(1) = n4(1)
           gbnormb(2) = (n4(1)+2.*n4(2))/sqr3
           gbnormb(3) = n4(4)/car
           gbdirnb(1)=1.5*b4(1)
           gbdirnb(2)=(b4(1)/2.+b4(2))*sqr3
           gbdirnb(3)=b4(4)*car
        else if ( iq_symm == 2 ) then   !  orthorhombic
           gbnormb(1) = n3(1)
           gbnormb(2) = n3(2) / bar
           gbnormb(3) = n3(3) / car
           gbdirnb(1) = b3(1)
           gbdirnb(2) = b3(2) * bar
           gbdirnb(3) = b3(3) * car
        else if ( iq_symm >= 3 .and. iq_symm <= 4 ) then   !  tetragonal
           gbnormb(1) = n3(1)
           gbnormb(2) = n3(2)
           gbnormb(3) = n3(3) / car
           gbdirnb(1) = b3(1)
           gbdirnb(2) = b3(2)
           gbdirnb(3) = b3(3) * car
        else if ( iq_symm == 0 ) then   !  cubic
           gbnormb(1) = n3(1)
           gbnormb(2) = n3(2)
           gbnormb(3) = n3(3)
           gbdirnb(1) = b3(1)
           gbdirnb(2) = b3(2)
           gbdirnb(3) = b3(3)
        end if  !  symmetry
        call vecnorm ( gbnormb )
        call vecnorm ( gbdirnb )

        rnorm = scalar ( gbnormb , gbdirnb )
        if ( abs(rnorm) > 1.e-5 ) then
           print* , ' plane & direction not perpendicular?!'
           call exit(1)
        end if

        if ( iq_choice == 1 ) then
           call axisang2rot( gbnormb , pi , matrixTwin ) ! 180 about the direction
           print*,'Matrix for Rotation Twin from Axis-Angle'
           do i = 1,3
              print"('[ ',3(1x,f8.3),' ]')",(matrixTwin(i,j),j=1,3)
           end do
           call euler( matrixTwin , 1 , 'B' , phi1(1) , capphi(1) , phi2(1) , ior , 0 )
           write(*,"('Eulers from matrixTwin: ',3(1x,f8.3))") phi1(1) , capphi(1) , phi2(1)
           phi1(2) = atan2( gbdirnb(2) , gbdirnb(1) ) + pi2
           ang2 = sqrt( gbdirnb(1)**2 + gbdirnb(2)**2)
           capphi(2) = atan2( ang2 , gbdirnb(3) )
           phi2(2) = 0.0
           write(*,"('check: Eulers: ',4(1x,f8.3))") phi1(2)*rad , capphi(2)*rad , 0.
           call euler( matrixB , 2 , 'B' , phi1(2)*rad , capphi(2)*rad , phi2(2)*rad , ior , 0 )
           print*,'Orientation matrix for twin from Eulers'
           do i = 1,3
              print"('[ ',3(1x,f8.3),' ]')",(matrixB(i,j),j=1,3)
           end do
           !      call euler(a,iopt,nomen,Bunge(1),Bunge(2),Bunge(3),ior,kerr)
           !cUFK      Note: I'm not sure Carlos isn't assuming a particular NOMEN here
           !  but, see above for Euler angle info that says it's Bunge
           DO I=1,3
              DO J=1,3
                 matrixA(I,J)=0.
                 DO K1 = 1 , 3
                    DO K2 = 1 , 3
                       matrixA(I,J) = matrixA(I,J) + ( matrixB(K1,I) * HPI(K1,K2) * matrixB(K2,J) )
                    END DO
                 END DO
              END DO
           END DO
           print*,'Orientation matrix for misor'
           do i = 1,3
              print"('matrixA[',3(1x,f8.3),' ]')" ,(matrixA(i,j),j=1,3)
           enddo
           call euler(matrixA,1,'B',phi1(2), capphi(2), phi2(2) , ior , 0 )
           call quatB( phi1(2)/rad , capphi(2)/rad , phi2(2)/rad , qr2 )
           write(*,"('check: quaternion: ',4f8.3,2i4)") qr2
           aquat(1,2) = qr2(1)
           aquat(2,2) = qr2(2)
           aquat(3,2) = qr2(3)
           aquat(4,2) = qr2(4)
           !     note about this method: I got the same results as noted in the
           !     spreadsheet that Nathalie Bozzolo gave me for Ti
           !     However when I check for cubic {111}<211> twins, I get the wrong
           !     answer!
!!!$            elseif ( iq_choice .eq. 1 ) then  !  make up matrix from b+n
!!!$               call vecnorm ( gbdirnb )
!!!$               call vecnorm ( gbnormb )
!!!$               DO I = 1 , 3
!!!$                  DO J = 1 , 3
!!!$                     matrix(i,j) = shear * 0.5
!!!$     $                    * ( ( gbdirnb(i) * gbnormb(j) )
!!!$     $                    - ( gbdirnb(j) * gbnormb(i) ) )
!!!$                  END DO
!!!$               END DO
!!!$               print*,' matrix: '
!!!$               do i = 1,3
!!!$                  print"('matrix[',3(1x,f8.3),' ]')"
!!!$     $                 ,(matrix(i,j),j=1,3)
!!!$               enddo
!!!$               omegaNorm = sqrt( matrix(1,2)**2 + matrix(2,3)**2
!!!$     $              + matrix(1,3)**2 )
!!!$               print*,'angle (rads, degrs): ' , omegaNorm
!!!$     $              , omegaNorm * rad
!!!$               DO I = 1 , 3
!!!$                  DO J = 1 , 3
!!!$                     matrix(i,j) = matrix(i,j) / omegaNorm
!!!$                  END DO
!!!$               END DO
!!!$               print*,' matrix==omega-hat: '
!!!$               do i = 1,3
!!!$                  print"('matrix[',3(1x,f8.3),' ]')"
!!!$     $                 ,(matrix(i,j),j=1,3)
!!!$               enddo
!!!$               !  matrix is now omega-hat
!!!$               DO I = 1 , 3
!!!$                  DO J = 1 , 3
!!!$                     matrixB(i,j) = 0.
!!!$                     DO k = 1 , 3
!!!$                        matrixB(i,j) = matrixB(i,j)
!!!$     $                       + matrix(i,k) * matrix(k,j)
!!!$                     END DO
!!!$                  END DO
!!!$               END DO
!!!$               print*,' matrixB==omega-hat^2: '
!!!$               do i = 1,3
!!!$                  print"('matrixB[',3(1x,f8.3),' ]')"
!!!$     $                 ,(matrixB(i,j),j=1,3)
!!!$               enddo
!!!$               !  matrixB is now omega-hat^2
!!!$               DO I = 1 , 3
!!!$                  DO J = 1 , 3
!!!$                     if ( i .eq. j ) then
!!!$                        matrixA(i,j) = 1.
!!!$                     else
!!!$                        matrixA(i,j) = 0.
!!!$                     end if
!!!$                     matrixA(i,j) = matrixA(i,j)
!!!$     $                    + ( matrix(i,j) * sin(omegaNorm) )
!!!$     $                    + ( matrixB(i,j) * ( 1. - cos(omegaNorm) ) )
!!!$                  END DO
!!!$               END DO
!!!$
!!!$               print*,'Orientation matrix for misor'
!!!$               do i = 1,3
!!!$                  print"('matrixA[',3(1x,f8.3),' ]')"
!!!$     $                 ,(matrixA(i,j),j=1,3)
!!!$               end do
!!!$               call euler(matrixA,1,'B',phi1(2), capphi(2), phi2(2)
!!!$     $              , ior , 0 )
!!!$               write(*,"('check: Eulers: ',4(1x,f8.3))") phi1(2)
!!!$     $              , capphi(2) , phi2(2)
!!!$               call quatB( phi1(2)/rad , capphi(2)/rad , phi2(2)/rad
!!!$     $              , qr2 )
!!!$               write(*,"('check: quaternion: ',4f8.3,2i4)") qr2
!!!$               aquat(1,2) = qr2(1)
!!!$               aquat(2,2) = qr2(2)
!!!$               aquat(3,2) = qr2(3)
!!!$               aquat(4,2) = qr2(4)

        elseif ( iq_choice == 2 ) then !  make up matrix from n only
           call axisang2rot( gbnormb , pi , matrixTwin ) ! 180 about the direction
           print*,'Matrix for Rotation Twin from Axis-Angle'
           do i = 1,3
              print"('[ ',3(1x,f8.3),' ]')",(matrixTwin(i,j),j=1,3)
           end do
           call euler( matrixTwin , 1 , 'B' , phi1(1) , capphi(1) , phi2(1) , ior , 0 )
           write(*,"('Eulers from matrixTwin: ',3(1x,f8.3))") phi1(1) , capphi(1) , phi2(1)
           phi1(1)=0.
           capphi(1)=0.
           phi2(1)=0.
           aquat(1,1)=0.
           aquat(2,1)=0.
           aquat(3,1)=0.
           aquat(4,1)=1.
           aquat(1,2) = gbnormb(1)
           aquat(2,2) = gbnormb(2)
           aquat(3,2) = gbnormb(3)
           aquat(4,2) = 0.
           do ijk=1,4
              qr2(ijk)=aquat(ijk,2)
           enddo
           write(*,"('quaternion = ',4(1x,f8.3),2i4)") qr2
           call q2eulB(phi1(2),capphi(2),phi2(2),qr2)
           phi1(2) = phi1(2)*180./pi
           capphi(2) = capphi(2)*180./pi
           phi2(2) = phi2(2)*180./pi
           if(phi1(2).gt.360.) phi1(2) = phi1(2)-360.
           if(phi1(2).lt.0.) phi1(2) = phi1(2)+360.
           if(phi2(2).gt.360.) phi2(2) = phi2(2)-360.
           if(phi2(2).lt.0.) phi2(2) = phi2(2)+360.
           write(*,"('check: Eulers: ',4(1x,f8.3))") phi1(2) , capphi(2) , phi2(2)
           call euler( matrixB , 2 , 'B' , phi1(2)*rad , capphi(2)*rad , phi2(2)*rad , ior , 0 )
           print*,'Matrix for Reflection Twin from Eulers'
           do i = 1,3
              print"('[ ',3(1x,f8.3),' ]')",(matrixB(i,j),j=1,3)
           end do
        elseif ( iq_choice == 3 ) then !  make up matrix from b only
           call axisang2rot( gbdirnb , pi , matrixTwin ) ! 180 about the direction
           print*,'Matrix for Rotation Twin from Axis-Angle'
           do i = 1,3
              print"('[ ',3(1x,f8.3),' ]')",(matrixTwin(i,j),j=1,3)
           end do
           call euler( matrixTwin , 1 , 'B' , phi1(1) , capphi(1) , phi2(1) , ior , 0 )
           write(*,"('Eulers from matrixTwin: ',3(1x,f8.3))") phi1(1) , capphi(1) , phi2(1)
           phi1(1) = 0.
           capphi(1) = 0.
           phi2(1) = 0.
           aquat(1,1) = 0.
           aquat(2,1) = 0.
           aquat(3,1) = 0.
           aquat(4,1) = 1.
           aquat(1,2) = gbdirnb(1)
           aquat(2,2) = gbdirnb(2)
           aquat(3,2) = gbdirnb(3)
           aquat(4,2) = 0.
           do ijk=1,4
              qr2(ijk)=aquat(ijk,2)
           enddo
           write(*,"('quaternion = ',4(1x,f8.3),2i4)") qr2
           call q2eulB(phi1(2),capphi(2),phi2(2),qr2)
           phi1(2) = phi1(2)*180./pi
           capphi(2) = capphi(2)*180./pi
           phi2(2) = phi2(2)*180./pi
           if(phi1(2).gt.360.) phi1(2) = phi1(2)-360.
           if(phi1(2).lt.0.) phi1(2) = phi1(2)+360.
           if(phi2(2).gt.360.) phi2(2) = phi2(2)-360.
           if(phi2(2).lt.0.) phi2(2) = phi2(2)+360.
           write(*,"('check: Eulers: ',4f8.3)") phi1(2) , capphi(2) , phi2(2)
           call euler( matrixB , 2 , 'B' , phi1(2)*rad , capphi(2)*rad , phi2(2)*rad , ior , 0 )
           print*,'Matrix for Rotation Twin from Eulers'
           do i = 1,3
              print"('[ ',3(1x,f8.3),' ]')",(matrixB(i,j),j=1,3)
           end do
        end if  !  iq_choice: 1 , 2 or 3

        tmp = atan2( gbdirnb(2) , gbdirnb(1) )  !  y, then x
        long = (180. / pi) * tmp
        if ( long .lt. 0. ) long = long + 360.
        if ( gbdirnb(3) .gt. 1. ) gbdirnb(3) = 1.
        if ( gbdirnb(3) .lt. -1. ) gbdirnb(3) = -1.
        lat = 90. - (180. / pi) * acos( gbdirnb(3) )
        if ( lat .lt. 0. ) then
           long = long + 180.
           if ( long .gt. 360. ) long = long - 360.
           write( 13 ,"(2(1x,g10.4),a,i8,a,i4)") long , -1.*lat , " 10 0 0 CM -B "
        else
           write( 13 ,"(2(1x,g10.4),a,i8,a,i4)") long , lat , ' 10 0 0 CM +B'
        end if

        tmp = atan2( gbnormb(2) , gbnormb(1) )  !  y, then x
        long = (180. / pi) * tmp
        if ( long .lt. 0. ) long = long + 360.
        if ( gbnormb(3) .gt. 1. ) gbnormb(3) = 1.
        if ( gbnormb(3) .lt. -1. ) gbnormb(3) = -1.
        lat = 90. - (180. / pi) * acos( gbnormb(3) )
        if ( lat .lt. 0. ) then
           long = long + 180.
           if ( long .gt. 360. ) long = long - 360.
           write( 13 ,"(2(1x,g10.4),a,i8,a,i4)") long , -1.*lat , " 10 0 0 CM -N "
        else
           write( 13 ,"(2(1x,g10.4),a,i8,a,i4)") long , lat , ' 10 0 0 CM +N'
        end if

        tmp = atan2( matrixTwin(2,3) , matrixTwin(1,3) )  !  y, then x
        long = (180. / pi) * tmp
        if ( long .lt. 0. ) long = long + 360.
        if ( matrixTwin(3,3) .gt. 1. ) matrixTwin(3,3) = 1.
        if ( matrixTwin(3,3) .lt. -1. ) matrixTwin(3,3) = -1.
        lat = 90. - (180. / pi) * acos( matrixTwin(3,3) )
        if ( lat .lt. 0. ) then
           long = long + 180.
           if ( long .gt. 360. ) long = long - 360.
           write( 13 ,"(2(1x,g10.4),a,i8,a,i4)") long , -1.*lat , " 10 0 0 CM -z "
        else
           write( 13 ,"(2(1x,g10.4),a,i8,a,i4)") long , lat , ' 10 0 0 CM +z'
           !  see PSTEXT and Appdx F for explanation of special characters
        end if
        tmp = atan2( matrixTwin(2,2) , matrixTwin(1,2) )  !  y, then x
        long = (180. / pi) * tmp
        if ( long .lt. 0. ) long = long + 360.
        if ( matrixTwin(3,2) .gt. 1. ) matrixTwin(3,2) = 1.
        if ( matrixTwin(3,2) .lt. -1. ) matrixTwin(3,2) = -1.
        lat = 90. - (180. / pi) * acos( matrixTwin(3,2) )
        if ( lat .lt. 0. ) then
           long = long + 180.
           if ( long .gt. 360. ) long = long - 360.
           write( 13 ,"(2(1x,g10.4),a,i8,a,i4)") long , -1.*lat , " 10 0 0 CM -y "
        else
           write( 13 ,"(2(1x,g10.4),a,i8,a,i4)") long , lat , ' 10 0 0 CM +y'
        end if
        tmp = atan2( matrixTwin(2,1) , matrixTwin(1,1) )  !  y, then x
        long = (180. / pi) * tmp
        if ( long .lt. 0. ) long = long + 360.
        if ( matrixTwin(3,1) .gt. 1. ) matrixTwin(3,1) = 1.
        if ( matrixTwin(3,1) .lt. -1. ) matrixTwin(3,1) = -1.
        lat = 90. - (180. / pi) * acos( matrixTwin(3,1) )
        if ( lat .lt. 0. ) then
           long = long + 180.
           if ( long .gt. 360. ) long = long - 360.
           write( 13 ,"(2(1x,g10.4),a,i8,a,i4)") long , -1.*lat , " 10 0 0 CM -x"
        else
           write( 13 ,"(2(1x,g10.4),a,i8,a,i4)") long , lat , ' 10 0 0 CM +x'
        end if
        close(13)

     elseif( iq2 .lt. 0 .or. iq2 .gt. 8 ) then
        stop 'invalid choice'
     endif                  !  if(iq2.eq.0  ... choice of input
     !     everything should now be set up to do the calculation

  elseif ( iq == 2 ) then
     call textur2( afile , iq_data_type )
     call textur3( bfile , iq_data_type )

  end if                     !  if(iq.lt.0 ...
  !     get the orientations, symms etc.

  !..........  COMBINE the two together ...........
  if( combine.eq.1 ) then
     do ijk = 1 , 3
        qtwo( ijk , 1 ) = aquat( ijk , 1 )
        qtwo( ijk , 2 ) = -1. * aquat( ijk , 2 )
        !  reversal required to compensate for effect of comquat
     end do
     qtwo( 4 , 1 ) = aquat( 4 , 1 )
     qtwo( 4 , 2 ) = aquat( 4 , 2 )
     call comquat(qtwo,qr1)
     axis(1) = qr1(1)
     axis(2) = qr1(2)
     axis(3) = qr1(3)
     call vecnorm(axis)
!!!$            print*,'output from combination'
!!!$            print"('Quat from COMQUAT = ',4(2x,f8.3))"
!!!$     1           ,(qr1( jkl ) , jkl = 1 , 4 )
!!!$            print"('Misor. Axis= ',3(2x,f8.3))",(axis(ijk),ijk=1,3)
     do ijk = 1 , 4
        aquat( ijk , 2 ) = qr1(ijk)
        aquat( ijk , 1 ) = 0.
     end do
     aquat( 4 , 1 ) = 1.
     iq2 = 2
     phi1(1) = 0.
     capphi(1) = 0.
     phi2(1) = 0.
     call q2eulB (phi1(2) , capphi(2) , phi2(2) , qr1 )
     phi1(2)=phi1(2)*180./pi
     capphi(2)=capphi(2)*180./pi
     phi2(2)=phi2(2)*180./pi
     if(phi1(2).gt.360.) phi1(2)=phi1(2)-360.
     if(phi1(2).lt.0.) phi1(2)=phi1(2)+360.
     if(phi2(2).gt.360.) phi2(2)=phi2(2)-360.
     if(phi2(2).lt.0.) phi2(2)=phi2(2)+360.
  end if
  !  ++++++  end of input +++++++++
  !  ++++++++++++++++++++++++++++++

  open(unit=1,file='rexgbs_output.txt',status='replace')
  print*,'output file rexgbs_output.txt contains list of misors.'
  write(1,*) ngrain,' = number of components'

  open(unit=2,file='rexgbs_misors.mdf',status='replace')
  print*, 'plot rexgbs_misors.mdf with RF_misor; e.g. ...'
  print*,'/code/RF_misor rexgbs_misors.mdf 1.5 0 0 999999999 0 0'
  print*,' for a cubic case, or this for HCP:'
  print*,'/code/RF_misor rexgbs_misors.mdf 1.5 0 1 999999999 0 0'
  write(2,*) ' i j q1 q2 q3 q4 E1 E2 E3 / for plotting with RFdist '

  open(unit=3,file='rexgbs_Eulers.wts',status='replace')
  print*, 'process rexgbs_Eulers.wts with [/code/texture/]wts2pop'
  print*,'and then use Draw_stereograms to plot the discrete PFs'
  write(3,*) 'rexgbs E1 E2 E3 1.0 i j q1 q2 q3 q4 '
  write(3,*) 'Evm    F11 F12 F13 F21 F22 D23 F31 F32 F33'
  write(3,*) ' 0.   1. 0. 0.  0. 1. 0.  0. 0. 1. '
  write(3,"('Bunge:Psi  Theta   phi   ,gr.wt.')")

  open(unit=4,file='rexgbs_misors_forward.mdf',status='replace')
  print*, 'plot rexgbs_misors_forward.mdf with RF_misor'
  print*,'this contains misor-A followed by misor-B only'
  write(4,*) ' i j q1 q2 q3 q4 E1 E2 E3 / for plotting with RFdist '

  open ( unit = 14 , file='min-misor-angles.txt' , status = 'replace' )
  write ( 14 , * ) ' i spin(i)  minimum_Misor axis1 axis2 axis3 j spin(j) AngleAfterTwinApplied ax1 ax2 ax3'
  open ( unit = 15 , file='misor-angles-pairs.txt' , status = 'replace' )
  write ( 15 , * ) ' i    j  spin(i)  spin(j)   Misor-angle axis1 axis2 axis3 AngleAfterTwinApplied ax1 ax2 ax3'

  ulimit = 1
  if ( iq_extra .eq. 1 ) ulimit = numsymm
  !  used below in additional loop

  print * , ('Type anything to start analysis ...')
  read *

  !  NOTES, late on 1 x 15
  !  perhaps preserve the loop structure but, for analyzing a sequence of individual boundaries/misorientations
  !  use a stride of two in the outer loop  and confine the inner loop to i+1
  innerLoopLimit = ngrain
  innerLoopStart = 1  !  but see below for adjustment based on outer loop value
  outerLoopLimit = ngrain
  OuterLoopIncrement = 1
  !  for text input as Euler pairs, note the special treatment of the loop variables
  if ( iq_data_type == 0 ) then
     outerLoopLimit = ngrain - 1
     OuterLoopIncrement = 2
  end if
  
  if ( iq <= 1 ) outerLoopLimit = ngrain - 1
!  if ( iq_data_type == 2 ) then
  if ( iq == 2 ) then
     innerLoopLimit = bngrain
  end if
  
  do 1000  i = 1 , outerLoopLimit , OuterLoopIncrement
     minMisor = 180.  !  track which pair yields the smallest misor
     jindexMin = 0
     !if ( iq_data_type == 1 ) innerLoopStart = (i+1)
     if ( iq <= 1 ) innerLoopStart = ( i + 1 )
     p1(1)=phi1(i)
     pp(1)=capphi(i)
     p2(1)=phi2(i)
     write(*,*)
     write(1,*)
     write(1,*) '..................................................'
     write(1,*) 'Grain number ' , i
     write(*,*) '1st Grain: Euler angles: ', p1(1) , pp(1) , p2(1)
     write(1,"('1st Grain: Euler angles: ',3(1x,f8.3))") p1(1),pp(1),p2(1)
     write(1,"(4f8.3,' = 1st gr quat')") (aquat(ijk,i),ijk=1,4)
     write(*,"(4f8.3,' = 1st gr quat')") (aquat(ijk,i),ijk=1,4)

     !  use this section when applying extra symmetry to 2nd g
!!!$         do ijk=1,4
!!!$            qtwo( ijk , 1 ) = aquat( ijk , i )
!!!$         enddo

     do ijk = 1 , 4
        qtmp1(ijk)= aquat( ijk , i )
     enddo
     do 997 iex = 1 , ulimit
        if ( ulimit .gt. 1 ) then
           print*,'Inner symmetry element no. ',iex
           write(1,*) 'Inner symmetry element number ' , iex
        endif
        !  this is the extra, inserted symmetry element
        call presymm( qtmp1 , iex , qtmp2 )
        !            call postsymm( qtmp1 , iex , qtmp2 )

        !     write(1,*) ' 1st  2nd  Angle  Axis123  Quat.  Sigma  Near'
        do ijk = 1 , 4
           qtwo( ijk , 1 ) = qtmp2(ijk)
           qdisq1(ijk) = aquat( ijk , i ) !  use w disquat2
        end do

        write(*,"(4f8.3,' = result of PRESYMM')") (qtmp2( ijk ),ijk=1,4)

        !  innerLoopLimit  controlled by 1 vs 2 lists
        !  do 990, j = (i+1) , ngrain
        if ( iq_data_type == 0 ) then
           innerLoopStart = i + 1
           innerLoopLimit = innerLoopStart
        end if
        do 990 j = innerLoopStart , innerLoopLimit

           write(*,*)
           write(*,*) ' against grain number ' , j
           write(*,"(4f8.3,' = 2nd gr quat')") (aquat(ijk,j),ijk=1,4)
           write(1,"('....')")

!!!$               do ijk = 1 , 4
!!!$                  qtmp1(ijk) = aquat( ijk , j )
!!!$               enddo
!!!$               do 997 iex = 1 , ulimit
!!!$                  if ( ulimit .gt. 1 ) then
!!!$                     print*,'Inner symmetry element no. ',iex
!!!$                  endif
!!!$!  this is the extra, inserted symmetry element
!!!$                  call presymm( qtmp1 , iex , qtmp2 )
!!!$!                  call postsymm( qtmp1 , iex , qtmp2 )
!!!$                  rnorm = sqrt(qtmp2(1)**2 + qtmp2(2)**2
!!!$     $                 + qtmp2(3)**2 + qtmp2(4)**2)
!!!$                  if( rnorm .lt. 1.e-5 )
!!!$     $                 stop 'Q2ROD: too small a quaternion!'
!!!$                  if( abs( rnorm - 1.0 ) .gt. 1.e-3 ) then
!!!$                     write(*,*) 'qtmp2 required normalization'
!!!$                     do ijk = 1 , 4
!!!$                        qtmp2(ijk) = qtmp2(ijk) / rnorm
!!!$                     enddo
!!!$                  endif

           if ( iq == 0 .OR. iq == 1 ) then
              do ijk = 1 , 4
                 qtwo( ijk , 2 ) = aquat( ijk , j )
                 !                     qtwo( ijk , 2 ) = qtmp2( ijk )
                 qdisq2(ijk) = aquat( ijk , j )
              end do
           elseif ( iq == 2 ) then  !  taking data from 2nd list
              do ijk = 1 , 4
                 qtwo( ijk , 2 ) = bquat( ijk , j )
                 !                     qtwo( ijk , 2 ) = qtmp2( ijk )
                 qdisq2(ijk) = bquat( ijk , j )
              end do
           end if
           call comquat(qtwo,qr1)
           axis(1) = qr1(1)
           axis(2) = qr1(2)
           axis(3) = qr1(3)
           call vecnorm(axis)
!!!$               print"('Quat from COMQUAT = ',4(2x,f8.3))"
!!!$     1              ,(qr1( jkl ) , jkl = 1 , 4 )
!!!$               print"('Misor. Axis= ',3(2x,f8.3))",(axis(ijk),ijk=1,3)

           do ijk = 1 , 4
              qgrain1(ijk)= aquat( ijk , i )
              if ( iq == 1 ) then
                 qgrain2(ijk)= aquat( ijk , j )
              elseif ( iq == 2 ) then
                 qgrain2(ijk)= bquat( ijk , j )
              end if
           end do
           call qrvec ( qgrain1 , axis , samp1 )
!!!$               print"('Qgrain1 = ',4(2x,f8.3))"
!!!$     1              ,( qgrain1( jkl ) , jkl = 1 , 4 )
!!!$               print"('misor axis in sample 1 = ',3(2x,f8.3))"
!!!$     $              , ( samp1(ijk) , ijk = 1 , 3 )
           call qrvec ( qgrain2 , axis , samp2 )
!!!$               print"('Qgrain2 = ',4(2x,f8.3))"
!!!$     1              ,( qgrain2( jkl ) , jkl = 1 , 4 )
!!!$               print"('misor axis in sample 2 = ',3(2x,f8.3))"
!!!$     $              , ( samp2(ijk) , ijk = 1 , 3 )

           if(iq_symm.eq.0) then !  only for cubic
              !     first use Dierk's subroutine MISORI

              p1(1)=phi1(i)
              pp(1)=capphi(i)
              p2(1)=phi2(i)
              if ( iq == 1 ) then
                 p1(2)=phi1(j)
                 pp(2)=capphi(j)
                 p2(2)=phi2(j)
              elseif ( iq == 2 ) then
                 p1(2) = bphi1(j)
                 pp(2) = bcapphi(j)
                 p2(2) = bphi2(j)
              end if
!!!$                  write(*,*) '2nd Grain: Euler angles: '
!!!$     $                 , p1(2) , pp(2) , p2(2)
              !write(1,*) '2nd Grain: Euler angles: ', p1(2) , pp(2) , p2(2)
              call misori( p1 , pp , p2 , tmp , 1 , ac_min )
              write(1,*) 'MISORI: angle= ',tmp ,' ,     axis= ',ac_min
!!!$                  write(*,*) 'MISORI: angle= ',tmp
!!!$     $                 ,' ,     axis= ',ac_min
              p1(1)=phi1(i)
              pp(1)=capphi(i)
              p2(1)=phi2(i)
              if ( iq == 1 ) then
                 p1(2)=phi1(j)
                 pp(2)=capphi(j)
                 p2(2)=phi2(j)
              elseif ( iq == 2 ) then
                 p1(2) = bphi1(j)
                 pp(2) = bcapphi(j)
                 p2(2) = bphi2(j)
              end if
              call misorinv(p1,pp,p2,tmp,1,ac_min)
              !     write(1,*) 'MISORInv: angle= ',tmp,' axis= ',ac_min
              write(*,*) 'MISORInv: angle= ',tmp,' axis= ',ac_min

              call  misquat(qtwo,tmp)
              write(1,*) 'MISQUAT: angle= ',tmp
              write(*,*) 'MISQUAT: angle= ',tmp
              !     &                 ,', IDs,angle,axis,quat,sig,nr'
              !     call  misquatinv(qtwo,tmp)
              !     write(1,*) 'MISQUATINV: angle= ',tmp,' ,IDs,angle,axis,quat,sig,nr'
              !     write(*,*) 'MISQUATINV: angle= ',tmp

           endif            !  cubic only

           if ( iq_symm .eq. 0 ) then
              call disquat( qr1 , index1 , index2 , qresult , index1a , index2a )
           elseif(iq_symm.eq.1) then
              call disquat_hex (qr1,index1,index2,qresult,index1a,index2a)
           elseif(iq_symm.eq.2) then
              call disquat_orth (qr1,index1,index2,qresult,index1a,index2a)
           elseif( iq_symm == 3 ) then
              call disquat_tet (qr1,index1,index2,qresult,index1a,index2a)
           elseif( iq_symm == 4 ) then
              call disquat_tet4onM (qr1,index1,index2,qresult,index1a,index2a)
           end if
           print*
!!!$               write(*,"('from DISQUAT, qresult = ',4(1x,f8.3)/4i4)")
!!!$     1              qresult, index1 , index2 , index1a , index2a
           if(qresult(4).gt.1.0) qresult(4) = 1.0
           if(qresult(4).lt.-1.0) qresult(4) = -1.0
           angle(i,j) = 360. / pi * acos( qresult(4) )
!           angle = 360. / pi * acos( qresult(4) )
!!!$               write(*,"( 'from DISQUAT, Angle = ',f8.3)") angle(i,j)

           qresulta(i,j,1)=qresult(1)
           qresulta(i,j,2)=qresult(2)
           qresulta(i,j,3)=qresult(3)
           qresulta(i,j,4)=qresult(4)
           axis(1)=qresult(1)
           axis(2)=qresult(2)
           axis(3)=qresult(3)
           call vecnorm(axis)
           gbaxis(1,i,j) = axis(1)
           gbaxis(2,i,j) = axis(2)
           gbaxis(3,i,j) = axis(3)

           call comquatinv( qtwo , qr1 )
           QresultForSampleFrame( i , j , : ) = qr1( : )
           !  for use in last section around # 2350
           
!!!$               print"('Quat from DISQUAT = ',4(2x,f8.3))"
!!!$     1              ,(qresulta(i,j,jkl),jkl=1,4)
!!!$               print"('      Axis= ',3(2x,f8.3))",(axis(ijk),ijk=1,3)

           do ijk = 1 , 4
              qgrain1(ijk)= aquat( ijk , i )
           enddo
           call qrvec ( qgrain1 , axis , samp1 )
!!!$               print"('Qgrain1 = ',4(2x,f8.3))"
!!!$     1              ,( qgrain1( jkl ) , jkl = 1 , 4 )
!!!$               print"('misor axis in sample 1 = ',3(2x,f8.3))"
!!!$     $              , ( samp1(ijk) , ijk = 1 , 3 )
           do ijk = 1 , 4
              if ( iq == 1 ) then
                 qgrain2(ijk) = aquat( ijk , j )
              elseif ( iq == 2 ) then
                 qgrain2(ijk) = bquat( ijk , j )
              end if
           enddo
           call qrvec ( qgrain2 , axis , samp2 )
!!!$               print"('Qgrain2 = ',4(2x,f8.3))"
!!!$     1              ,( qgrain2( jkl ) , jkl = 1 , 4 )
!!!$               print"('misor axis in sample 2 = ',3(2x,f8.3))"
!!!$     $              , ( samp2(ijk) , ijk = 1 , 3 )
           call presymm( qgrain1 , index1 , qgrain3 )
           call qrvec ( qgrain3 , axis , samp1 )
!!!$               print"('Qgrain3 = ',4(2x,f8.3))"
!!!$     1              ,( qgrain3( jkl ) , jkl = 1 , 4 )
!!!$               print"('misor axis in sample 1, after symm. = '
!!!$     $                ,3(2x,f8.3))"
!!!$     $              , ( samp1(ijk) , ijk = 1 , 3 )
           call presymm( qgrain2 , index2 , qgrain4 )
           call qrvec ( qgrain4 , axis , samp1 )
!!!$               print"('Qgrain4 = ',4(2x,f8.3))"
!!!$     1              ,( qgrain4( jkl ) , jkl = 1 , 4 )
!!!$               print"('misor axis in sample 2, after symm. = '
!!!$     $                ,3(2x,f8.3))"
!!!$     $              , ( samp2(ijk) , ijk = 1 , 3 )

           do ijk=1,4
              qr2(ijk) = qresulta(i,j,ijk)
           enddo
           call q2rod(rfvector,qr2)
           write( * , "('Quat = ',4(1x,f9.4))" ) qresult
           print"('RF_vec = ',3(2x,f8.3))" , ( rfvector(ijk) , ijk = 1 , 3 )
           write(1,"('RF_vec = ',3(2x,f8.3))")  ( rfvector(ijk) , ijk = 1 , 3 )
           write( * , "('Angle = ', f8.3 , ' degrees')" ) angle(i,j)
           !write( * , "('Angle = ', f8.3 , ' degrees')" ) angle
           write( * , "('Axis = ',3(1x,f8.3))" ) ( axis(ijk) , ijk = 1 , 3 )
           ! special effort to "undo" a twin (Ag on Ni, Aug 17)
           !call ComposeQuat2( qresult , QuatReverseTwin , QuatAfterTwinApplied )
           !write( * , "('Quat after applying reverse twin = ',4(1x,f9.4))" ) ( QuatAfterTwinApplied(ijk) , ijk = 1 , 4 )
           call comquat2( qresult , QuatTwin , QuatAfterTwinApplied )
           write( * , "('Quat after applying reverse twin via comquat2 = ',4(1x,f9.4))" ) QuatAfterTwinApplied
           if( QuatAfterTwinApplied(4) > 1.0 ) QuatAfterTwinApplied(4) = 1.0d0
           if( QuatAfterTwinApplied(4) < -1.0 ) QuatAfterTwinApplied(4) = -1.0d0
           AngleAfterTwinApplied = 0.
           AngleAfterTwinApplied = 360. / pi * acos( QuatAfterTwinApplied(4) )
           write(*,"( ' Angle after undoing twin = ',f8.3)") AngleAfterTwinApplied
           AxisAfterTwinApplied(1) = QuatAfterTwinApplied(1); AxisAfterTwinApplied(2) = QuatAfterTwinApplied(2)
           AxisAfterTwinApplied(3) = QuatAfterTwinApplied(3); call vecnorm(AxisAfterTwinApplied)
           write( * , "('Axis After Twin Applied = ',3(1x,f9.4))" ) AxisAfterTwinApplied


           if ( iq_symm .eq. 1 ) then   !  hex
              MillerAxis(1) = axis(1) / 1.5
              MillerAxis(2) = ( axis(2) / sqr3 ) - ( MillerAxis(1) * 0.5 )
              MillerAxis(3) = ( -1. * MillerAxis(2) ) - MillerAxis(1)
              MillerAxis(4) = axis(3) / car
              call  vsort( MillerAxis , MillerAxisSorted , 4 )
              if ( MillerAxisSorted(4) > 1.e-5 ) then
                 factor = MillerAxisSorted(4)
              else 
                 if ( MillerAxisSorted(3) > 1.e-5 ) then
                    factor = MillerAxisSorted(3)
                 else
                    if ( MillerAxisSorted(2) > 1.e-5 ) then
                       factor = MillerAxisSorted(2)
                    else 
                       if ( MillerAxisSorted(1) > 1.e-5 ) then
                          factor = MillerAxisSorted(1)
                       end if
                    end if
                 end if
              end if
              do ijk = 1 , 4
                 MillerAxis(ijk) = MillerAxis(ijk) / factor
              end do
              write( * , "('Axis in Miller indices = ',4(1x,f8.3))" ) ( MillerAxis(ijk) , ijk = 1 , 4 )
              !))))
           else if ( iq_symm .eq. 2 ) then   !  orthorhombic
              MillerAxis(1) = axis(1)
              MillerAxis(2) = axis(2) / bar
              MillerAxis(3) = axis(3) / car
              call  vsort( MillerAxis , MillerAxisSorted , 3 )
              if ( MillerAxisSorted(3) > 1.e-5 ) then
                 factor = MillerAxisSorted(3)
              else
                 if ( MillerAxisSorted(2) > 1.e-5 ) then
                    factor = MillerAxisSorted(2)
                 else 
                    if ( MillerAxisSorted(1) > 1.e-5 ) then
                       factor = MillerAxisSorted(1)
                    end if
                 end if
              end if
              do ijk = 1 , 3
                 MillerAxis(ijk) = MillerAxis(ijk) / factor
              end do
              write( * , "('Axis in Miller indices = ',3(1x,f8.3))" ) ( MillerAxis(ijk) , ijk = 1 , 3 )
              !))))
           else if ( iq_symm >=3 .and. iq_symm <=4 ) then   !  tetragonal
              MillerAxis(1) = axis(1)
              MillerAxis(2) = axis(2)
              MillerAxis(3) = axis(3) / car
              call  vsort( MillerAxis , MillerAxisSorted , 3 )
              if ( MillerAxisSorted(3) > 1.e-5 ) then
                 factor = MillerAxisSorted(3)
              else
                 if ( MillerAxisSorted(2) > 1.e-5 ) then
                    factor = MillerAxisSorted(2)
                 else 
                    if ( MillerAxisSorted(1) > 1.e-5 ) then
                       factor = MillerAxisSorted(1)
                    end if
                 end if
              end if
              do ijk = 1 , 3
                 MillerAxis(ijk) = MillerAxis(ijk) / factor
              end do
              write( * , "('Axis in Miller indices = ',3(1x,f8.3))" ) ( MillerAxis(ijk) , ijk = 1 , 3 )
           end if
           !  if you want to check the result (for cubics) uncomment this block
           if( iq_symm == 0 ) then
              call disquat2(qdisq1,qdisq2,index1 ,index2,qresult,index1a,index2a)
              if(qresult(4).gt.1.0) qresult(4)=1.0
              if(qresult(4).lt.-1.0) qresult(4)=-1.0
              angle(i,j)=360./pi*acos(qresult(4))
              !angle = 360./pi*acos(qresult(4))
!!!$                  write(*,"( 'from DISQUAT2, Angle = ',f8.3)")
!!!$     $                 angle(i,j)
!!!$               write(*,"('from DISQUAT2, qresult = ',4(1x,f8.3),4i4)")
!!!$     1              qresult,index1,index2,index1a,index2a
           endif

           if( iq_symm == 0 ) then  ! only if CUBIC
              lsig(i,j)=0
              rnear(i,j)=99.
              if(angle(i,j).lt.15.0) then
              !if( angle .lt. 15.0 ) then
                 lsig(i,j)=1
                 rnear(i,j)=angle(i,j)
                 !rnear(i,j) = angle
                 !     call it sigma=1 if low angle bdy
              end if

              thetmin = 66.0
              if( lsig(i,j) == 0 ) then
                 do 150, ijk = 1 , numcsl
                    brandon=15.0/sqrt(float(nsig(ijk)))
                    call cslquat(qresult,ijk,theta,qr2)
                    crit=theta/brandon
                    !      write(*,"('sig,qr2,theta = ',i4,8f8.3)") nsig(ijk),qr2,theta,crit
                    if(crit.lt.thetmin) then
                       !      write(*,"('sig,qr2,theta = ',i4,8f8.3)") nsig(ijk),qr2,theta
                       thetmin = crit
                       lsig(i,j) = nsig(ijk)
                       rnear(i,j) = crit
                    end if
150              end do
              end if  !   if(lsig(i,j).eq.0)
              print*
              print"('indices: ',2i4,'   angle= ',f8.3)",i,j,angle(i,j)
              !print"('nearest Sigma, distance = 'i4,f8.3)" ,lsig(i,j),rnear(i,j)
              !     write(1,"(' IDs, angle, axis, quat, sigma, nearness')")
              write(1,"('symm op or grain indices: ',2(1x,i4),'   angle= ',1x,f8.3)") i,j,angle(i,j)
              if ( abs(angle(i,j)-60.) < 0.1 ) write(1,"('**Angle is close to 60 **')") 
              !write(1,"('Angle = ',f8.3)") angle
              write(1,"('Axis = ',3(1x,f8.3))") (axis(ijk),ijk=1,3)
              write(1,"('quaternion = ',4(1x,f8.3))") (qresulta(i,j,jkl),jkl=1,4)
              write(1,"('sigma value, near = ',i5,(1x,f8.3))") lsig(i,j),rnear(i,j)

           else             !  non cubic, leave out Sigma
              write(1,*) " i j Angle(degr) Axis(3) Quaternion(4)"
              write(1,890) i , j , angle(i,j) , (axis(ijk),ijk=1,3) &
              !write(1,890) i , j , angle , (axis(ijk),ijk=1,3) &
                   &, (qresulta(i,j,jkl),jkl=1,4)
           end if            !  cubic?

           if ( angle(i,j) < minMisor ) then
           !if ( angle < minMisor ) then
              minMisor = angle(i,j)
              jindexMin = j
           end if
           !  check for pair that gives the smallest misor angle
           !write ( 15 , * ) 'i j  aspin(i)  bspin(j) Misor-angle Axis: '&
           if ( iq_symm == 0 ) then
              write ( 15 , * ) &
                   !&, i , j , aspin(i) , bspin(j), angle
                   &, i , j , aspin(i) , bspin(j), angle(i,j), (axis(ijk),ijk=1,3), AngleAfterTwinApplied, AxisAfterTwinApplied
           else
              write ( 15 , * ) i , j , aspin(i) , bspin(j), angle(i,j), (axis(ijk),ijk=1,3)
           end if

           do jkl = 1 , 4
              qtmp(jkl) = qresulta(i,j,jkl)
           enddo
           call q2eulB( d1 , d2 , d3 , qtmp )  !  d3 was d2, corrected 16 viii 11
           write(2,"(2i6,8(1x,f10.4))") i,j ,(qtmp(jkl),jkl=1,4),d1*rad,d2*rad,d3*rad
           !     write out result for RFdist input

           goto 2100
           !     jump around reverse order calculation if above not commented
           !     *****************  now for the reverse order!

           call comquatinv(qtwo,qr1)
           call disquat(qr1,index1,index2,qresult,index1a,index2a)
           write(*,"('qresult= ',4f8.3,2i4)") qresult,index1a,index2a
           if(qresult(4).gt.1.0) qresult(4)=1.0
           if(qresult(4).lt.-1.0) qresult(4)=-1.0
           angle(i,j)=360./pi*acos(qresult(4))
           !angle = 360./pi*acos(qresult(4))
           write(*,"('Angle= ',f8.3)") angle(i,j)
           !write(*,"('Angle= ',f8.3)") angle
           lsig(i,j)=0
           rnear(i,j)=99.
           qresulta(i,j,1)=qresult(1)
           qresulta(i,j,2)=qresult(2)
           qresulta(i,j,3)=qresult(3)
           qresulta(i,j,4)=qresult(4)
           if(angle(i,j).lt.15.0) then
              lsig(i,j)=1
              rnear(i,j)=angle(i,j)
              !rnear(i,j)=angle
              !     call it sigma=1 if low angle bdy
           endif
           axis(1)=qresult(1)
           axis(2)=qresult(2)
           axis(3)=qresult(3)
           call vecnorm(axis)  !  make unit length
           gbaxis(1,i,j)=axis(1)
           gbaxis(2,i,j)=axis(2)
           gbaxis(3,i,j)=axis(3)

           thetmin=66.0
           if(lsig(i,j).eq.0) then
              do 250, ijk=1,numcsl
                 brandon=15.0/sqrt(float(nsig(ijk)))
                 call cslquat(qresult,ijk,theta,qr2)
                 crit=theta/brandon
                 if(crit.lt.thetmin) then
                    !     write(*,"('qr2,theta= ',8f8.3)") qr2,theta
                    thetmin=crit
                    lsig(i,j)=nsig(ijk)
                    rnear(i,j)=crit
                 endif
250           end do
           endif

           write(*,890) i,j,angle(i,j),(axis(ijk),ijk=1,3) &
           !write(*,890) i,j,angle , (axis(ijk),ijk=1,3) &
                &,(qresulta(i,j,jkl),jkl=1,4),lsig(i,j),rnear(i,j)
           write(1,890) i,j,angle(i,j),(axis(ijk),ijk=1,3) &
           !write(1,890) i,j,angle , (axis(ijk),ijk=1,3) &
                &,(qresulta(i,j,jkl),jkl=1,4),lsig(i,j),rnear(i,j)
890        format(2i4,8f8.3,i4,f8.3)

           !     *****************  
2100       continue         !  jumps around reverse order calculation to here

990     end do  !  INNER LOOP

997  end do   !  apply extra symmetry to 1st quat

     print * , ' for i= ', i ,' the smallest misor = ' , minMisor,' for number ', jindexMin
     ! 997     enddo   !  apply extra symmetry to 2nd quat
     write ( 14 , * ) i , aspin(i) , minMisor, (axis(ijk),ijk=1,3) , jindexMin , bspin(jindexMin)  &
          & , AngleAfterTwinApplied, AxisAfterTwinApplied

1000 end do   !  OUTER LOOP

  write(1,*)
  write(1,*) 'Next section intended as input for RFplot==rfp'
  write(1,*) 'which you can find in /dihedral/.'
  write(1,*)

  write(*,*) 'Last section of rexgbs.out intended as '
  write(*,*) 'input for RFplot==rfp, '
  print*,'  which you can find in /dihedral/.'
  print*,'on the screen, R1,R2,R3; in the file R3,R2,R1 '

  print*
  print*,'Output also available as rexgbs_misors.mdf, '
  print*,' which you can plot with /code/MDF/RF_misor.'
  print*

  innerLoopLimit = ngrain
  innerLoopStart = 1
  outerLoopLimit = ngrain
  OuterLoopIncrement = 1
  !  for text input as Euler pairs, note the special treatment of the loop variables
  if ( iq_data_type == 0 ) then
     outerLoopLimit = ngrain - 1
     OuterLoopIncrement = 2
  end if
  if ( iq <= 1 ) outerLoopLimit = ngrain - 1
!  if ( iq_data_type == 2 ) then
  if ( iq == 2 ) then
     innerLoopLimit = bngrain
  end if

  !  redefine the vectors to make sure
  vec100(1) = 0.; vec100(2) = 0.; vec100(3) = 1.
  vec110(1) = 0.; vec110(2) = 1.; vec110(3) = 1.
  vec111(1) = 1.; vec111(2) = 1.; vec111(3) = 1.
  call vecnorm(vec100) ; call vecnorm(vec110) ; call vecnorm(vec111)
  
  open( unit = 17 , file = 'OutputForAnalysis.txt' , status = 'replace' )
  write( 17 , * ) ' i j sigma weight angle1 rf1 rf2 rf3 near100 near110 near111' &
       & , ' MisorSampleX MisorSampleY MisorSampleZ'
  !QtwoForSampleFrame( 4 , 2 ) , QresultForSampleFrame(4) , AxisForSampleFrame(3)
  do 2000, ij = 1 , outerLoopLimit , OuterLoopIncrement
     if ( iq_data_type == 0 ) then
        innerLoopStart = ij + 1
        innerLoopLimit = innerLoopStart
     end if
     do jk = innerLoopLimit , innerLoopLimit
        rfvector(1) = qresulta(ij,jk,1)/qresulta(ij,jk,4)
        rfvector(2) = qresulta(ij,jk,2)/qresulta(ij,jk,4)
        rfvector(3) = qresulta(ij,jk,3)/qresulta(ij,jk,4)
        !            print*,'qresulta: ',(qresulta(i,j,jkl),jkl=1,4)
        call comquatinv( QtwoForSampleFrame , QresultForSampleFrame )
        AxisForSampleFrame( 1 ) = QresultForSampleFrame( ij , jk , 1 ) / QresultForSampleFrame( ij , jk , 4 )
        AxisForSampleFrame( 2 ) = QresultForSampleFrame( ij , jk , 2 ) / QresultForSampleFrame( ij , jk , 4 )
        AxisForSampleFrame( 3 ) = QresultForSampleFrame( ij , jk , 3 ) / QresultForSampleFrame( ij , jk , 4 )
        call donorm( AxisForSampleFrame )
        !  tested this 29 i 16 with these combos hkl uvw
        !  111 -110  & 111 -211, which gave 30/111 in xtal, 001 in sample
        !  -211 111  & -110 111, which gave 30/111 in xtal, 100 in sample
        !  these were the expected results
        !print * , 'Axis from Sample frame misorientation = ', AxisForSampleFrame
        write(*,892) ij,jk,lsig(ij,jk)  &
             , angle(ij,jk)  &
             !,angle  &
             , qresulta(ij,jk,1)/qresulta(ij,jk,4)  &
             , qresulta(ij,jk,2)/qresulta(ij,jk,4)  &
             , qresulta(ij,jk,3)/qresulta(ij,jk,4)  &
             , grwt(ij)
        if(qresulta(ij,jk,4).ne.0.) then
           write(1,892) ij,jk,lsig(ij,jk)   &
                , angle(ij,jk)               &
                !,angle               &
                , qresulta(ij,jk,3)/qresulta(ij,jk,4)  &
                , qresulta(ij,jk,2)/qresulta(ij,jk,4)  &
                , qresulta(ij,jk,1)/qresulta(ij,jk,4)  &
                , grwt(ij)
           !    ,1.0,0.,0.,1.0,0.,0.
892        format( 3(1x,i4),1x,4(1x,f9.3),1x,3(1x,f9.3),1x,4(f12.5) )
893        format( 3(1x,i4),1x,2(1x,f9.3),1x,3(1x,f9.3),1x,3(f12.5),1x,3(f12.5) )
           call vecnorm( rfvector )
           rtmp100 = (scalar( rfvector , vec100 ) )**4
           rtmp110 = (scalar( rfvector , vec110 ) )**4
           rtmp111 = (scalar( rfvector , vec111 ) )**4
           write( 17 , 893 ) ij , jk , lsig(ij,jk) , grwt(ij)   &
                , angle(ij,jk)               &
                !,angle               &
                , qresulta(ij,jk,1)/qresulta(ij,jk,4)  &
                , qresulta(ij,jk,2)/qresulta(ij,jk,4)  &
                , qresulta(ij,jk,3)/qresulta(ij,jk,4)  &
                , rtmp100 , rtmp110 , rtmp111 &
                ,  AxisForSampleFrame
        end if
        !  converts the quaternion to Rodrigues vector
     end do
2000 end do
  close ( unit = 17 )
  write(*,*) 'Rotate a vector from one lattice to the other [yes=1]?'
  read(*,*) iq
  if(iq.eq.1) then
     write(*,*) 'Enter the vector to be rotated (transformed)'
     read(*,*) tmpvec
     call donorm(tmpvec)
     call rfrvecR ( rfvector , tmpvec , tmpout )
     write(*,"('input vector normalized',3(1x,f8.4))") tmpvec
     write(*,"('output vector from rfrvecR: ',3(1x,f8.4))") tmpout
     print * , 'Repeat with QRVEC, INVQRVEC ...'
     call qrvec ( qresult , tmpvec , tmpout )
     write(*,"('output from QRVEC: ',3(1x,f8.4))") tmpout
     call invqrvec ( qresult , tmpvec , tmpout )
     write(*,"('output from INVQRVEC: ',3(1x,f8.4))") tmpout
  end if

  close(1)
  close(2)

  call exit(0)

end program rexgbs

! .............................................................
subroutine uvwhkl2quat( gbnorma , gbdirna , quatn , Bunge , a )

  implicit none

  double precision :: scalar
  double precision :: gbnorma(3),gbdirna(3),quatn(4),Bunge(3)
  double precision :: degrad , rtmp , t1(3) , a(3,3) , d1 , d2 , d3
  integer :: i,iopt,kerr,ior,j
  character :: nomen*1
  double precision :: tmp
  double precision :: axis(3)
  double precision, parameter :: pi  = 4. * atan ( 1.0 ) 

  logical verbose
  common /general/ verbose

  !  CODE::
  print*
  degrad = 180./pi
  call vecnorm(gbnorma)
  call vecnorm(gbdirna)
  rtmp = scalar(gbnorma,gbdirna)
  if(abs(rtmp).gt.1e-4) then
     stop 'uvwhkl2quat: sorry, uvw is not perp. to hkl'
  endif
  call vecpro2(gbnorma,gbdirna,t1)
  call vecnorm(t1)
  do 20, i = 1,3
     a(3,i) = gbnorma(i)  ! 3rd ROW
     a(2,i) = t1(i)       ! 2nd ROW
     a(1,i) = gbdirna(i)  ! 1st ROW
20 enddo
  !  gives the transpose of what you find in the books
  !  because this is what Euler needs as input
  iopt = 1  !  to get angles from matrix
  nomen = 'B'
  kerr = 0
  ior = 1                   ! grain number
  call euler(a,iopt,nomen,Bunge(1),Bunge(2),Bunge(3),ior,kerr)
  !  note: should use Changsoo's more reliable code
  if(Bunge(1).lt.0.) Bunge(1) = Bunge(1) + 360.
  if(Bunge(2).lt.0.) Bunge(2) = Bunge(2) + 360.
  if(Bunge(3).lt.0.) Bunge(3) = Bunge(3) + 360.

  !  now we transpose A
  rtmp = a(1,3)
  a(1,3) = a(3,1)
  a(3,1) = rtmp
  rtmp = a(2,3)
  a(2,3) = a(3,2)
  a(3,2) = rtmp
  rtmp = a(1,2)
  a(1,2) = a(2,1)
  a(2,1) = rtmp

  rtmp = 1. + a(1,1) + a(2,2) + a(3,3)
  if ( rtmp .lt. 0. ) rtmp = 0.
  if ( abs(rtmp) .lt. 1.e-4 ) then
     !  these lines deal with the case that the angle = 180 degr
     axis(1) = sqrt( ( a(1,1) + 1.) * 0.5 )
     axis(2) = sqrt( ( a(2,2) + 1.) * 0.5 )
     axis(3) = sqrt( ( a(3,3) + 1.) * 0.5 )
     !quatn(4) = tmp * 0.25
     quatn(4) = rtmp * 0.25  !  CHECK THIS!
     quatn(3) = axis(3)
     quatn(2) = axis(2)
     quatn(1) = axis(1)
  else   !  normal situation
     tmp = 2. * sqrt ( rtmp )
     quatn(4) = tmp * 0.25
     quatn(3) = ( a(1,2) - a(2,1) ) / tmp
     quatn(2) = ( a(3,1) - a(1,3) ) / tmp
     quatn(1) = ( a(2,3) - a(3,2) ) / tmp
     !  this implements Morawiec's formula to directly from Eulers to quat
     !  see the spreadsheet for 27-750, Hwk on orientations
  endif

  if (verbose) then
     print*, 'HKLUVW2QUAT:  Bunge angles: ',Bunge(1),Bunge(2),Bunge(3)
     print*,'HKLUVW2QUAT:  Quat : ',(quatn(i),i=1,4)
  endif
  if (verbose) then
     print*,'HKLUWV2QUAT: Orientation matrix:'
     do i = 1,3
        print"('[ ',3(1x,f8.3),' ]')",(a(i,j),j=1,3)
     enddo
  endif

  return
end subroutine uvwhkl2quat

! ________________________________________

subroutine vecnorm(vec)
  double precision vec(3),rnorm
  rnorm=0.
  do 10, i=1,3
     rnorm=rnorm+vec(i)**2
10 end do
  if(rnorm.le.0.0) return
  rnorm=sqrt(rnorm)
  do 20, i=1,3
     vec(i)=vec(i)/rnorm
20 end do
  return
end subroutine vecnorm

! ________________________________________

subroutine quatB(p1,p,p2,q)
  !  convert Bunge Euler angles (radians) to quaternion;
  !  direct conversion from angles, see Altmann's book
  !  the rotation is a vector transformation (active rotation)

  double precision :: p1,p,p2,q(4)
  double precision :: c1,c,c2,s1,s,s2,d(3,3),tmp(4),sin2,cos2,rtmp,rnorm

  PI=3.14159265
  !  form cosine, sine of Phi, and sum & diff of phi1, phi2
  S=DSIN(0.5d0*P)
  C=DCOS(0.5d0*P)
  S1=DSIN(0.5d0*(P1-P2))
  C1=DCOS(0.5d0*(P1-P2))
  S2=DSIN(0.5d0*(P1+P2))
  C2=Dcos(0.5d0*(P1+P2))
  !      write(*,*) 's1,c1,s,c,s2,c2'
  !      write(*,*) s1,c1
  !      write(*,*) s,c
  !      write(*,*) s2,c2
  q(1) = real ( s*c1 )
  q(2) = real ( s*s1 )
  q(3) = real ( c*s2 )
  q(4) = real ( c*c2 )
  return
end subroutine quatB
!  ____________________________

subroutine rod2q(d,quat)
  implicit none

  double precision d(3),quat(4),rtmp,rnorm,rnorm2
  !  converts Rodrigues vector to a quaternion

  !  CODE::
  rtmp=d(1)**2+d(2)**2+d(3)**2
  !	rnorm=sqrt(rtmp)
  rnorm2=sqrt(1.+rtmp)
  quat(4)=1./rnorm2
  quat(3)=d(3)/rnorm2
  quat(2)=d(2)/rnorm2
  quat(1)=d(1)/rnorm2
  !  cute - this works for all rotation angles!
  return
end subroutine rod2q

!  ____________________________

subroutine q2rod ( d , quat )

  implicit none

  integer i
  double precision d(3),quat(4),rtmp,rnorm,rnorm2,stheta,ttheta
  double precision large
  data large / 1.e33 /
  !  converts a quaternion to a Rodrigues vector

  !  CODE::
  rnorm = sqrt(quat(1)**2+quat(2)**2+quat(3)**2+quat(4)**2)
  !if( rnorm .lt. 1.e-5 ) stop 'Q2ROD: too small a quaternion!'
  if( rnorm < 1.e-5 ) then
     d(1) = 0.
     d(2) = 0.
     d(3) = 0.
     return
  end if
  if( abs( rnorm - 1.0 ) .gt. 1.e-3 ) then
     write(*,*) 'Q2ROD: quat required normalization'
     do i=1,4
        quat(i)=quat(i)/rnorm
     enddo
  endif

  if( abs(quat(4)) .gt. 1.e-2 ) then
     d(1) = quat(1)/quat(4)
     d(2) = quat(2)/quat(4)
     d(3) = quat(3)/quat(4)
  else
     !  we calculate sine(theta/2) from quat1-3
     !  and scale by tangent/sine in order to make
     !  the system deal with large numbers!

     rtmp = sqrt(quat(1)**2+quat(2)**2+quat(3)**2)
     !  this will be close to unity
     stheta = asin(rtmp)
     if ( abs(stheta) .gt. 1.e-5 ) then
        ttheta = tan(stheta) ! and this will be large
        d(1) = quat(1)*ttheta/stheta
        d(2) = quat(2)*ttheta/stheta
        d(3) = quat(3)*ttheta/stheta
     else
        d(1) = large
        d(2) = large
        d(3) = large
     endif
  endif

  return
end subroutine q2rod

! _____________________________

subroutine ComposeQuat2(qq1,qq2,qresult)
  double precision qq1(4),qq2(4),qresult(4)

  !  algorithm for forming resultant quaternion as Q1 x Q2
  !  where Q2 follows Q1, i.e. is the second rotation;
  !  note change of signs to get inverse of second orientation

  qresult(1) = qq1(1)*qq2(4) + qq1(4)*qq2(1) + qq1(2)*qq2(3) - qq1(3)*qq2(2)
  qresult(2) = qq1(2)*qq2(4) + qq1(4)*qq2(2) + qq1(3)*qq2(1) - qq1(1)*qq2(3)
  qresult(3) = qq1(3)*qq2(4) + qq1(4)*qq2(3) + qq1(1)*qq2(2) - qq1(2)*qq2(1)
  qresult(4) = qq1(4)*qq2(4) - qq1(1)*qq2(1) - qq1(2)*qq2(2) - qq1(3)*qq2(3)

  !	write(*,*) 'qresult ',qresult

  return
end subroutine ComposeQuat2

! _____________________________

subroutine comquat2(qq1,qq2,qresult)
  double precision qq1(4),qq2(4),qresult(4)

  !       PI=3.14159265

  !  algorithm for forming resultant quaternion as Q1 x Q2^-1
  !  where Q2 follows Q1, i.e. is the second rotation;
  !  note change of signs to get inverse of second orientation

  qresult(1)=qq1(1)*qq2(4)-qq1(4)*qq2(1) +qq1(2)*qq2(3)-qq1(3)*qq2(2)
  qresult(2)=qq1(2)*qq2(4)-qq1(4)*qq2(2) +qq1(3)*qq2(1)-qq1(1)*qq2(3)
  qresult(3)=qq1(3)*qq2(4)-qq1(4)*qq2(3) +qq1(1)*qq2(2)-qq1(2)*qq2(1)
  qresult(4)=qq1(4)*qq2(4)+qq1(1)*qq2(1) +qq1(2)*qq2(2)+qq1(3)*qq2(3)

  !	write(*,*) 'qresult ',qresult

  return
end subroutine comquat2

! _____________________________

subroutine comquat(qq,qresult)
  double precision qq(4,2),qresult(4)

  !       PI=3.14159265

  !  algorithm for forming resultant quaternion as Q1 x Q2^-1
  !  where Q2 follows Q1, i.e. is the second rotation;
  !  note change of signs to get inverse of second orientation

  qresult(1)=qq(1,1)*qq(4,2)-qq(4,1)*qq(1,2) +qq(2,1)*qq(3,2)-qq(3,1)*qq(2,2)
  qresult(2)=qq(2,1)*qq(4,2)-qq(4,1)*qq(2,2) +qq(3,1)*qq(1,2)-qq(1,1)*qq(3,2)
  qresult(3)=qq(3,1)*qq(4,2)-qq(4,1)*qq(3,2) +qq(1,1)*qq(2,2)-qq(2,1)*qq(1,2)
  qresult(4)=qq(4,1)*qq(4,2)+qq(1,1)*qq(1,2) +qq(2,1)*qq(2,2)+qq(3,1)*qq(3,2)

  write(*,*) 'COMQUAT: qresult ',qresult

  return
end subroutine comquat

! _____________________________

subroutine comquatinv(qq,qresult)
  double precision qq(4,2),qresult(4)

  !       PI=3.14159265

  !  algorithm for forming resultant quaternion as Q1^1 x Q2
  !  where Q2 follows Q1, i.e. is the second rotation;
  !  note change of signs to get inverse of second orientation

  qresult(1)=(-1.*qq(1,1)*qq(4,2))+qq(4,1)*qq(1,2) +qq(2,1)*qq(3,2)-qq(3,1)*qq(2,2)
  qresult(2)=(-1.*qq(2,1)*qq(4,2))+qq(4,1)*qq(2,2) +qq(3,1)*qq(1,2)-qq(1,1)*qq(3,2)
  qresult(3)=(-1.*qq(3,1)*qq(4,2))+qq(4,1)*qq(3,2) +qq(1,1)*qq(2,2)-qq(2,1)*qq(1,2)
  qresult(4)=qq(4,1)*qq(4,2)+qq(1,1)*qq(1,2) +qq(2,1)*qq(2,2)+qq(3,1)*qq(3,2)

  !	write(*,*) 'qresult ',qresult

  return
end subroutine comquatinv

! _____________________________

subroutine presymm(qq,lindex,qresult)
  use var
  implicit none

  double precision qq(4),qresult(4)
  integer lindex

  !     include 'common.f'

!!$  integer numsymm,numsamp
!!$  double precision quatsymm(4,48),quatsamp(4,4)
!!$  common/a1/ quatsymm,numsymm,quatsamp,numsamp

  !  CODE::

!!!$c  algorithm for forming resultant quaternion
!!!$c  from applying a symmetry operator, QUATSYMM(n,lindex)
!!!$c  to the first quaternion/rotation
!!!$
!!!$c  If the symm oper == O, then
!!!$c  Qresult = (O x QQ) where QQ= Q1 x Q2^-1
!!!$c  i.e. Qresult = (O x Q1) x Q2^-1
!!!$c  note change of signs to get inverse of first orientation
!!!$
!!!$	if(lindex.gt.numsymm) stop 'error in presymm, lindex>numsymm'
!!!$	if(lindex.lt.1) stop 'error in presymm, lindex<1'
!!!$	qresult(1) = quatsymm(4,lindex)*qq(1) + quatsymm(1,lindex)*qq(4)
!!!$     &	+ quatsymm(3,lindex)*qq(2) - quatsymm(2,lindex)*qq(3)
!!!$	qresult(2) = quatsymm(4,lindex)*qq(2) + quatsymm(2,lindex)*qq(4)
!!!$     &	+ quatsymm(1,lindex)*qq(3) - quatsymm(3,lindex)*qq(1)
!!!$	qresult(3) = quatsymm(4,lindex)*qq(3) + quatsymm(3,lindex)*qq(4)
!!!$     &	+ quatsymm(2,lindex)*qq(1) - quatsymm(1,lindex)*qq(2)
!!!$	qresult(4) = quatsymm(4,lindex)*qq(4) - quatsymm(1,lindex)*qq(1)
!!!$     &	- quatsymm(2,lindex)*qq(2) - quatsymm(3,lindex)*qq(3)
  ! old version

  !  If the symm oper == O, and the input quaternion, QQ, is an active
  !  rotation, then Qresult = (O x QQ)
  !  the "PRE" in the name refers to writing the operator before the
  !  rotation/orientation in vector/tensor notation:  Q' = O x Q
  !  or in conventional quaternion notation, Qresult = QQ Ä Qsymm

  !  Thus, for active rotations (standard definition of orientation in
  !  mechanics) this is suitable for applying SAMPLE symmetry

  if(lindex.gt.numsymm) stop 'error in PRESYMM, lindex>numsamp'
  if(lindex.lt.1) stop 'error in presymm, lindex<1'
  qresult(1)=qq(1)*quatsymm(4,lindex) + qq(4)*quatsymm(1,lindex) - qq(2)*quatsymm(3,lindex) + qq(3)*quatsymm(2,lindex)
  qresult(2)=qq(2)*quatsymm(4,lindex) + qq(4)*quatsymm(2,lindex) - qq(3)*quatsymm(1,lindex) + qq(1)*quatsymm(3,lindex)
  qresult(3)=qq(3)*quatsymm(4,lindex) + qq(4)*quatsymm(3,lindex) - qq(1)*quatsymm(2,lindex) + qq(2)*quatsymm(1,lindex)
  qresult(4)=qq(4)*quatsymm(4,lindex) - qq(1)*quatsymm(1,lindex) - qq(2)*quatsymm(2,lindex) - qq(3)*quatsymm(3,lindex)
  !  taken from MILL2EUL

  !	write(*,*) 'qresult ',qresult

  return
end subroutine presymm
! _____________________________

subroutine presamp(qq,lindex,qresult)
  use var
  implicit none
  
  double precision qq(4),qresult(4)
  integer lindex

!!$  integer numsymm,numsamp
!!$  double precision quatsymm(4,48),quatsamp(4,4)
!!$  common/a1/ quatsymm,numsymm,quatsamp,numsamp

  !  CODE::

!!!$c  algorithm for forming resultant quaternion
!!!$c  from applying a symmetry operator, QUATSYMM(n,lindex)
!!!$c  to the first quaternion/rotation
!!!$
!!!$c  If the symm oper == O, then
!!!$c  Qresult = (O x QQ) where QQ= Q1 x Q2^-1
!!!$c  i.e. Qresult = (O x Q1) x Q2^-1
!!!$c  note change of signs to get inverse of first orientation
!!!$
!!!$	if( lindex .gt. numsamp) stop 'error in presymm, lindex>numsamp'
!!!$	if( lindex .lt. 1 ) stop 'error in presamp, lindex<1'
!!!$	qresult(1) = quatsamp(4,lindex)*qq(1) + quatsamp(1,lindex)*qq(4)
!!!$     &	+ quatsamp(3,lindex)*qq(2) - quatsamp(2,lindex)*qq(3)
!!!$	qresult(2) = quatsamp(4,lindex)*qq(2) + quatsamp(2,lindex)*qq(4)
!!!$     &	+ quatsamp(1,lindex)*qq(3) - quatsamp(3,lindex)*qq(1)
!!!$	qresult(3) = quatsamp(4,lindex)*qq(3) + quatsamp(3,lindex)*qq(4)
!!!$     &	+ quatsamp(2,lindex)*qq(1) - quatsamp(1,lindex)*qq(2)
!!!$	qresult(4) = quatsamp(4,lindex)*qq(4) - quatsamp(1,lindex)*qq(1)
!!!$     &	- quatsamp(2,lindex)*qq(2) - quatsamp(3,lindex)*qq(3)
  ! old version

  !  If the symm oper == O, and the input quaternion, QQ, is an active
  !  rotation, then Qresult = (O x QQ)
  !  the "PRE" in the name refers to writing the operator before the
  !  rotation/orientation in vector/tensor notation:  Q' = O x Q
  !  or in conventional quaternion notation, Qresult = QQ Ä Qsymm

  !  Thus, for active rotations (standard definition of orientation in
  !  mechanics) this is suitable for applying SAMPLE symmetry

  if(lindex.gt.numsamp) stop 'error in PRESAMP, lindex>numsamp'
  if(lindex.lt.1) stop 'error in presamp, lindex<1'
  qresult(1)=qq(1)*quatsamp(4,lindex) + qq(4)*quatsamp(1,lindex) - qq(2)*quatsamp(3,lindex) + qq(3)*quatsamp(2,lindex)
  qresult(2)=qq(2)*quatsamp(4,lindex) + qq(4)*quatsamp(2,lindex) - qq(3)*quatsamp(1,lindex) + qq(1)*quatsamp(3,lindex)
  qresult(3)=qq(3)*quatsamp(4,lindex) + qq(4)*quatsamp(3,lindex) - qq(1)*quatsamp(2,lindex) + qq(2)*quatsamp(1,lindex)
  qresult(4)=qq(4)*quatsamp(4,lindex) - qq(1)*quatsamp(1,lindex) - qq(2)*quatsamp(2,lindex) - qq(3)*quatsamp(3,lindex)
  !  taken from MILL2EUL

  !	write(*,*) 'qresult ',qresult

  return
end subroutine presamp

! _____________________________

subroutine postsymm(qq,lindex,qresult)
  use var
  implicit none
  
  double precision qq(4),qresult(4)
  integer lindex

!!$  integer numsymm,numsamp
!!$  double precision quatsymm(4,48),quatsamp(4,4)
!!$  common/a1/ quatsymm,numsymm,quatsamp,numsamp


!!!$c  algorithm for forming resultant quaternion
!!!$c  from applying a symmetry operator, QUATSYMM(n,lindex)
!!!$c  to the second quaternion
!!!$
!!!$c  If the symm oper == O, then
!!!$c  Qresult = (QQ x O^-1) where QQ= Q1 x Q2^-1
!!!$c  i.e. Qresult = Q1 x (O x Q2)^-1
!!!$c  note change of signs to get inverse of first orientation
!!!$
!!!$	if(lindex.gt.numsymm) stop 'error in postsymm, lindex>numsymm'
!!!$	if(lindex.lt.1) stop 'error in postsymm, lindex<1'
!!!$	qresult(1)=qq(1)*quatsymm(4,lindex) - qq(4)*quatsymm(1,lindex) + qq(2)*quatsymm(3,lindex) - qq(3)*quatsymm(2,lindex)
!!!$	qresult(2)=qq(2)*quatsymm(4,lindex) - qq(4)*quatsymm(2,lindex) + qq(3)*quatsymm(1,lindex) - qq(1)*quatsymm(3,lindex)
!!!$	qresult(3)=qq(3)*quatsymm(4,lindex) - qq(4)*quatsymm(3,lindex) + qq(1)*quatsymm(2,lindex) - qq(2)*quatsymm(1,lindex)
!!!$	qresult(4)=qq(4)*quatsymm(4,lindex) + qq(1)*quatsymm(1,lindex) + qq(2)*quatsymm(2,lindex) + qq(3)*quatsymm(3,lindex)
  !  old version of POSTSYMM

  !  If the symm oper == O, and the input quaternion, QQ, is an active
  !  rotation, then Qresult = (QQ x O)
  !  the "POST" in the name refers to writing the operator after the
  !  rotation/orientation in vector/tensor notation:  Q' = Q x O
  !  or in conventional quaternion notation, Qresult = Qsymm Ä QQ

  !  Thus, for active rotations (standard definition of orientation in
  !  mechanics) this is suitable for applying CRYSTAL symmetry

  if(lindex.gt.numsymm) stop 'error in presymm, lindex>numsymm'
  if(lindex.lt.1) stop 'error in presymm, lindex<1'
  qresult(1)=quatsymm(4,lindex)*qq(1) + quatsymm(1,lindex)*qq(4) + quatsymm(3,lindex)*qq(2) - quatsymm(2,lindex)*qq(3)
  qresult(2)=quatsymm(4,lindex)*qq(2) + quatsymm(2,lindex)*qq(4) + quatsymm(1,lindex)*qq(3) - quatsymm(3,lindex)*qq(1)
  qresult(3)=quatsymm(4,lindex)*qq(3) + quatsymm(3,lindex)*qq(4) + quatsymm(2,lindex)*qq(1) - quatsymm(1,lindex)*qq(2)
  qresult(4)=quatsymm(4,lindex)*qq(4) - quatsymm(1,lindex)*qq(1) - quatsymm(2,lindex)*qq(2) - quatsymm(3,lindex)*qq(3)
  !  taken from MILL2EUL

  !     write(*,*) 'postsymm: qresult ',qresult

  return
end subroutine postsymm

! _____________________________

subroutine disquat2(qq1,qq2,index1,index2,qresult,index1a,index2a)

  implicit none

  double precision qq1(4),qq2(4),qqn(4),qresult(4),quint(4),quintn(4),qmax
  double precision q4neg
  integer index1,index2,index1a,index2a,index3,iq4neg,i,j,k
  integer ijk , kjl , jjm , jkl

  !     algorithm for finding disorientation corresponding to
  !     input quaternion, QQ (from COMQUAT):  resultant quat. is in QRESULT
  !     indices for symmetry operators are in INDEX1, INDEX2

  logical keep
  integer ie , maxeul , eulcount
  parameter ( maxeul = 99000 )
  double precision euler(maxeul,3) , epseul
  double precision d1 , d2 , d3
  double precision eps
  double precision, parameter :: pi  = 4. * atan ( 1.0 ) 
  double precision rad
  double precision rodr(3) , qtemp(4)
  double precision qgrain1(4) , qgrain2(4) , axis(3) , samp1(3) , samp2(3)

  ! CODE::
  eps = 1.0e-3
  rad = pi / 180.
  index1 = 0
  index2 = 0
  index1a = 0
  index2a = 0
  index3 = 0
  qmax = -1.0

  epseul = 0.1
  eulcount = 0              ! no. of new sets of Euler angles

!!!$      print"('DISQUAT2 input 1: ',4(1x,f8.4))",qq1
!!!$      print"('DISQUAT2 input 2: ',4(1x,f8.4))",qq2
  do 200, i=1,24
     !         call presymm(qq1,i,quint)
     call postsymm(qq1,i,quint)
     !  now matches DISQUAT
     do 195, jjm=1,2
        if(jjm.eq.2) then
           do 170, k=1,4
              quint(k) = -1.0 * quint(k)
170        end do
           !     can take negative of Q because identity = +/-(0,0,0,1)
        endif

        do 190, j=1,24
           !     loop over the two sides of the grain boundary
           call postsymm(qq2,j,quintn)
           !               call presymm(qq2,j,quintn)

           call comquat2(quint,quintn,qqn)
           !               print"('DISQUAT2 qqn: ',4(1x,f8.4))",qqn
           axis(1) = qqn(1)
           axis(2) = qqn(2)
           axis(3) = qqn(3)
           call vecnorm(axis)
!!!$                  print*,' i   jjm  j ',i , jjm , j
!!!$               print"('Quat from COMQUAT2 = ',4(2x,f8.3))"
!!!$     1              ,(qqn( jkl ) , jkl = 1 , 4 )
!!!$               print"('Misor. Axis= ',3(2x,f8.3))",(axis(ijk),ijk=1,3)

           do ijk = 1 , 4
              qgrain1(ijk)= quint( ijk )
              qgrain2(ijk)= quintn( ijk )
           enddo
           call qrvec ( qgrain1 , axis , samp1 )
!!!$               print"('Qgrain1 = ',4(2x,f8.3))"
!!!$     1              ,( qgrain1( jkl ) , jkl = 1 , 4 )
!!!$               print"('misor axis in sample 1 = ',3(2x,f8.3))"
!!!$     $              , ( samp1(ijk) , ijk = 1 , 3 )
           call qrvec ( qgrain2 , axis , samp2 )
!!!$               print"('Qgrain2 = ',4(2x,f8.3))"
!!!$     1              ,( qgrain2( jkl ) , jkl = 1 , 4 )
!!!$               print"('misor axis in sample 2 = ',3(2x,f8.3))"
!!!$     $              , ( samp2(ijk) , ijk = 1 , 3 )

           do 185, kjl=1,2
              if(kjl.eq.2) then
                 do 180, k=1,4
                    qresult(k) = -1.0 * qqn(k)
180              end do
                 !     can take negative of Q because identity = +/-(0,0,0,1)
              endif

              qresult(1)=qqn(1)
              qresult(2)=qqn(2)
              qresult(3)=qqn(3)
              qresult(4)=qqn(4)
              do 210, ijk=1,2
                 if(ijk.gt.1) qresult(4) = -1.0 * qresult(4)
                 !     take the inverse rotation, i.e. apply the switching symmetry
                 !               print"('DISQUAT2 qresult: ',4(1x,f8.4))",qresult

                 call q2rod(rodr,qresult)
                 call q2eulB(d1,d2,d3,qresult)
                 d1 = d1 / rad
                 if( d1 .lt. 0. ) d1 = d1 + 360.
                 d2 = d2 / rad
                 if( d2 .lt. 0. ) d2 = d2 + 360.
                 d3 = d3 / rad
                 if( d3 .lt. 0. ) d3 = d3 + 360.
!!!$                    write(3,"(4(1x,g10.3),3i4,4(1x,g10.3))")
!!!$     $                   , d1 , d2 , d3 , 1. , i , j , kjl , qresult
                 keep = .false.
                 !                    if( i .eq. 1 .and. j .eq. 1 .and. ijk .eq. 2 )
                 !                    if( kjl .lt. 2 .and. ijk .lt. 2 )
                 if( ijk .eq. 2 )  keep = .true. !  debug
                 !                    keep = .true.
                 if ( eulcount .gt. 0 ) then
                    do ie = 1 , eulcount
                       if ( abs(d1-euler(ie,1)) .lt. epseul  &
                            .and.  abs(d2-euler(ie,2)) .lt. epseul  &
                            .and. abs(d3-euler(ie,3)).lt.epseul) then
                          keep = .false.
                       endif
                    enddo
                 endif
                 if ( keep .and. (eulcount .lt. maxeul) ) then
                    eulcount = eulcount + 1
                    euler(eulcount,1) = d1
                    euler(eulcount,2) = d2
                    euler(eulcount,3) = d3
!!!$                       print"('DISQUAT2: I/J/KJL/RodrV/Euler: '
!!!$     $                    ,3i4,6(1x,g10.3))"
!!!$     $                   ,i,j,kjl,rodr, d1 , d2 , d3
                    write(3,"(4(1x,g10.3),3i4,4(1x,g10.3))") d1 , d2 , d3 , 1. , i , j , kjl , qresult
                    write(4,"(2(1x,i4),8(1x,g10.3))")  i , j , qresult , 1. , d1 , d2 , d3
                    !      write(4,*) ' i j q1 q2 q3 q4 E1 E2 E3 / for plotting with RF_misor '
                 endif

!!!$      write(*,"('i,j,kjl,q4neg,qresult ',3i3,5f8.3)")
!!!$     $                   i,j,kjl,q4neg,qresult
                 if(qresult(1).ge.0.0) then
!!!$                        if((qresult(1)-qresult(2)).lt.eps) then
!!!$                           if((qresult(2)-qresult(3)).lt.eps) then
!!!$                              if((qresult(3)-qresult(4)).lt.eps) then
                    if((qresult(1) .le. qresult(2))) then
                       if((qresult(2) .le. qresult(3))) then
                          if((qresult(3) .le. qresult(4))) then
                             !     these IFs ensure that 0<=q1<q2<q3<q4 to place it in the fundamental zone

                             if(qresult(4).gt.qmax) then
                                !     so now we test to see if we have a smaller angle
                                qmax=qresult(4)
!!!$                                    write(*,
!!!$     $       "('FOUND LARGER q4: i,j,kjl,q4neg,qresult ',3i3,8f8.3)")
!!!$     $                                   i,j,kjl,q4neg,qresult
                                index1=i
                                index2=j

                                index1a = 0
                                if( ijk .eq. 2 ) index1a = 1
                                index2a = 0
                                if( kjl .eq. 2 ) index2a = 1
                                index3 = 0
                                if( jjm .eq. 2 ) index3 = 1
                                !     record the result if we have found a quat in the fund. zone, and a smaller angle
                             endif
                          endif
                       endif
                    endif
                 endif
210           end do
185        end do
190     end do
195  end do
200 end do

  if(index1.eq.0.or.index2.eq.0) stop 'DISQUAT2 failed!'

!!!$      print"('DISQUAT2 input 1: ',4(1x,f8.4))",qq1
!!!$      print"('DISQUAT2 input 2: ',4(1x,f8.4))",qq2
  !      call presymm(qq1,index1,quint)
  call postsymm(qq1,index1,quint)
  if ( index3 .gt. 0 ) then
     quint(1) = -1.0 * quint(1)
     quint(2) = -1.0 * quint(2)
     quint(3) = -1.0 * quint(3)
     quint(4) = -1.0 * quint(4)
  endif

  call postsymm(qq2,index2,quintn)

  call comquat2( quint , quintn , qqn )
!!!$      print"('DISQUAT2 interm. result: ',4(1x,f8.4))",qqn

  do k=1,4
     if(index2a.eq.1) then
        !  if we had to take the negative of the Q, then do so
        !  to the answer

        qresult(k) = -1.0 * qqn(k)
        !  can take negative of Q because identity = +/-(0,0,0,1)
     else
        qresult(k)=qqn(k)
     endif
  enddo

  if(index1a.eq.1) then
     qresult(4) = -1. * qresult(4)
     !  restore the 4th component if disorient found for ijk=1
  endif

!!!$      print"('DISQUAT2 result: ',4(1x,f8.4))",qresult
!!!$      print"('DISQUAT2 result: ',4i4)",index1 , index2 ,index1a,index2a

  return
end subroutine disquat2


!     _____________________________

subroutine disquat(qq,index1,index2,qresult,index1a,index2a)

  implicit none

  double precision qq(4),qqn(4),qresult(4),quint(4),quintn(4),qmax
  integer index1,index2,index1a,index2a,index3,iq4neg,i,j,k
  integer kjl , ijk
  double precision switch , rodr(3)

  !  algorithm for finding disorientation corresponding to
  !  input quaternion, QQ (from COMQUAT):  resultant quat. is in QRESULT
  !  indices for symmetry operators are in INDEX1, INDEX2

  logical keep
  integer ie , maxeul , eulcount
  parameter ( maxeul = 12000 )
  double precision euler(maxeul,3) , epseul
  double precision d1 , d2 , d3
  double precision eps , rad
  double precision, parameter :: pi  = 4.d0 * atan ( 1.0d0 ) 

  ! CODE::
  eps = 1.0e-3
  rad = pi / 180.
  index1 = 0
  index2 = 0
  index1a = 0
  index2a = 0
  qmax = -1.

  epseul = 0.1
  eulcount = 0 ! no. of new sets of Euler angles

  do 200, i = 1 , 24
     call presymm(qq,i,quint)

     do 190, j=1,24
        !  loop over the two sides of the grain
        call postsymm(quint,j,quintn)

        do 185, kjl=1,2
           switch=1.0
           if(kjl.eq.2) then
              switch=-1.0
           endif
           do 180, k=1,4
              qqn(k)=quintn(k)*switch
180        end do
           !  can take negative of Q because identity = +/-(0,0,0,1)

           do 210, ijk=1,2
              qresult(1)=qqn(1)
              qresult(2)=qqn(2)
              qresult(3)=qqn(3)
              qresult(4)=qqn(4)
              iq4neg=ijk
              if(ijk.eq.2) qresult(4)=-1.*qqn(4)
              !  take the inverse rotation, i.e. apply the switching symmetry

              call q2rod(rodr,qresult)
              call q2eulB(d1,d2,d3,qresult)
              d1 = d1 / rad
              if( d1 .lt. 0. ) d1 = d1 + 360.
              d2 = d2 / rad
              if( d2 .lt. 0. ) d2 = d2 + 360.
              d3 = d3 / rad
              if( d3 .lt. 0. ) d3 = d3 + 360.
!!!$  write(3,"(4(1x,g10.3),3i4,4(1x,g10.3))")
!!!$  $                   , d1 , d2 , d3 , 1. , i , j , kjl , qresult
              keep = .false.
              !     if( i .eq. 1 .and. j .eq. 1 .and. ijk .eq. 2 )
              !     if( kjl .lt. 2 .and. ijk .lt. 2 )
              if( ijk .eq. 2 )  keep = .true. !  debug
              if ( eulcount .gt. 0 ) then
                 do ie = 1 , eulcount
                    if ( abs(d1-euler(ie,1)) .lt. epseul  &
                         .and.  abs(d2-euler(ie,2)) .lt. epseul  &
                         .and. abs(d3-euler(ie,3)).lt.epseul) then
                       keep = .false.
                    endif
                 enddo
              endif
              if ( keep .and. (eulcount .le. maxeul) ) then
                 eulcount = eulcount + 1
                 euler(eulcount,1) = d1
                 euler(eulcount,2) = d2
                 euler(eulcount,3) = d3
!!!$                       print"('DISQ_HEX: I/J/KJL/RodrV/Euler: '
!!!$     $                    ,3i4,6(1x,g10.3))"
!!!$     $                   ,i,j,kjl,rodr, d1 , d2 , d3
                 write(3,"(4(1x,g10.3),3i4,4(1x,g10.3))") d1 , d2 , d3 , 1. , i , j , kjl , qresult
                 write(4,"(2(1x,i4),8(1x,g10.3))") i , j , qresult , 1. , d1 , d2 , d3
                 !      write(4,*) ' i j q1 q2 q3 q4 E1 E2 E3 / for plotting with RFdist '
              endif

              if( qresult(1) .ge. 0.0 ) then
                 !			  if((qresult(1)-qresult(2)).lt.1e-5) then
                 !			    if((qresult(2)-qresult(3)).lt.1e-5) then
                 !			      if((qresult(3)-qresult(4)).lt.1e-5) then
                 if(qresult(1).le.qresult(2)) then
                    if(qresult(2).le.qresult(3)) then
                       if(qresult(3).le.qresult(4)) then
                          !   these IFs ensure that 0<=q1<q2<q3<q4 to place it in the fundamental zone
                          if(qresult(4).gt.qmax) then
                             !  so now we test to see if we have a smaller angle
                             qmax=qresult(4)
                             index1=i
                             index2=j
                             !			        index3=0
                             index1a=0
                             if(ijk.eq.2) index1a=1
                             index2a=0
                             if(kjl.eq.2) index2a=1
                             !  record the result if we have found a quat in the fund. zone, and a smaller angle
                          endif
                       endif
                    endif
                 endif
              endif
210        end do
185     end do
190  end do
200 end do

  if(index1.eq.0.or.index2.eq.0) stop 'DISQUAT failed!'

  call presymm(qq,index1,quint)

  call postsymm(quint,index2,qqn)
  !		 write(*,"('DISQUAT: qqn = ',4f8.3)") qqn
  do k=1,4
     if(index2a.eq.1) then
        !  if we had to take the negative of the Q, then do so
        !  to the answer

        qresult(k) = -1.0 * qqn(k)
        !  can take negative of Q because identity = +/-(0,0,0,1)
     else
        qresult(k)=qqn(k)
     endif
  enddo

  if(index1a.eq.1) then
     qresult(4) = -1. * qresult(4)
     !  restore the 4th component if disorient found for ijk=1
  endif

  !	write(*,*) 'DISQUAT: qmax= ',qmax
  return
end subroutine disquat

! _____________________________

subroutine disquat_hex(qq,index1,index2,qresult,index1a,index2a)

  !  *****************************************   HCP !!
  !  this version for HCP,  ADR  25 vi 05
  !  the shape of the fundamental zone (FZ) for hcp is
  !  a prismatic box, with opt and bottom faces at +/- pi/2*6 = 30 degrees
  !  such that the limits on R3 are independent of R1 and R2
  !  and square prism faces at a distance 1 (+/- pi/2*4 = 45 degrees) from the center

  !  algorithm for finding disorientation corresponding to
  !  input quaternion, QQ (from COMQUAT):  resultant quat. is in QRESULT
  !  indices for symmetry operators are in INDEX1, INDEX2

  implicit none

  double precision qq(4),qqn(4),qresult(4),quint(4),quintn(4),qmax
  integer index1,index2,index1a,index2a,index3,iq4neg,i,j,k
  integer ijk , kjl , ik
  double precision switch
  double precision :: rtmp,eps
  double precision tan60,sin60,cos60,tan30,sin30,cos30
  double precision rodr(3)
  logical tests(20)
  double precision rad , qtmp(4) , d1 , d2 , d3 , rnorm

  logical keep
  integer ie , maxeul , eulcount
  parameter ( maxeul = 12000 )
  double precision euler(maxeul,3) , epseul
  double precision, parameter :: pi  = 4.d0 * atan ( 1.0d0 ) 

  !  CODE::

  eps = 1.0e-5
  epseul = 0.1
  rad = pi / 180.
  tan60 = tan ( rad * 60. )
  sin60 = sin ( rad * 60. )
  cos60 = cos ( rad * 60. )
  tan30 = tan ( rad * 30. )
  sin30 = sin ( rad * 30. )
  cos30 = cos ( rad * 30. )

  index1=0
  index2=0
  index1a=0
  index2a=0
  qmax=0.

  eulcount = 0 ! no. of new sets of Euler angles

  !	do 200, i=1,24
  do 200, i = 1 , 12  !   HCP
     call presymm(qq,i,quint)
     rnorm = sqrt(quint(1)**2 + quint(2)**2 + quint(3)**2 + quint(4)**2)
     if( abs( rnorm - 1.0 ) .gt. 1.e-3 ) then
        write(*,*) 'quint required normalization'
        do ik = 1 , 4
           quint(ik) = quint(ik) / rnorm
        enddo
     endif

     !       do 190, j=1,24
     do 190, j = 1 , 12 !   HCP
        !  loop over the two sides of the grain (boundary)
        call postsymm(quint,j,quintn)
        rnorm = sqrt(quintn(1)**2 + quintn(2)**2 + quintn(3)**2 + quintn(4)**2)
        if( abs( rnorm - 1.0 ) .gt. 1.e-3 ) then
           write(*,*) 'quint required normalization'
           do ik = 1 , 4
              quintn(ik) = quintn(ik) / rnorm
           enddo
        endif


        do 185, kjl = 1 , 2
           switch = 1.0
           if( kjl .eq. 2 ) then
              switch = -1.0
           endif
           do 180, k = 1 , 4
              qqn(k) = quintn(k) * switch
180        end do
           !  can take negative of Q because identity = +/-(0,0,0,1)
           rnorm = sqrt(qqn(1)**2 + qqn(2)**2 + qqn(3)**2 + qqn(4)**2)
           if( abs( rnorm - 1.0 ) .gt. 1.e-3 ) then
              write(*,*) 'qqn required normalization'
              do ik = 1 , 4
                 qqn(ik) = qqn(ik) / rnorm
              enddo
           endif


           do  ijk = 1 , 2
              qresult(1) = qqn(1)
              qresult(2) = qqn(2)
              qresult(3) = qqn(3)
              qresult(4) = qqn(4)
              iq4neg = ijk
              if( ijk .eq. 2 ) qresult(4) = -1. * qqn(4)
              !  take the inverse rotation, i.e. apply the switching symmetry
              rnorm = sqrt(qresult(1)**2 + qresult(2)**2 + qresult(3)**2 + qresult(4)**2)
              if( abs( rnorm - 1.0 ) .gt. 1.e-3 ) then
                 write(*,*) 'qresult required normalization'
                 do ik = 1 , 4
                    qresult(ik) = qresult(ik) / rnorm
                 enddo
              endif

              call q2rod(rodr,qresult)
              call q2eulB(d1,d2,d3,qresult)
              d1 = d1 / rad
              if( d1 .lt. 0. ) d1 = d1 + 360.
              d2 = d2 / rad
              if( d2 .lt. 0. ) d2 = d2 + 360.
              d3 = d3 / rad
              if( d3 .lt. 0. ) d3 = d3 + 360.
!!!$                    write(3,"(4(1x,g10.3),3i4,4(1x,g10.3))")
!!!$     $                   , d1 , d2 , d3 , 1. , i , j , kjl , qresult
              keep = .false.
              !                    if( i .eq. 1 .and. j .eq. 1 .and. ijk .eq. 2 )
              !                    if( kjl .lt. 2 .and. ijk .lt. 2 )
              if( ijk .eq. 2 )  keep = .true. !  debug
              !                    keep = .true.
              if ( eulcount .gt. 0 ) then
                 do ie = 1 , eulcount
                    if ( abs(d1-euler(ie,1)) .lt. epseul  &
                         .and.  abs(d2-euler(ie,2)) .lt. epseul  &
                         .and. abs(d3-euler(ie,3)).lt.epseul) then
                       keep = .false.
                    endif
                 enddo
              endif
              if ( keep .and. (eulcount .le. maxeul) ) then
                 eulcount = eulcount + 1
                 euler(eulcount,1) = d1
                 euler(eulcount,2) = d2
                 euler(eulcount,3) = d3
!!!$                       print"('DISQ_HEX: I/J/KJL/RodrV/Euler: '
!!!$     $                    ,3i4,6(1x,g10.3))"
!!!$     $                   ,i,j,kjl,rodr, d1 , d2 , d3
                 write(3,"(4(1x,g10.3),3i4,4(1x,g10.3))") d1 , d2 , d3 , 1. , i , j , kjl , qresult
                 write(4,"(2(1x,i4),8(1x,g10.3))")  i , j , qresult , 1. , d1 , d2 , d3
                 !      write(4,*) ' i j q1 q2 q3 q4 E1 E2 E3 / for plotting with RFdist '
              endif

              do k = 1 , 20
                 tests(k) = .false.
              enddo

              if(rodr(2).ge. ( -1. * eps ) ) then
                 tests(1) = .true.
                 if(rodr(1).le. ( 1.0 + eps ) ) then
                    tests(2) = .true.

                    !			  if((qresult(1)-qresult(2)).lt.1e-5) then
                    !			    if((qresult(2)-qresult(3)).lt.1e-5) then
                    !			      if((qresult(3)-qresult(4)).lt.1e-5) then
                    !			  if(qresult(2).le.qresult(1)) then

                    if(rodr(3).ge. (-1.*eps) ) then
                       tests(3) = .true.
                       if(rodr(3).le. (tan30+eps) ) then
                          !  these two IFs put the point inside the FZ for R3
                          tests(4) = .true.

                          rtmp = rodr(1)*tan30 - rodr(2)
                          if(rtmp.ge. (-1.*eps) ) then
                             tests(5) = .true.
                             !  line from origin at 30 degr
                             rtmp = rodr(1)*cos30 + rodr(2)*sin30
                             if(rtmp.le. (1.+eps) ) then
                                tests(6) = .true.
                                !  1st sector
                                rtmp = rodr(1) * cos60 + rodr(2)*sin60
                                if(rtmp.le. (1.+eps) ) then
                                   tests(7) = .true.
                                   !  2nd sector
                                   if((qresult(4)-qmax).gt.eps) then
                                      tests(8) = .true.
                                      qmax=qresult(4)
                                      !  so now we test to see if we have a smaller angle
                                      index1=i
                                      index2=j
                                      !       index3=0
                                      index1a=0
                                      if(ijk.eq.2) index1a=1
                                      index2a=0
                                      if(kjl.eq.2) index2a=1
                                      !  record the result if we have found a quat in the fund. zone, and a smaller angle
                                   endif
                                endif
                             endif
                          endif
                       endif
                    endif
                 endif
              endif
!!!$		    print*,'TESTS: ',(tests(k),k=1,8)
!!!$		    if(tests(8)) print*,'ALL = TRUE'
210        enddo
185     end do  !  see commenting above
190  end do
200 end do

  if(index1.eq.0.or.index2.eq.0) then
     write(*,"('DISQUAT: qq (input) = ',4f8.3)") qq
     stop 'DISQUAT_HEX failed!'
  endif

  call presymm(qq,index1,quint)

  call postsymm(quint,index2,qqn)
  !		 write(*,"('DISQUAT: qqn = ',4f8.3)") qqn
  do 280, k=1,4
     if(index2a.eq.1) then
        !  if we had to take the negative of the Q, then do so
        !  to the answer

        qresult(k) = -1.0 * qqn(k)
        !  can take negative of Q because identity = +/-(0,0,0,1)
     else
        qresult(k)=qqn(k)
     endif
280 end do

  if(index1a.eq.1) then
     qresult(4) = -1. * qresult(4)
     !  restore the 4th component if disorient found for ijk=1
  endif

  write(*,*) 'DISQUAT-HEX: qmax = ' , qmax
  write(*,*) 'DISQUAT-HEX: eulcount = ', eulcount
  return
end subroutine disquat_hex

! _____________________________

subroutine disquat_ang_hex(qq,angle)

  !  *****************************************   HCP !!
  !  this version for HCP,  ADR  25 vi 05
  !  returning only an angle [degrees], not a quaternion

  !  the shape of the fundamental zone (FZ) for hcp is
  !  a prismatic box, with opt and bottom faces at +/- pi/2*6 = 30 degrees
  !  such that the limits on R3 are independent of R1 and R2
  !  and square prism faces at a distance 1 (+/- pi/2*4 = 45 degrees) from the center

  !  algorithm for finding disorientation corresponding to
  !  input quaternion, QQ (from COMQUAT):  resultant quat. is in QRESULT
  !  indices for symmetry operators are in INDEX1, INDEX2

  implicit none

  double precision qq(4),qqn(4),qresult(4),quint(4),quintn(4),qmax
  integer index1,index2,index1a,index2a,index3,iq4neg,i,j,k
  integer ijk,kjl
  double precision switch
  double precision :: rtmp,eps
  double precision tan60,sin60,cos60,tan30,sin30,cos30
  double precision rodr(3)
  logical tests(20)
  double precision angle
  double precision, parameter :: pi  = 4.d0 * atan ( 1.0d0 ) 

  !  CODE::

  eps = -1.0e-5
  tan60 = tan (pi*60./180.)
  sin60 = sin (pi*60./180.)
  cos60 = cos (pi*60./180.)
  tan30 = tan (pi*30./180.)
  sin30 = sin (pi*30./180.)
  cos30 = cos (pi*30./180.)

  index1=0
  index2=0
  index1a=0
  index2a=0
  qmax = 0.

  !	do 200, i=1,24
  do 200, i=1,12  !   HCP
     call presymm(qq,i,quint)

     !       do 190, j=1,24
!!!$	   do 190, j=1,12	!   HCP
!!!$c  loop over the two sides of the grain
!!!$	      call postsymm(quint,j,quintn)

     do 185, kjl=1,2
        switch=1.0
        if(kjl.eq.2) then
           switch=-1.0
        endif
        do 180, k=1,4
!!!$		    qqn(k)=quintn(k)*switch
           qqn(k)=quint(k)*switch
180     end do
        !  can take negative of Q because identity = +/-(0,0,0,1)

        do  ijk=1,2
           qresult(1)=qqn(1)
           qresult(2)=qqn(2)
           qresult(3)=qqn(3)
           qresult(4)=qqn(4)
           iq4neg=ijk
           if(ijk.eq.2) qresult(4)=-1.*qqn(4)
           !  take the inverse rotation, i.e. apply the switching symmetry

!!!$		    call q2rod(rodr,qresult)
!!!$		    print"('I/J/KJL/Rodr V: ',3i4,3(1x,g10.3))",i,j,kjl,rodr
!!!$
!!!$		    do k = 1,20
!!!$		       tests(k) = .false.
!!!$		    enddo

!!!$		    if(rodr(2).ge.0.0) then
!!!$		       tests(1) = .true.
!!!$		       if(rodr(1).le.1.0) then
!!!$			  tests(2) = .true.

           !			  if((qresult(1)-qresult(2)).lt.1e-5) then
           !			    if((qresult(2)-qresult(3)).lt.1e-5) then
           !			      if((qresult(3)-qresult(4)).lt.1e-5) then
           !			  if(qresult(2).le.qresult(1)) then

!!!$			     if(rodr(3).ge.0.) then
!!!$				tests(3) = .true.
!!!$				if(rodr(3).le.tan30) then
!!!$!  these two IFs put the point inside the FZ for R3
!!!$				   tests(4) = .true.
!!!$
!!!$				   rtmp = rodr(1)*tan30 - rodr(2)
!!!$				   if(rtmp.ge.0.) then
!!!$				      tests(5) = .true.
!!!$c  line from origin at 30 degr
!!!$				      rtmp = rodr(1)*cos30 + rodr(2)*sin30
!!!$				      if(rtmp.le.1.) then
!!!$					 tests(6) = .true.
!!!$c  1st sector
!!!$					 rtmp = rodr(1)*cos60 + rodr(2)*sin60
!!!$					 if(rtmp.le.1.) then
!!!$					    tests(7) = .true.
           !  2nd sector
           if((qresult(4)-qmax).gt.eps) then
!!!$			     tests(8) = .true.
              qmax=qresult(4)
              !  so now we test to see if we have a smaller angle
!!!$					    index1=i
!!!$					    index2=j
              !       index3=0
!!!$					    index1a=0
!!!$					    if(ijk.eq.2) index1a=1
!!!$					    index2a=0
!!!$					    if(kjl.eq.2) index2a=1
              !  record the result if we have found a quat in the fund. zone, and a smaller angle
           endif
!!!$				      endif
!!!$				   endif
!!!$				endif
!!!$			     endif
!!!$			  endif
!!!$		       endif
!!!$		    endif
!!!$		    print*,'TESTS: ',(tests(k),k=1,8)
!!!$		    if(tests(8)) print*,'ALL = TRUE'
        enddo
185  end do
!!!$ 190	   continue
200 end do

!!!$	if(index1.eq.0.or.index2.eq.0) stop 'DISQUAT failed!'
!!!$
!!!$	call presymm(qq,index1,quint)
!!!$
!!!$	call postsymm(quint,index2,qqn)
!!!$c		 write(*,"('DISQUAT: qqn = ',4f8.3)") qqn
!!!$	do 280, k=1,4
!!!$	   if(index2a.eq.1) then
!!!$c  if we had to take the negative of the Q, then do so
!!!$c  to the answer
!!!$
!!!$	      qresult(k)=qqn(k)*-1.0
!!!$c  can take negative of Q because identity = +/-(0,0,0,1)
!!!$	   else
!!!$	      qresult(k)=qqn(k)
!!!$	   endif
!!!$ 280	continue
!!!$
!!!$	if(index1a.eq.1) then
!!!$	   qresult(4)=-1.*qresult(4)
!!!$c  restore the 4th component if disorient found for ijk=1
!!!$	endif

  !	write(*,*) 'DISQUAT_ANG: qmax= ',qmax
  if(qmax.gt.1.) qmax = 1.
  angle = 360. / pi * acos(qmax)

  return
end subroutine disquat_ang_hex


!     _____________________________


subroutine disquat_orth(qq,index1,index2,qresult ,index1a,index2a)

  implicit none

  double precision :: qq(4),qqn(4),qresult(4),quint(4),quintn(4),qmax
  integer :: index1,index2,index1a,index2a,index3,iq4neg
  integer :: i,j,k,ip,jp , ijk , kjl
  double precision :: switch , rodr(3) , d1 , d2 , d3 , rad
  !  algorithm for finding disorientation corresponding to
  !  input quaternion, QQ (from COMQUAT):  resultant quat. is in QRESULT
  !  indices for symmetry operators are in INDEX1, INDEX2

  !  edited for orthorhombic crystal symmetry, Jul 07

  logical :: keep
  integer :: ie , maxeul , eulcount
  parameter ( maxeul = 6000 )
  double precision :: euler(maxeul,3)
  double precision, parameter :: PI = 4.d0 * atan( 1.d0 )
  double precision :: epseul = 0.1

  !  CODE::
  rad = pi / 180.
  index1 = 0
  index2 = 0
  index1a = 0
  index2a = 0
  qmax = -1.  !  allow for 180 since max angle can be > 90
  ip = 1
  jp = 1

  eulcount = 0 ! no. of new sets of Euler angles

  do 200, i = 1 , 4  !  ortho
     call presymm(qq,i,quint)

     do 190, j = 1 , 4        !  ortho
        !  loop over the two sides of the grain
        call postsymm(quint,j,quintn)

        do 185, kjl=1,2
           switch=1.0
           if(kjl.eq.2) then
              switch=-1.0
           endif
           do k = 1 , 4
              qqn(k)=quintn(k)*switch
           end do
           !  can take negative of Q because identity = +/-(0,0,0,1)

           do 210, ijk=1,2
              qresult(1)=qqn(1)
              qresult(2)=qqn(2)
              qresult(3)=qqn(3)
              qresult(4)=qqn(4)
              iq4neg=ijk
              if(ijk.eq.2) qresult(4) = -1. * qqn(4)
              !  take the inverse rotation, i.e. apply the switching symmetry

              call q2rod(rodr,qresult)
              call q2eulB(d1,d2,d3,qresult)
              d1 = d1 / rad
              if( d1 .lt. 0. ) d1 = d1 + 360.
              d2 = d2 / rad
              if( d2 .lt. 0. ) d2 = d2 + 360.
              d3 = d3 / rad
              if( d3 .lt. 0. ) d3 = d3 + 360.
!!!$                    write(*,"('result: ',4(1x,g10.3),4i4,4(1x,g10.4))")
!!!$     $                 , d1 , d2 , d3 , 1. , i , j
!!!$     $                 , kjl , ijk , qresult
              keep = .false.
              !                    if( i .eq. 1 .and. j .eq. 1 .and. ijk .eq. 2 )
              !                    if( kjl .lt. 2 .and. ijk .lt. 2 )
              if( ijk .eq. 2 )  keep = .true. !  debug
              !                    keep = .true.
              if ( eulcount .gt. 0 ) then
                 do ie = 1 , eulcount
                    if ( abs(d1-euler(ie,1)) .lt. epseul .and.  abs(d2-euler(ie,2)) .lt. epseul  &
                         .and. abs(d3-euler(ie,3)).lt.epseul) then
                       keep = .false.
                    endif
                 enddo
              endif
              if ( keep .and. (eulcount .le. maxeul) ) then
                 eulcount = eulcount + 1
                 euler(eulcount,1) = d1
                 euler(eulcount,2) = d2
                 euler(eulcount,3) = d3
!!!$                       print"('DISQ_HEX: I/J/KJL/RodrV/Euler: '
!!!$     $                    ,3i4,6(1x,g10.3))"
!!!$     $                   ,i,j,kjl,rodr, d1 , d2 , d3
                 write(3,"(4(1x,g10.3),3i4,4(1x,g10.3))") d1 , d2 , d3 , 1. , i , j , kjl , qresult
                 write(4,"(2(1x,i4),8(1x,g10.3))") i , j , qresult , 1. , d1 , d2 , d3
                 !      write(4,*) ' i j q1 q2 q3 q4 E1 E2 E3 / for plotting with RFdist '
              endif

              if(qresult(1).ge.0.0) then
                 if(qresult(2).ge.0.0) then
                    if(qresult(3).ge.0.0) then
                       !   these IFs ensure that 0<=q1, 0<=q2, 0<=q3 to place it in the fundamental zone
                       !  which is a cube of size=1., between 0 and 1 in all 3 axes
                       if(qresult(4).gt.qmax) then
                          !  so now we test to see if we have a smaller angle, which closes the box
                          qmax=qresult(4)
                          index1=i
                          index2=j
                          !               index3=0
                          index1a=0
                          if(ijk.eq.2) index1a=1
                          index2a=0
                          if(kjl.eq.2) index2a=1
                          !  record the result if we have found a quat in the fund. zone, and a smaller angle
                       endif
                    endif
                 endif
              endif
210        end do
185     end do
190  end do
200 end do

  if(index1.eq.0.or.index2.eq.0) stop 'DISQUAT_ORTH failed!'

  call presymm(qq,index1,quint)

  call postsymm(quint,index2,qqn)
  !       write(*,"('DISQUAT: qqn = ',4f8.3)") qqn
  do 280, k=1,4
     if(index2a.eq.1) then
        !  if we had to take the negative of the Q, then do so
        !  to the answer

        qresult(k) = -1.0 * qqn(k)
        !  can take negative of Q because identity = +/-(0,0,0,1)
     else
        qresult(k)=qqn(k)
     endif
280 end do

  if(index1a.eq.1) then
     qresult(4) = -1. * qresult(4)
     !  restore the 4th component if disorient found for ijk=1
  endif

  write(*,*) 'DISQUAT_ORTH: qmax= ',qmax
  write(*,"('DISQUAT_ORTH: qresult = ',4f8.3)") qresult
  return
end subroutine disquat_orth

!     _____________________________

!     _____________________________


subroutine disquat_tet( qq , index1 , index2 , qresult , index1a , index2a )

  implicit none

  double precision :: qq(4),qqn(4),qresult(4),quint(4),quintn(4),qmax
  integer :: index1,index2,index1a,index2a,index3,iq4neg
  integer :: i,j,k,ip,jp , ijk , kjl
  double precision :: switch , rodr(3) , d1 , d2 , d3 , rad
  !  algorithm for finding disorientation corresponding to
  !  input quaternion, QQ (from COMQUAT):  resultant quat. is in QRESULT
  !  indices for symmetry operators are in INDEX1, INDEX2

  !  edited for orthorhombic crystal symmetry, Jul 07

  logical :: keep
  integer :: ie , maxeul , eulcount
  parameter ( maxeul = 6000 )
  double precision :: euler(maxeul,3)
  double precision, parameter :: PI = 4.d0 * atan( 1.d0 )
  double precision :: epseul = 0.1

  !  CODE::
  rad = pi / 180.
  index1 = 0
  index2 = 0
  index1a = 0
  index2a = 0
  qmax = -1.  !  allow for 180 since max angle can be > 90
  ip = 1
  jp = 1

  eulcount = 0 ! no. of new sets of Euler angles

  do 200, i = 1 , 8  !  tet 4/mmm
     call presymm(qq,i,quint)

     do 190, j = 1 , 8        !  tet  4/mmm
        !  loop over the two sides of the grain
        call postsymm(quint,j,quintn)

        do 185, kjl=1,2
           switch=1.0
           if(kjl.eq.2) then
              switch=-1.0
           endif
           do k = 1 , 4
              qqn(k)=quintn(k)*switch
           end do
           !  can take negative of Q because identity = +/-(0,0,0,1)

           do 210, ijk=1,2
              qresult(1)=qqn(1)
              qresult(2)=qqn(2)
              qresult(3)=qqn(3)
              qresult(4)=qqn(4)
              iq4neg=ijk
              if(ijk.eq.2) qresult(4) = -1. * qqn(4)
              !  take the inverse rotation, i.e. apply the switching symmetry

              call q2rod(rodr,qresult)
              call q2eulB(d1,d2,d3,qresult)
              d1 = d1 / rad
              if( d1 .lt. 0. ) d1 = d1 + 360.
              d2 = d2 / rad
              if( d2 .lt. 0. ) d2 = d2 + 360.
              d3 = d3 / rad
              if( d3 .lt. 0. ) d3 = d3 + 360.
!!!$                    write(*,"('result: ',4(1x,g10.3),4i4,4(1x,g10.4))")
!!!$     $                 , d1 , d2 , d3 , 1. , i , j
!!!$     $                 , kjl , ijk , qresult
              keep = .false.
              !                    if( i .eq. 1 .and. j .eq. 1 .and. ijk .eq. 2 )
              !                    if( kjl .lt. 2 .and. ijk .lt. 2 )
              if( ijk .eq. 2 )  keep = .true. !  debug
              !                    keep = .true.
              if ( eulcount .gt. 0 ) then
                 do ie = 1 , eulcount
                    if ( abs(d1-euler(ie,1)) .lt. epseul .and.  abs(d2-euler(ie,2)) .lt. epseul  &
                         .and. abs(d3-euler(ie,3)).lt.epseul) then
                       keep = .false.
                    endif
                 enddo
              endif
              if ( keep .and. (eulcount .le. maxeul) ) then
                 eulcount = eulcount + 1
                 euler(eulcount,1) = d1
                 euler(eulcount,2) = d2
                 euler(eulcount,3) = d3
!!!$                       print"('DISQ_HEX: I/J/KJL/RodrV/Euler: '
!!!$     $                    ,3i4,6(1x,g10.3))"
!!!$     $                   ,i,j,kjl,rodr, d1 , d2 , d3
                 write(3,"(4(1x,g10.3),3i4,4(1x,g10.3))") d1 , d2 , d3 , 1. , i , j , kjl , qresult
                 write(4,"(2(1x,i4),8(1x,g10.3))") i , j , qresult , 1. , d1 , d2 , d3
                 !      write(4,*) ' i j q1 q2 q3 q4 E1 E2 E3 / for plotting with RFdist '
              endif

              if(qresult(1).ge.0.0) then
                 if(qresult(2).ge.0.0) then
                    if(qresult(3).ge.0.0) then
                       !   these IFs ensure that 0<=q1, 0<=q2, 0<=q3 to place it in the fundamental zone
                       !  which is a cube of size=1., between 0 and 1 in all 3 axes
                       if(qresult(4).gt.qmax) then
                          !  so now we test to see if we have a smaller angle, which closes the box
                          qmax=qresult(4)
                          index1=i
                          index2=j
                          !               index3=0
                          index1a=0
                          if(ijk.eq.2) index1a=1
                          index2a=0
                          if(kjl.eq.2) index2a=1
                          !  record the result if we have found a quat in the fund. zone, and a smaller angle
                       endif
                    endif
                 endif
              endif
210        end do
185     end do
190  end do
200 end do

  if(index1.eq.0.or.index2.eq.0) stop 'DISQUAT_TET failed!'

  call presymm(qq,index1,quint)

  call postsymm(quint,index2,qqn)
  !       write(*,"('DISQUAT: qqn = ',4f8.3)") qqn
  do 280, k=1,4
     if(index2a.eq.1) then
        !  if we had to take the negative of the Q, then do so
        !  to the answer

        qresult(k) = -1.0 * qqn(k)
        !  can take negative of Q because identity = +/-(0,0,0,1)
     else
        qresult(k)=qqn(k)
     endif
280 end do

  if(index1a.eq.1) then
     qresult(4) = -1. * qresult(4)
     !  restore the 4th component if disorient found for ijk=1
  endif

  write(*,*) 'DISQUAT_TET: qmax= ',qmax
  write(*,"('DISQUAT_TET: qresult = ',4f8.3)") qresult
  return
end subroutine disquat_tet

!     _____________________________


subroutine disquat_tet4onM( qq , index1 , index2 , qresult , index1a , index2a )

  implicit none

  double precision :: qq(4),qqn(4),qresult(4),quint(4),quintn(4),qmax
  integer :: index1,index2,index1a,index2a,index3,iq4neg
  integer :: i,j,k,ip,jp , ijk , kjl
  double precision :: switch , rodr(3) , d1 , d2 , d3 , rad
  !  algorithm for finding disorientation corresponding to
  !  input quaternion, QQ (from COMQUAT):  resultant quat. is in QRESULT
  !  indices for symmetry operators are in INDEX1, INDEX2

  !  edited for tetragonal 4/m crystal symmetry, Mar 15
  !  only the 4 90 degr operators

  logical :: keep
  integer :: ie , maxeul , eulcount
  parameter ( maxeul = 6000 )
  double precision :: euler(maxeul,3)
  double precision, parameter :: PI = 4.d0 * atan( 1.d0 )
  double precision :: epseul = 0.1

  !  CODE::
  rad = pi / 180.
  index1 = 0
  index2 = 0
  index1a = 0
  index2a = 0
  qmax = -1.  !  allow for 180 since max angle can be > 90
  ip = 1
  jp = 1

  eulcount = 0 ! no. of new sets of Euler angles

  do 200, i = 1 , 4  !  tet
     call presymm(qq,i,quint)

     do 190, j = 1 , 4        !  tet
        !  loop over the two sides of the grain
        call postsymm(quint,j,quintn)

        do 185, kjl=1,2
           switch=1.0
           if(kjl.eq.2) then
              switch=-1.0
           endif
           do k = 1 , 4
              qqn(k)=quintn(k)*switch
           end do
           !  can take negative of Q because identity = +/-(0,0,0,1)

           do 210, ijk=1,2
              qresult(1)=qqn(1)
              qresult(2)=qqn(2)
              qresult(3)=qqn(3)
              qresult(4)=qqn(4)
              iq4neg=ijk
              if(ijk.eq.2) qresult(4) = -1. * qqn(4)
              !  take the inverse rotation, i.e. apply the switching symmetry

              call q2rod(rodr,qresult)
              call q2eulB(d1,d2,d3,qresult)
              d1 = d1 / rad
              if( d1 .lt. 0. ) d1 = d1 + 360.
              d2 = d2 / rad
              if( d2 .lt. 0. ) d2 = d2 + 360.
              d3 = d3 / rad
              if( d3 .lt. 0. ) d3 = d3 + 360.
!!!$                    write(*,"('result: ',4(1x,g10.3),4i4,4(1x,g10.4))")
!!!$     $                 , d1 , d2 , d3 , 1. , i , j
!!!$     $                 , kjl , ijk , qresult
              keep = .false.
              !                    if( i .eq. 1 .and. j .eq. 1 .and. ijk .eq. 2 )
              !                    if( kjl .lt. 2 .and. ijk .lt. 2 )
              if( ijk .eq. 2 )  keep = .true. !  debug
              !                    keep = .true.
              if ( eulcount .gt. 0 ) then
                 do ie = 1 , eulcount
                    if ( abs(d1-euler(ie,1)) .lt. epseul .and.  abs(d2-euler(ie,2)) .lt. epseul  &
                         .and. abs(d3-euler(ie,3)).lt.epseul) then
                       keep = .false.
                    endif
                 enddo
              endif
              if ( keep .and. (eulcount .le. maxeul) ) then
                 eulcount = eulcount + 1
                 euler(eulcount,1) = d1
                 euler(eulcount,2) = d2
                 euler(eulcount,3) = d3
!!!$                       print"('DISQ_HEX: I/J/KJL/RodrV/Euler: '
!!!$     $                    ,3i4,6(1x,g10.3))"
!!!$     $                   ,i,j,kjl,rodr, d1 , d2 , d3
                 write(3,"(4(1x,g10.3),3i4,4(1x,g10.3))") d1 , d2 , d3 , 1. , i , j , kjl , qresult
                 write(4,"(2(1x,i4),8(1x,g10.3))") i , j , qresult , 1. , d1 , d2 , d3
                 !      write(4,*) ' i j q1 q2 q3 q4 E1 E2 E3 / for plotting with RFdist '
              endif

              if(qresult(1).ge.0.0) then
                 if(qresult(2).ge.0.0) then
                    !if(qresult(3).ge.0.0) then
                       !   these IFs ensure that 0<=q1, 0<=q2 to place it in the fundamental zone
                       if(qresult(4).gt.qmax) then
                          !  so now we test to see if we have a smaller angle, which closes the box
                          qmax=qresult(4)
                          index1=i
                          index2=j
                          !               index3=0
                          index1a=0
                          if(ijk.eq.2) index1a=1
                          index2a=0
                          if(kjl.eq.2) index2a=1
                          !  record the result if we have found a quat in the fund. zone, and a smaller angle
                       end if
                    !end if - not using qresult(3)
                 end if
              end if
210        end do
185     end do
190  end do
200 end do

  if(index1.eq.0.or.index2.eq.0) stop 'DISQUAT_TET4onM failed!'

  call presymm(qq,index1,quint)

  call postsymm(quint,index2,qqn)
  !       write(*,"('DISQUAT: qqn = ',4f8.3)") qqn
  do 280, k=1,4
     if(index2a.eq.1) then
        !  if we had to take the negative of the Q, then do so
        !  to the answer

        qresult(k) = -1.0 * qqn(k)
        !  can take negative of Q because identity = +/-(0,0,0,1)
     else
        qresult(k)=qqn(k)
     endif
280 end do

  if(index1a.eq.1) then
     qresult(4) = -1. * qresult(4)
     !  restore the 4th component if disorient found for ijk=1
  endif

  write(*,*) 'DISQUAT_TET4onM: qmax= ',qmax
  write(*,"('DISQUAT_TET4onM: qresult = ',4f8.3)") qresult
  return
end subroutine disquat_tet4onM

!  ______________________

subroutine cslquat(qq,numsig,theta,qresult)

  implicit none

  double precision :: qq(4),theta,qresult(4),tmp(2)
  double precision :: quatcsl
  double precision, parameter :: pi  = 4. * atan ( 1.0 ) 
  integer :: numsig, nsig, numcsl ! index of CSL boundary - must be supplied

  common/a2/quatcsl(4,50),numcsl,nsig(50)

  !  algorithm for forming resultant quaternion
  !  and determining angle based on fourth component

  !  note change of signs to get inverse of second orientation

  qresult(1)=qq(1)*quatcsl(4,numsig)-qq(4)*quatcsl(1,numsig) +qq(2)*quatcsl(3,numsig)-qq(3)*quatcsl(2,numsig)
  qresult(2)=qq(2)*quatcsl(4,numsig)-qq(4)*quatcsl(2,numsig) +qq(3)*quatcsl(1,numsig)-qq(1)*quatcsl(3,numsig)
  qresult(3)=qq(3)*quatcsl(4,numsig)-qq(4)*quatcsl(3,numsig) +qq(1)*quatcsl(2,numsig)-qq(2)*quatcsl(1,numsig)
  qresult(4)=qq(4)*quatcsl(4,numsig)+qq(1)*quatcsl(1,numsig) +qq(2)*quatcsl(2,numsig)+qq(3)*quatcsl(3,numsig)

  !	write(*,*) 'qresult ',qresult
  if(qresult(4).gt.1.0) qresult(4)=1.0
  if(qresult(4).lt.-1.0) qresult(4)=-1.0
  theta=acos(qresult(4))*360./pi
  !	write(*,*) 'theta= ',theta

  return
end subroutine cslquat

! _____________________________

subroutine misquat(qq,thetamin)
  double precision qq(4,2),thetamin,qresult(4),tmp(2),rquat(4)
  double precision disor
  double precision, parameter :: pi  = 4. * atan ( 1.0 )
  double precision qmax,q1max,q2max

  !  algorithm for forming resultant quaternion
  !  and determining minimum angle taken from Sutton & Baluffi

  !  note that the resultant quaternion is not returned
  !  because it is not in the fundamental zone

  !  note change of signs to get inverse of second orientation

  qresult(1)=qq(1,1)*qq(4,2)-qq(4,1)*qq(1,2) +qq(2,1)*qq(3,2)-qq(3,1)*qq(2,2)
  qresult(2)=qq(2,1)*qq(4,2)-qq(4,1)*qq(2,2) +qq(3,1)*qq(1,2)-qq(1,1)*qq(3,2)
  qresult(3)=qq(3,1)*qq(4,2)-qq(4,1)*qq(3,2) +qq(1,1)*qq(2,2)-qq(2,1)*qq(1,2)
  qresult(4)=qq(4,1)*qq(4,2)+qq(1,1)*qq(1,2)+qq(2,1)*qq(2,2)+qq(3,1)*qq(3,2)

  !	write(*,*) 'qresult ',qresult
  qmax=0.
  iqindex=0
  do 10, i=1,4
     qresult(i)=abs(qresult(i))
     if(qresult(i).gt.qmax) then
        qmax=qresult(i)
        iqindex=i
     endif
10 end do

  q1max=0.
  iq1index=0
  !  find the next highest q component
  do 20, i=1,4
     if(i.eq.iqindex) goto 20
     if(qresult(i).gt.q1max) then
        q1max=qresult(i)
        iq1index=i
     endif
20 end do

  disor=amax1(qmax,(qmax+q1max)/sqrt(2.),(qresult(1)+qresult(2)+qresult(3)+qresult(4))/2.)
  if(disor.gt.1.0) disor=1.0
  if(disor.lt.-1.0) disor=-1.0
  thetamin=acos(disor)*360./pi
  !	write(*,*) 'thetamin ',thetamin

  q2max=0.
  iq2index=0
  !  find the next highest q component
  do 30, i=1,4
     if(i.eq.iqindex.or.i.eq.iq1index) goto 30
     if(qresult(i).gt.q2max) then
        q2max=qresult(i)
        iq2index=i
     endif
30 end do

  rquat(4)=qmax
  rquat(3)=q1max
  rquat(2)=q2max
  do 40, i=1,4
     if(i.ne.iqindex.and. i.ne.iq1index.and. i.ne.iq2index)  rquat(1)=qresult(i)
40 end do
  !  supply the sorted quaternion
  !  CAUTION:  note that q1<q2<q3<q4
  !  whereas typical Rodrigues sorting is R1>R2>R3

  return
end subroutine misquat

! -----------------------------------------

!  subroutine ZUR BERECHNUNG VON ORIENTIERUNGSBEZIEHUNGEN ZWISCHEN
!  ZWEI VORGEGEBENEN ORIENTIERUNGEN

!  subroutine provided by Dierk Raabe (to Paul Lee, 1996)

subroutine misori(p1,p,p2,al_min,nsymm,ac_min)

  implicit none

  !  p1,p,p2 are (2 sets of) Bunge Euler angles
  !  nsymm controls the sample symmetry
  !  al_min is the disorientation angle

  character meu,isy
  integer ac_min(3)
  integer i,j,k,i31max,nsymm,i31,i1,i24,i25,iz
  integer iomega,i_min,ii
  double precision D(3,3,2),SYM(3,3,24),P1(2),P(2),P2(2)
  double precision DM(3,3),DRR(3,3,24),A(3),IA(3)
  double precision DR(3,3)
  double precision tmp, tmp1, tmp2, tmp3 ,al_min ,fakt ,f
  double precision spur,sp,omega,alpha,x
  double precision c1,s1,pp1,pp,pp2,c,c2,s2,s
  double precision, parameter :: pi  = 4. * atan ( 1.0 ) 

  al_min=360.
  ac_min(1)=0
  ac_min(2)=0
  ac_min(3)=0
  !  ERSTELLEN DER SYMMETRIEMATRIZEN
  !  **************************************************************
  DO I=1,3
     DO J=1,3
        DO K=1,24
           SYM(I,J,K)=0.
        end do
     end do
  end do

  !   1
  SYM(1,1,1)=1.
  SYM(2,2,1)=1.
  SYM(3,3,1)=1.
  !   5
  SYM(1,1,2)=1.
  SYM(2,3,2)=-1.
  SYM(3,2,2)=1.
  !   2
  SYM(1,1,3)=1.
  SYM(2,2,3)=-1.
  SYM(3,3,3)=-1.
  !   11
  SYM(1,1,4)=1.
  SYM(2,3,4)=1.
  SYM(3,2,4)=-1.
  !   7
  SYM(1,3,5)=-1.
  SYM(2,2,5)=1.
  SYM(3,1,5)=1.
  !   12
  SYM(1,3,6)=1.
  SYM(2,2,6)=1.
  SYM(3,1,6)=-1.
  !   3
  SYM(1,1,7)=-1.
  SYM(2,2,7)=1.
  SYM(3,3,7)=-1.
  !   4
  SYM(1,1,8)=-1.
  SYM(2,2,8)=-1.
  SYM(3,3,8)=1.
  !   13
  SYM(1,2,9)=1.
  SYM(2,1,9)=-1.
  SYM(3,3,9)=1.
  !   6
  SYM(1,2,10)=-1.
  SYM(2,1,10)=1.
  SYM(3,3,10)=1.
  !   20
  SYM(1,2,11)=-1.
  SYM(2,3,11)=1.
  SYM(3,1,11)=-1.
  !   23
  SYM(1,3,12)=1.
  SYM(2,1,12)=-1.
  SYM(3,2,12)=-1.
  !   19
  SYM(1,2,13)=-1.
  SYM(2,3,13)=-1.
  SYM(3,1,13)=1.
  !   21
  SYM(1,3,14)=-1.
  SYM(2,1,14)=1.
  SYM(3,2,14)=-1.
  !   18
  SYM(1,2,15)=1.
  SYM(2,3,15)=-1.
  SYM(3,1,15)=-1.
  !   22
  SYM(1,3,16)=-1.
  SYM(2,1,16)=-1.
  SYM(3,2,16)=1.
  !   9
  SYM(1,2,17)=1.
  SYM(2,3,17)=1.
  SYM(3,1,17)=1.
  !   8
  SYM(1,3,18)=1.
  SYM(2,1,18)=1.
  SYM(3,2,18)=1.
  !   17
  SYM(1,2,19)=1.
  SYM(2,1,19)=1.
  SYM(3,3,19)=-1.
  !   24
  SYM(1,1,20)=-1.
  SYM(2,3,20)=1.
  SYM(3,2,20)=1.
  !   10
  SYM(1,3,21)=1.
  SYM(2,2,21)=-1.
  SYM(3,1,21)=1.
  !   15
  SYM(1,1,22)=-1.
  SYM(2,3,22)=-1.
  SYM(3,2,22)=-1.
  !   16
  SYM(1,3,23)=-1.
  SYM(2,2,23)=-1.
  SYM(3,1,23)=-1.
  !   14
  SYM(1,2,24)=-1.
  SYM(2,1,24)=-1.
  SYM(3,3,24)=-1.
  !  **********************************
  !  ERSTELLEN DER BEIDEN DMATRIZEN
  ! 205    FORMAT(A1)

  meu='E'
  isy='E'
  I31MAX=nsymm

  DO 31 I31=1,I31MAX
     DO 1 I1=1,2
        FAKT=PI/180.
        PP1=P1(I1)*FAKT
        PP=P(I1)*FAKT
        PP2=P2(I1)*FAKT
        C1=COS(PP1)
        C=COS(PP)
        C2=COS(PP2)
        S1=SIN(PP1)
        S=SIN(PP)
        S2=SIN(PP2)
        D(1,1,I1)=C1*C2-S1*S2*C
        D(1,2,I1)=S1*C2+C1*S2*C
        D(1,3,I1)=S2*S
        D(2,1,I1)=-C1*S2-S1*C2*C
        D(2,2,I1)=-S1*S2+C1*C2*C
        D(2,3,I1)=C2*S
        D(3,1,I1)=S1*S
        D(3,2,I1)=-C1*S
        D(3,3,I1)=C
        !         GOTO 1
1    END DO
     !  **********************************

     DO 24 I24=1,3
        DO 25 I25=1,3
           DM(I24,I25)=D(I24,I25,1)
           DR(I24,I25)=D(I24,I25,2)
25      END DO
24   END DO
     F=PI/180.

!!$     print*
!!$     print*,'1st matrix:'
!!$     do i24 = 1,3
!!$        print"('[',3(2x,f8.3),' ]')",(DM(I24,I25), i25 = 1,3)
!!$     enddo
!!$
!!$     print*
!!$     print*,'2nd matrix:'
!!$     do i24 = 1,3
!!$        print"('[',3(2x,f8.3),' ]')",(dr(I24,I25), i25 = 1,3)
!!$     enddo

     !  ****************************************************
     !  BESTIMMUNG DER INVERSEN MATRIX ZUR ORIENTIERUNG 1:DM
     !  ****************************************************

     DO I=1,3
        DO J=1,3
           DM(I,J)=D(J,I,1)
        end do
     end do

     !  *************************
     !  FELDER GLEICH NULL SETZEN
     !  *************************
     DO I=1,3
        DO J=1,3
           DO K=1,24
              DR(I,J)=0.
              DRR(I,J,K)=0.
           end do
        end do
     end do

     !  ***********************************************************
     !  MATRIZENMULTIPLIKATION DER MATRIZEN D(I,J,2) UND DM=DR(I,J)
     !  ***********************************************************

     DO 5 I=1,3
        DO J=1,3
           DO K=1,3
              DR(I,J)=D(I,K,2)*DM(K,J)+DR(I,J)
           end do
        end do
5    END DO
     !  DM is associated with the first orientation
     !   and has been transposed so that the misorientation is calculated properly
     !  thus, we have product matrix as gB X gA^-1

     !  ***********************************************
     !  MATRIZENMULITIPLIKATION DER MATRIZEN SYM UND DR
     !  ***********************************************

     DO 7 IZ=1,24
        DO 8 I=1,3
           DO J=1,3
              DO 9 K=1,3
                 DRR(I,J,IZ)=SYM(I,K,IZ)*DR(K,J)+DRR(I,J,IZ)
9             END DO
           END DO
8       END DO

!!!$	 print*
!!!$	 print*,' Symmetry operator number ',iz
!!!$	 print*,'Product matrix for gA X gB^-1:'
!!!$	 do i24 = 1,3
!!!$	    print"('[',3(2x,f8.3),' ]')",(drr(i24,i25,iz), i25 = 1,3)
!!!$	 enddo
!!!$	 tmp = (drr(1,1,iz)+drr(2,2,iz)+drr(3,3,iz))
!!!$	 print*,'Trace = ',tmp
!!!$	 tmp2 = (tmp-1)/2.
!!!$	 if(tmp2.lt.-1.) tmp2 = -1.
!!!$	 if(tmp2.gt.1.) tmp2 = 1.
!!!$	 tmp3 = acos(tmp2)*180./pi
!!!$	 print*,' angle = ',tmp3

        !  *******************************
        !  BESTIMMUNG DES ROTATIONSWINKELS
        !  *******************************
        SPUR=0.
        DO I=1,3
           SPUR=SPUR+DRR(I,I,IZ)
        end do
        SP=(SPUR-1.)/2.
        IF(SP.GE.1.)SP=1.-1.E-5
        IF(SP.LE.-1.)SP=-1.+1.E-5
        OMEGA=PI/2.-ASIN(SP)
        IOMEGA=INT(OMEGA*180./PI+0.5)
        ALPHA=OMEGA*180./PI
        IF(IOMEGA.EQ.180.or.iomega.eq.0)GOTO 12
        X=2.*SIN(OMEGA)

        !  ********************************************************
        !  BESTIMMUNG DER ROTATIONSACHSE, INCL SONDERFALL OMEGA = 180
        !  ********************************************************
        A(1)=100.*(DRR(2,3,IZ)-DRR(3,2,IZ))/X
        A(2)=100.*(DRR(3,1,IZ)-DRR(1,3,IZ))/X
        A(3)=100.*(DRR(1,2,IZ)-DRR(2,1,IZ))/X
        GOTO 13
12      A(1)=SQRT((DRR(1,1,IZ)+1.)/2.)*100.
        A(2)=SQRT((DRR(2,2,IZ)+1.)/2.)*100.
        A(3)=SQRT((DRR(3,3,IZ)+1.)/2.)*100.
13      DO 11 I=1,3
           IF(A(I).EQ.0.)GOTO 19
           IA(I)=INT(A(I)+0.5*A(I)/ABS(A(I)))
           GOTO 11
19         IA(I)=0
11      END DO

        if(abs(alpha).lt.al_min) then
           al_min=abs(alpha)
           ac_min(1)=ia(1)
           ac_min(2)=ia(2)
           ac_min(3)=ia(3)
        end if

7    END DO
     IF(I31.EQ.1)P1(2)=180.+P1(2)
     IF(I31.EQ.2)P1(2)=360.-P1(2)
     IF(I31.EQ.2)P2(2)=90.-P2(2)
     IF(I31.EQ.3)P1(2)=180.+P1(2)


31 END DO
  i_min=100
  if(abs(ac_min(1)).lt.i_min.and.abs(ac_min(1)).gt.0) i_min=abs(ac_min(1))
  if(abs(ac_min(2)).lt.i_min.and.abs(ac_min(2)).gt.0) i_min=abs(ac_min(2))
  if(abs(ac_min(3)).lt.i_min.and.abs(ac_min(3)).gt.0) i_min=abs(ac_min(3))

  do ii=1,3
     ac_min(ii)=nint(float(ac_min(ii))/float(i_min))
  end do

  if(p1(1).gt.360.) p1(1)=p1(1)-360.
  if(p1(1).lt.0.) p1(1)=p2(1)+360.
  if(p2(1).gt.360.) p2(1)=p1(1)-360.
  if(p2(1).lt.0.) p2(1)=p2(1)+360.
  if(p1(2).gt.360.) p1(2)=p1(2)-360.
  if(p1(2).lt.0.) p1(2)=p2(2)+360.
  if(p2(2).gt.360.) p2(2)=p1(2)-360.
  if(p2(2).lt.0.) p2(2)=p2(2)+360.

  return
end subroutine misori

!  -------------------

subroutine donorm(a)
  double precision a(3),b(3),rsum,rnorm
  rsum=0.

  do 100, i=1,3
     rsum=rsum+a(i)**2
100 enddo
  if (rsum .lt. 1.e-5 ) stop 'zero length axis?!'

  rnorm=sqrt(rsum)

  do 200, i=1,3
     b(i)=a(i)/rnorm
     a(i)=b(i)
200 enddo

  return
end subroutine donorm

!  ______________________

subroutine rfrvecR( rf , din , dout )

  implicit none

  !  uses RF-vector to rotate a vector from DIN to DOUT;
  !  Eq. 1.36 in Sutton & Balluffi
  !  REVERSE/NEGATIVE ROTATION; based on rfrvec;  7 vi 03

  double precision :: scalar
  double precision :: rf(3),din(3),dout(3)
  double precision :: tmp,tmp2,cross(3),rftmp(3)
  integer :: i

  do i=1,3
     rftmp(i)=-1.*rf(i)
  enddo
  tmp=2.*scalar(rftmp,din)
  tmp2=rftmp(1)**2+rftmp(2)**2+rftmp(3)**2
  call vecpro2(din,rftmp,cross)
  dout(1)=(((tmp*rftmp(1))+(din(1)*(1.-tmp2))) -(2.*cross(1)))/(1.+tmp2)
  dout(2)=(((tmp*rftmp(2))+(din(2)*(1.-tmp2))) -(2.*cross(2)))/(1.+tmp2)
  dout(3)=(((tmp*rftmp(3))+(din(3)*(1.-tmp2))) -(2.*cross(3)))/(1.+tmp2)
  return
end subroutine rfrvecR

!     ______________________

subroutine q2eulB ( p1 , p , p2 , qq )

  implicit none
  !     converts quaternion to Bunge Euler angles
  !     based on Altmann's solution for Euler->quat

  double precision qq(4) , p1 , p , p2
  double precision sum , diff , tmp
  double precision, parameter :: pi  = 4. * atan ( 1.0 ) 

  !  CODE:
  if((abs(qq(2)).lt.1e-35).and.(abs(qq(1)).lt.1e-35)) then
     diff=pi/4.
  else
     diff=atan2(qq(2),qq(1))
  endif
  !     ATAN2 is unhappy is both arguments are exactly zero!
  !     
  if((abs(qq(3)).lt.1e-35).and.(abs(qq(4)).lt.1e-35)) then
     sum=pi/4.
  else
     sum=atan2(qq(3),qq(4))
  endif
  !     
  p1 = (diff+sum)
  p2 = (sum-diff)
  tmp = sqrt(qq(3)**2+qq(4)**2)
  if ( tmp > 1.0 ) tmp = 1.0
  if ( tmp < -1.0 ) tmp = -1.0
  p = 2.*acos(tmp)
  !     write(*,*) ' quaternion input= ',qq
  !     write(*,*) 'Bunge angles output= ',p1,p,p2
  !     write(*,*) ' angles [degrees]= ',180.*p1/pi,180.*p/pi,180.*p2/pi
  return
end subroutine q2eulB

!  --------------------

function scalar(a,b)

  implicit none
  !  return SCALAR PRODUCT of A and B
  integer :: i
  double precision :: scalar,a(3),b(3)
  scalar=0.d0
  do 100, i=1,3
     scalar = scalar+a(i)*b(i)
100 end do
  return
end function scalar

! ____________________________

subroutine vecpro2(v1,v2,vout)

  implicit none
  !  vector product
  double precision v1(3),v2(3),vout(3)
  vout(3)=v1(1)*v2(2)-v1(2)*v2(1)
  vout(1)=v1(2)*v2(3)-v1(3)*v2(2)
  vout(2)=v1(3)*v2(1)-v1(1)*v2(3)
  return
end subroutine vecpro2

! ________________________________________


subroutine textur2( afile , iq_data_type )

  !  second version, for use in REXGBS, WTS2POP
  !  not for use in REX2D
  use var
  use strings

  implicit none

  character, intent(IN) :: afile*200
  !integer, intent(IN) :: iq_data_type , iq2
  integer, intent(IN) :: iq_data_type  !  as of 15 i 16
  integer :: i,j,k,ijk ,ng
  double precision :: rad ,d1 ,d2 ,d3 , d4 , d5 , d6 , rtmp

  double precision :: aquat(4,norients),phi1(norients),capphi(norients),phi2(norients),grwt(norients)
  integer :: aspin(norients) , ngrain
  common/a3/aquat ,phi1 ,capphi ,phi2 ,grwt , aspin , ngrain

!!!$	common/a3/aquat(4,norients),ngrain
!!!$     &	,phi1(norients),capphi(norients),phi2(norients),grwt(norients)
  !	common/a4/ textitl,eulan,iper

  !	include 'common.f'
  !  following code taken from LApp68.f for reading in TEXIN

  double precision :: evm,fk1(9),nstate
  character*80 :: textitl
  character :: eulan*25,nomen,iper*6
  double precision :: ph,th,ps,dph,dth,dps,sph,cph,cth,sth,sps,cps,qr(4)
  double precision :: theta,r1,r2,r3
  integer :: lmn
  integer :: tmp1 , tmp2 , tmp3
  character*120 :: inputline

  integer :: ios  !  for I/O errors
  integer :: nargs
  character(len=100) , dimension(12) :: args
  
!  CODE::
  rad = 57.29578
  !  symmetry now read in by textur4

!!$  write(*,*) 'reading CSL information from "quat.csl.data"'
!!$  open(unit=3,file = 'quat.csl.data',status='old')
!!$  read(3,*)               ! skip over title line
!!$  read(3,*) numcsl ! read number of CSL types
!!$  if(numcsl.gt.50) then
!!$     numcsl=50
!!$     write(*,*) 'too many CSL types, cutting back to 50!'
!!$  else
!!$     write(*,*) 'reading ',numcsl,' boundary types'
!!$  endif
!!$  read(3,*)               ! skip over comment
!!$  do 1900, i=1,numcsl
!!$     read(3,*) nsig(i),theta,lmn,r1,r2,r3,(quatcsl(ijk,i),ijk=1,4)
!!$     !     write(*,*) i,nsig(i),theta,lmn,r1,r2,r3,(quatcsl(ijk,i),ijk=1,4)
!!$1900 end do
!!$  close(unit=3)

  write(*,*) 'TEXTUR2: reading grain orientations from ', trim(afile)
  print*, 'TEXTUR2: reading with iq_data_type = ', iq_data_type
  !     open(unit=3,file = 'texin',status='old')
  open( unit=3 , file=afile , status='old' )

  if( iq_data_type == 1 ) then  !  TEXIN type inputs
     !print*, 'TEXTUR2: reading with iq_data_type = ', iq_data_type
     read( 3 , 7 ) textitl
     print*,'TEXTUR2: Header Line: ',textitl
7    format(a)
     read( 3 , * ) !  ignore header line
     read( 3 , * ) !  ignore these values: evm,fk1,nstate
70   format(10f7.3,i3)
     read(3,71) eulan,iper
71   format(a25,49x,a6)
     nomen=eulan(1:1)
     if(nomen.eq.'B') then
        write(*,*) 'detected Bunge angles'
     elseif(nomen.eq.'K') then
        write(*,*) 'detected Kocks angles'
     else
        print*,'Unknown Euler angle type, quitting'
        call exit(1)
     end if
  elseif( iq_data_type == 2 ) then  !  Sukbin lists
     read( 3 , 7 ) textitl
     print*,'TEXTUR2: Header Line: ',textitl
     read(3,*)        !  discard 1st zero quat
  elseif( iq_data_type == 0 ) then  !  text file with 1 header line
     print*, 'TEXTUR2: CAREFUL!  a file with Mac line endings may not read in properly or fully!!'
     read( 3 , 7 ) textitl
     print*,'TEXTUR2: Header Line: ',textitl
  end if

  if( iq_data_type == 2 .OR. iq_data_type == 1 ) then
     do 529 ng = 1 , norients

        if( iq_data_type == 1 ) then  !  TEXIN type inputs
           !     do 520 i=1,18
           !     520  state(i,ng)=0.
           !     read(3,50,err=530,end=530) d1,d2,d3,grwt(ng),(state(i,ng),i=1,nstate)
           read( 3 , * , end=530 , err=530 ) d1 , d2 , d3 , grwt(ng)
           aspin( ng ) = 0  !  no info!
           write(*,*) 'angles..  ',d1,d2,d3
           !     if((d1.eq.999.).and.(d2.eq.999.).and.(d3.eq.999.)) goto 530
50         format(6f8.2,24f6.2)
           !     read(3,*,end=530) d1,d2,d3,grwt(ng)
           !     .: for old unformatted TEXIN files

           if(nomen.eq.'B') then
              phi1(ng)=d1
              capphi(ng)=d2
              phi2(ng)=d3
           elseif(nomen.eq.'K') then
              phi1(ng)=d1+90.
              capphi(ng)=d2
              phi2(ng)=90.-d3
           end if
           if(phi1(ng).gt.360.) phi1(ng) = phi1(ng)-360.
           if(phi1(ng).lt.0.) phi1(ng) = phi1(ng)+360.
           if(phi2(ng).gt.360.) phi2(ng) = phi2(ng)-360.
           if(phi2(ng).lt.0.) phi2(ng)=phi2(ng)+360.

20         dth=d2
           dps=d1

           call  quatB(phi1(ng)/rad,capphi(ng)/rad,phi2(ng)/rad,qr)

           !     use the new conversion routine based on Altmann if Bunge
           !     else, use longer routine

        elseif( iq_data_type == 2 ) then
           read( 3 , * , end=530 ) tmp1 , tmp2 , tmp3 , ( qr(j) , j = 1 , 4 )
           !  note that Sukbin writes out the sequence number (original spin, I think)
           !    and then the new, randomized spin number, then the volume, then the quat

           aspin( ng ) = tmp1
           call q2eulB ( phi1(ng) , capphi(ng) , phi2(ng) , qr )
           !$$$  print*
           !$$$  print"(4(i7),' q: ',4(1x,f8.3))"
           !$$$  &           ,ng,tmp1,tmp2,tmp3,(qr(j),j=1,4)
           !     if( tmp2.lt.1 .or. tmp2.gt.qvalue ) then
           !if( tmp2 < llimit .or. tmp2 > qvalue ) then
           !   print*,'warning for line ', tmp1 ,', spin number out of range '
           ! call exit(1)
           !else
           !call qnorm(qr)
           !     3 integers (sequence_no, spin_number, volume) then quaternion

           !do  j = 1,4
           !   aquat( j , tmp2 ) = qr( j )
           !     note use of tmp2 as index
           !end do
           !     print*,'no, quat..  ',tmp2,(aquat(j,tmp2),j=1,4)

        end if

        call qnorm(qr)
        do j=1,4
           aquat(j,ng)=qr(j)
        end do
        ngrain = ng

529  end do
530  continue

  end if  !  iq_data_type = 1 or 2

  if( iq_data_type == 0 ) then
     do ng = 1 , 99999 , 2
        i = ng
        !read( 3 , "(a)" , end=1530 ) inputline
        call readline( 3 , inputline , ios )
        if ( ios /= 0 ) exit
        print*, 'inputline: ', inputline
        call compact( inputline )
        call parse( inputline , ' ' , args , nargs )
        if ( ng == 1 .AND. nargs <= 6 ) then
           print*,' TEXTUR2: only 6 numbers on the input line, no WEIGHT input'
        elseif ( ng == 1 .AND. nargs < 6 ) then
           print*,' TEXTUR2: not enough numbers (<6) on the input line, stopping'
           call exit(1)
        else
           call value( args(7) , rtmp , ios )
        end if
        call value( args(1) , d1 , ios )
        call value( args(2) , d2 , ios )
        call value( args(3) , d3 , ios )
        call value( args(4) , d4 , ios )
        call value( args(5) , d5 , ios )
        call value( args(6) , d6 , ios )
        !read( inputline , * ) d1 , d2 , d3 , d4 , d5 , d6 , rtmp
        !read( 3 , "(7f10.3)" , end=1530 ) d1 , d2 , d3 , d4 , d5 , d6 , rtmp
        ngrain = ng + 1  !  corrected 26 xi 15, but check old MacBook
        !print * , i, ngrain, d1 , d2 , d3 , d4 , d5 , d6
        !ngrain = (ng+1)/2
        phi1( i )   = d1  !  assumed Bunge
        capphi( i ) = d2
        phi2( i )   = d3
        phi1( i+1 )   = d4
        capphi( i+1 ) = d5
        phi2( i+1 )   = d6
        grwt( i ) = rtmp
        grwt( i+1 ) = rtmp
        call  quatB(phi1( i )/rad,capphi( i )/rad,phi2( i )/rad,qr)
        call qnorm(qr)
        do j=1,4
           aquat( j , i ) = qr(j)
        end do
        call  quatB(phi1( i+1 )/rad,capphi( i+1 )/rad,phi2( i+1 )/rad,qr)
        call qnorm(qr)
        do j=1,4
           aquat( j , i+1 ) = qr(j)
        end do

     end do
  end if   !   iq_data_type == 0

1530 continue
  print * , 'TEXTUR2:  number of grains = ngrain = ' , ngrain
!!$  do i = 1 , ngrain
!!$     print*,'i, anglesA ',i , phi1( i ) , capphi( i ) , phi2( i )
!!$  end do
  
  return
end subroutine textur2

! ________________________________________


subroutine textur3( bfile , iq_data_type )

  !  reads in the 2nd data file
  use var
  implicit none

  integer :: iq_data_type
  integer i,j,k,ijk ,ng
  double precision rad ,d1 ,d2 ,d3

  double precision :: bquat(4,norients),bphi1(norients),bcapphi(norients),bphi2(norients),bgrwt(norients)
  integer :: bspin(norients) , bngrain
  common/a4/ bquat , bphi1 , bcapphi , bphi2 , bgrwt , bspin , bngrain

  !  following code taken from LApp68.f for reading in TEXIN

  double precision evm,fk1(9),nstate
  character*80 textitl
  character eulan*25,nomen,iper*6
  double precision ph,th,ps,dph,dth,dps,sph,cph,cth,sth,sps,cps,qr(4)
  double precision theta,r1,r2,r3
  integer lmn
  character bfile*200

  integer :: tmp1 , tmp2 , tmp3

!  CODE::
  rad=57.29578
  !  symmetry now read in by textur4

  write(*,*) 'reading 2nd set of grain orientations from ',bfile
  !     open(unit=3,file = 'texin',status='old')
  open(unit=3, file=bfile, status='old')

  if( iq_data_type == 1 ) then  !  TEXIN type inputs
     read(3,7) textitl
     print*,'Header Line: ',textitl
7    format(a)
     read ( 3 , * )
     read ( 3 , * ) !  ignore these: evm,fk1,nstate
70   format(10f7.3,i3)
     read( 3 , "(a)" ) eulan !  ignore these: iper
71   format(a25,49x,a6)
     nomen=eulan(1:1)
     if(nomen.eq.'B') then
        write(*,*) 'detected Bunge angles'
     elseif(nomen.eq.'K') then
        write(*,*) 'detected Kocks angles'
     else
        print*,'Unknown Euler angle type, quitting'
        call exit(1)
     end if
  elseif( iq_data_type == 2 ) then  !  Sukbin lists
     read(3,*)        !  1 header line
     read(3,*)        !  discard 1st zero quat
  end if

  do 529 ng=1,norients

     if( iq_data_type == 1 ) then  !  TEXIN type inputs
        !     do 520 i=1,18
        !     520  state(i,ng)=0.
        !     read(3,50,err=530,end=530) d1,d2,d3,grwt(ng),(state(i,ng),i=1,nstate)
        read(3,*,end=530,err=530) d1,d2,d3,bgrwt(ng)
        bspin( ng ) = 0  !  no info!
        write(*,*) 'angles..  ',d1,d2,d3
        !     if((d1.eq.999.).and.(d2.eq.999.).and.(d3.eq.999.)) goto 530
50      format(6f8.2,24f6.2)
        !     read(3,*,end=530) d1,d2,d3,grwt(ng)
        !     .: for old unformatted TEXIN files

        if(nomen.eq.'B') then
           bphi1(ng) = d1
           bcapphi(ng) = d2
           bphi2(ng) = d3
        elseif(nomen.eq.'K') then
           bphi1(ng) = d1 + 90.
           bcapphi(ng) = d2
           bphi2(ng) = 90. - d3
        endif
        if(bphi1(ng).gt.360.) bphi1(ng) = bphi1(ng) - 360.
        if(bphi1(ng).lt.0.  ) bphi1(ng) = bphi1(ng) + 360.
        if(bphi2(ng).gt.360.) bphi2(ng) = bphi2(ng) - 360.
        if(bphi2(ng).lt.0.  ) bphi2(ng) = bphi2(ng) + 360.

20      dth=d2
        dps=d1

        call  quatB(bphi1(ng)/rad , bcapphi(ng)/rad , bphi2(ng)/rad , qr )

        !     use the new conversion routine based on Altmann if Bunge
        !     else, use longer routine

     elseif( iq_data_type == 2 ) then
        read( 3 , * , end = 530 ) tmp1 , tmp2 , tmp3 , ( qr(j) , j = 1 , 4 )
!  note that SUkbin writes out the sequence number (original spin, I think)
!    and then the new, randomized spin number, then the volume, then the quat
        bspin( ng ) = tmp1  !  no info!
        call q2eulB ( bphi1(ng) , bcapphi(ng) , bphi2(ng) , qr )
        !$$$  print*
        !$$$  print"(4(i7),' q: ',4(1x,f8.3))"
        !$$$  &           ,ng,tmp1,tmp2,tmp3,(qr(j),j=1,4)
        !     if( tmp2.lt.1 .or. tmp2.gt.qvalue ) then
        !if( tmp2 < llimit .or. tmp2 > qvalue ) then
        !   print*,'warning for line ', tmp1 ,', spin number out of range '
           ! call exit(1)
        !else
        !call qnorm(qr)
        !     3 integers (sequence_no, spin_number, volume) then quaternion

        !do  j = 1,4
        !   aquat( j , tmp2 ) = qr( j )
        !     note use of tmp2 as index
        !end do
     end if

     call qnorm(qr)
     do j=1,4
        bquat(j,ng) = qr(j)
     end do

529 end do
530 continue
  !bngrain = ng - 1
  bngrain = ng
  print * , 'TEXTUR3:  number of grains, 2nd list = bngrain = ' , bngrain

  return
end subroutine textur3

! ________________________________________

subroutine textur4(iq_symm)
  ! 1 x 15, not needed any longer
  
  !  second version, for use in REXGBS, WTS2POP
  !  not for use in REX2D

  integer numsymm,numsamp
  integer iq_symm,n_max
  double precision quatsymm(4,48),quatsamp(4,4)
  common/a1/ quatsymm,numsymm,quatsamp,numsamp

  !  CODE::
  if(iq_symm < 0 .or. iq_symm > 4 ) stop 'bad symmetry class?'

  if( iq_symm == 0 ) then
     n_max = 24
  elseif( iq_symm == 1 ) then
     n_max = 12  !  Hex
  elseif( iq_symm == 2 ) then
     n_max = 4   !  Orth
  elseif( iq_symm == 3 ) then
     n_max = 8   !  Tetr  4/mmm
  elseif( iq_symm == 4 ) then
     n_max = 4   !  Tetr, 4/m
  end if
  !  confusing - wts2pop uses 1 to 3
  !     print*, 'reading symmetry operators from "quat.symm.cubic"'
  if( iq_symm == 0 ) then
     open( unit=3,file = 'quat.symm.cubic',status='old')
  elseif( iq_symm == 1 ) then
     open( unit=3,file = 'quat.symm.hex',status='old')
  elseif( iq_symm == 2 ) then
     open( unit=3,file = 'quat.symm.orth',status='old')
  elseif( iq_symm == 3 ) then
     open( unit=3,file = 'quat.symm.tet',status='old')
  elseif( iq_symm == 4 ) then
     open( unit=3,file = 'quat.symm.tet4onM',status='old')
  end if

  read(3,*)   ! skip over title line
  read(3,*) numsymm  ! read number of symmetry operators
  if( numsymm > n_max ) then
     numsymm = n_max
     print*, 'too many symmetry operators, cutting back to ',n_max
  else
     !     print*, 'reading ',numsymm,' symmetry operators'
  endif

  do 1800, i=1,numsymm
     read(3,*) (quatsymm(ijk,i),ijk=1,4)
     !     print*, 'symm. no. ',i,(quatsymm(ijk,i),ijk=1,4)
1800 end do
  close(unit=3)

  !  now we hard-wire the four operators into the code
  !  so that choices by user are easier to deal with
  numsamp = 4
  do i = 1,numsamp
     do ijk = 1,4
        quatsamp(ijk,i) = 0.
     enddo
  enddo
  quatsamp(4,1) = 1.
  quatsamp(1,2) = 1.
  quatsamp(2,3) = 1.
  quatsamp(3,4) = 1.

  return
end subroutine textur4

! ________________________________________
!  1 x 15, not needed any longer
subroutine textur6
  use var
  implicit none

  !  second version, for use in REXGBS, WTS2POP
  !  not for use in REX2D
  !  edited because textur4 now reads in symmetry operators

  !integer numsymm,numsamp
  integer :: i,ijk , lmn
  !double precision quatsymm(4,48),quatsamp(4,4)
  !double precision quatcsl(4,50),
  double precision :: theta,r1,r2,r3
  !integer numcsl,nsig(50),lmn

  !common/a2/quatcsl,numcsl,nsig

  !  CODE::
  write(*,*) 'reading CSL information from "quat.csl.data"'
  open(unit=3,file = 'quat.csl.data',status='old')
  read(3,*) ! skip over title line
  read(3,*) numcsl ! read number of CSL types
  if(numcsl.gt.50) then
     numcsl=50
     write(*,*) 'too many CSL types, cutting back to 50!'
  else
     write(*,*) 'reading ',numcsl,' boundary types'
  endif
  read(3,*) ! skip over comment
  do 1900, i=1,numcsl
     read(3,*) nsig(i),theta,lmn,r1,r2,r3,(quatcsl(ijk,i),ijk=1,4)
     !		write(*,*) i,nsig(i),theta,lmn,r1,r2,r3,(quatcsl(ijk,i),ijk=1,4)
1900 end do
  close(unit=3)

  return
end subroutine textur6


!  -------------------

subroutine euler(a,iopt,nomen,d1,d2,d3,ior,kerr)

  implicit none

  !                          Last revision 20nov90 UFK
  !      common a(3,3),grvol(1152),epsga(5),ist1,ist2,sqr3,sqrh,ident(3,3)
  !  SPECIAL VERSION WITHOUT COMMON BLOCK
  integer :: iopt, ior, kerr, kor
  double precision :: a(3,3),d1,d2,d3,th,sth,cth,sph,cph,sps,cps ,ps,ph,dth,dph,dps
  double precision :: rad
  character nomen
  double precision, parameter :: pi  = 4. * atan ( 1.0 ) 

  save kor

  !  CODE:
  rad = 57.29578

  ! *** this subroutine calculates the euler angles associated with the
  !     transformation matrix a(i,j) if iopt=1 and viceversa if iopt=2
  ! ***  Note that a is sample (rows) in terms of crystal (columns);
  !      -- opposite of standard g (e.g.Bunge) - this is Canova's
  ! ***  Note that in this version, the Euler angles are defined symmetrically:
  !      so that interchanging phi and psi means transposing a.
  !      ("Kocks" nomen: defined going from +X to +Y in both COD and SOD)
  ! ***  However, other angle conventions are translated, according to
  !  nomen="K" - Kocks (as internally) -- also sometimes "N"...
  !        "R" - Roe          (Psi=psi, Phi=180-phi)
  !        "B" - Bunge        (phi1=90+psi, Phi=Theta, phi2=90-phi)
  !        any other - Canova (Theta first, phiC=90+phi, omega=90-psi)
  ! ***   Note: only in symm. notation does a point with all Euler angles
  !             between 0 and 90 deg appear in the +x,+y quadrant!
  !       If you want to see an individual point in this quadrant and are:
  !       using Roe,    the third angle must be between 90 and 180;
  !             Bunge,      first                                 ;
  !             Canova,     second                                .
  ! ***  Input and output Euler angles d1,d2,d3 in degrees

  goto(5,20),iopt
5 th = acos(a(3,3))
  sth = sin(th)
  if( abs(a(3,3)) .ge. 0.999999 ) goto 10

  if((abs(a(2,3)/sth).lt.1e-35).and. (abs(a(1,3)/sth).lt.1e-35)) then
     ps=pi/4.
  else
     ps=atan2(a(2,3)/sth,a(1,3)/sth)
  endif

  if((abs(a(3,2)/sth).lt.1e-35).and. (abs(a(3,1)/sth).lt.1e-35)) then
     ph=pi/4.
  else
     ph=atan2(a(3,2)/sth,a(3,1)/sth)
  endif
  !   if it bombs out here, both a-components in one arg. are zero
  !       (this should not be possible, but has happened, probably fixed)
  !  ADR: i 01 - above is the fix!!

  go to 15

10 if((abs(a(1,2)/sth).lt.1e-35).and. (abs(a(1,1)/sth).lt.1e-35)) then
     ps=pi/4.
  else
     ps=0.5*atan2(a(1,2),-a(1,1))
  endif

  ph=-ps
  !       The above still have the problem that they give too many
  !       equivalents of the same grain in DIOROUT and density file. Therefore:
  if(kerr.eq.1.and.kor.ne.ior) then
     print*,'NOTE: grain ',ior,' has Theta < 1 deg.: sometimes problems'
     print*
     kor=ior
  endif

15 dth=th*rad
  dph=ph*rad
  dps=ps*rad
  d1=dps
  d2=dth
  if(nomen.eq.'K'.or.nomen.eq.'N') then
     d3=dph
  elseif(nomen.eq.'R') then
     d3=180.-dph
  elseif(nomen.eq.'B') then
     d1=dps+90.
     d3=90.-dph
  else
     d1=dth
     d2=dph+90.
     d3=90.-dps
  endif
  if(d1.ge.360.) d1=d1-360.
  if(d3.ge.360.) d3=d3-360.
  if(d1.lt.0.) d1=d1+360.
  if(d3.lt.0.) d3=d3+360.
  return
  !  *************************************
20 dth=d2
  dps=d1
  if(nomen.eq.'K') then
     dph=d3
  elseif(nomen.eq.'R') then
     dph=180.-d3
  elseif(nomen.eq.'B') then
     dps=d1-90.
     dph=90.-d3
  else
     dth=d1
     dph=d2-90.
     dps=90.-d3
  endif
  ph=dph/rad
  th=dth/rad
  ps=dps/rad
  sph=sin(ph)
  cph=cos(ph)
  sth=sin(th)
  cth=cos(th)
  sps=sin(ps)
  cps=cos(ps)
  a(1,1)=-sps*sph-cph*cps*cth
  a(2,1)=cps*sph-cph*sps*cth
  a(3,1)=cph*sth
  a(1,2)=sps*cph-sph*cps*cth
  a(2,2)=-cph*cps-sph*sps*cth
  a(3,2)=sth*sph
  a(1,3)=sth*cps
  a(2,3)=sps*sth
  a(3,3)=cth
  !  note that this is the transpose of the standard definition
  !  because of the way in which Canova used the matrix
  return
end subroutine euler


! -----------------------------------------

!  subroutine ZUR BERECHNUNG VON ORIENTIERUNGSBEZIEHUNGEN ZWISCHEN
!  ZWEI VORGEGEBENEN ORIENTIERUNGEN

!  subroutine provided by Dierk Raabe (to Paul Lee, 1996)

subroutine misorinv(p1,p,p2,al_min,nsymm,ac_min)

  !  p1,p,p2 are (2 sets of) Bunge Euler angles
  !  nsymm controls the sample symmetry
  !  al_min is the disorientation angle

  character meu,isy
  integer ac_min(3)
  double precision :: D(3,3,2),SYM(3,3,24),P1(2),P(2),P2(2),DM(3,3),DRR(3,3,24),A(3),IA(3),DR(3,3)
  double precision :: tmp, tmp1, tmp2, tmp3
  double precision :: PP, PP1, PP2, F, FAKT, SPUR, OMEGA, ALPHA, C, S, C1, C2, S1, S2, SP, X, al_min
  double precision, parameter :: pi  = 4. * atan ( 1.0 ) 

  al_min=360.
  ac_min(1)=0
  ac_min(2)=0
  ac_min(3)=0
  !  ERSTELLEN DER SYMMETRIEMATRIZEN
  !  **************************************************************
  DO I=1,3
     DO J=1,3
        DO K=1,24
           SYM(I,J,K)=0.
        end do
     end do
  end do
  !   1
  SYM(1,1,1)=1.
  SYM(2,2,1)=1.
  SYM(3,3,1)=1.
  !   5
  SYM(1,1,2)=1.
  SYM(2,3,2)=-1.
  SYM(3,2,2)=1.
  !   2
  SYM(1,1,3)=1.
  SYM(2,2,3)=-1.
  SYM(3,3,3)=-1.
  !   11
  SYM(1,1,4)=1.
  SYM(2,3,4)=1.
  SYM(3,2,4)=-1.
  !   7
  SYM(1,3,5)=-1.
  SYM(2,2,5)=1.
  SYM(3,1,5)=1.
  !   12
  SYM(1,3,6)=1.
  SYM(2,2,6)=1.
  SYM(3,1,6)=-1.
  !   3
  SYM(1,1,7)=-1.
  SYM(2,2,7)=1.
  SYM(3,3,7)=-1.
  !   4
  SYM(1,1,8)=-1.
  SYM(2,2,8)=-1.
  SYM(3,3,8)=1.
  !   13
  SYM(1,2,9)=1.
  SYM(2,1,9)=-1.
  SYM(3,3,9)=1.
  !   6
  SYM(1,2,10)=-1.
  SYM(2,1,10)=1.
  SYM(3,3,10)=1.
  !   20
  SYM(1,2,11)=-1.
  SYM(2,3,11)=1.
  SYM(3,1,11)=-1.
  !   23
  SYM(1,3,12)=1.
  SYM(2,1,12)=-1.
  SYM(3,2,12)=-1.
  !   19
  SYM(1,2,13)=-1.
  SYM(2,3,13)=-1.
  SYM(3,1,13)=1.
  !   21
  SYM(1,3,14)=-1.
  SYM(2,1,14)=1.
  SYM(3,2,14)=-1.
  !   18
  SYM(1,2,15)=1.
  SYM(2,3,15)=-1.
  SYM(3,1,15)=-1.
  !   22
  SYM(1,3,16)=-1.
  SYM(2,1,16)=-1.
  SYM(3,2,16)=1.
  !   9
  SYM(1,2,17)=1.
  SYM(2,3,17)=1.
  SYM(3,1,17)=1.
  !   8
  SYM(1,3,18)=1.
  SYM(2,1,18)=1.
  SYM(3,2,18)=1.
  !   17
  SYM(1,2,19)=1.
  SYM(2,1,19)=1.
  SYM(3,3,19)=-1.
  !   24
  SYM(1,1,20)=-1.
  SYM(2,3,20)=1.
  SYM(3,2,20)=1.
  !   10
  SYM(1,3,21)=1.
  SYM(2,2,21)=-1.
  SYM(3,1,21)=1.
  !   15
  SYM(1,1,22)=-1.
  SYM(2,3,22)=-1.
  SYM(3,2,22)=-1.
  !   16
  SYM(1,3,23)=-1.
  SYM(2,2,23)=-1.
  SYM(3,1,23)=-1.
  !   14
  SYM(1,2,24)=-1.
  SYM(2,1,24)=-1.
  SYM(3,3,24)=-1.
  !  **********************************
  !  ERSTELLEN DER BEIDEN DMATRIZEN
  ! 205    FORMAT(A1)

  meu='E'
  isy='E'
  I31MAX=nsymm

  DO 31 I31=1,I31MAX
     DO 1 I1=1,2
        FAKT=PI/180.
        PP1=P1(I1)*FAKT
        PP=P(I1)*FAKT
        PP2=P2(I1)*FAKT
        C1=COS(PP1)
        C=COS(PP)
        C2=COS(PP2)
        S1=SIN(PP1)
        S=SIN(PP)
        S2=SIN(PP2)
        D(1,1,I1)=C1*C2-S1*S2*C
        D(1,2,I1)=S1*C2+C1*S2*C
        D(1,3,I1)=S2*S
        D(2,1,I1)=-C1*S2-S1*C2*C
        D(2,2,I1)=-S1*S2+C1*C2*C
        D(2,3,I1)=C2*S
        D(3,1,I1)=S1*S
        D(3,2,I1)=-C1*S
        D(3,3,I1)=C
        GOTO 1
1    END DO
     !  **********************************

     !       DO 24 I24=1,3
     !	  DO 25 I25=1,3
     !       DM(I24,I25)=D(I24,I25,1)
     !	     DR(I24,I25)=D(I24,I25,2)
     ! 25	  END DO
     !24     END DO
     F=PI/180.

     !  ****************************************************
     !  BESTIMMUNG DER INVERSEN MATRIX ZUR ORIENTIERUNG 1:DM
     !  ****************************************************

     DO 3 I=1,3
        DO J=1,3
           DM(I,J)=D(J,I,1)
        end do
3    end do

     !  *************************
     !  FELDER GLEICH NULL SETZEN
     !  *************************

     DO 4 I=1,3
        DO J=1,3
           DO K=1,24
              DR(I,J)=0.
              DRR(I,J,K)=0.
           end do
        end do
4    end do

     !  ***********************************************************
     !  MATRIZENMULTIPLIKATION DER MATRIZEN D(I,J,2) UND DM=DR(I,J)
     !  ***********************************************************

     !       DO 5 I=1,3
     !       DO 5 J=1,3
     !       DO 6 K=1,3
     !6      DR(I,J)=D(I,K,2)*DM(K,J)+DR(I,J)
     !5      END DO

     !  ***********************************************
     !  MATRIZENMULITIPLIKATION DER MATRIZEN SYM UND DR
     !  ***********************************************

     DO 7 IZ=1,24

        DO I=1,3
           DO J=1,3
              dr(i,j) = 0.  !  just to make sure
              DO K=1,3
                 !		   DR(I,J)=D(i,k,2)*SYM(k,j,IZ) + DR(I,J)
                 DR(I,J) = SYM(i,k,IZ)*D(k,j,2) + DR(I,J)
                 !  apply symmetry before combining
              enddo
           enddo
        enddo
        !8      END DO

        DO I=1,3
           DO J=1,3
              drr(i,j,iz) = 0.  !  just to make sure
              DO K=1,3
                 drr(I,J,iz)=dr(k,i)*DM(j,k)+drr(I,J,iz)
                 !  and use transpose of DM, because it is already transposed
                 !  so as to have the product matrix as gB X gA^-1 (but both transposed!)
              enddo
           enddo
        enddo
        !5      END DO

!!!$	 print*
!!!$	 print*,' Symmetry operator number ',iz
!!!$	 print*,'Product matrix for gB X gA^-1:'
!!!$	 do i24 = 1,3
!!!$	    print"('[',3(2x,f8.3),' ]')",(drr(i24,i25,iz), i25 = 1,3)
!!!$	 enddo
!!!$	 tmp = (drr(1,1,iz)+drr(2,2,iz)+drr(3,3,iz))
!!!$	 print*,'Trace = ',tmp
!!!$	 tmp2 = (tmp-1)/2.
!!!$	 if(tmp2.lt.-1.) tmp2 = -1.
!!!$	 if(tmp2.gt.1.) tmp2 = 1.
!!!$	 tmp3 = acos(tmp2)*180./pi
!!!$	 print*,' angle = ',tmp3

        !  *******************************
        !  BESTIMMUNG DES ROTATIONSWINKELS
        !  *******************************
        SPUR=0.
        DO 10 I=1,3
           SPUR=SPUR+DRR(I,I,IZ)
10      end do
        SP=(SPUR-1.)/2.
        IF(SP.GE.1.)SP=1.-1.E-5
        IF(SP.LE.-1.)SP=-1.+1.E-5
        OMEGA=PI/2.-ASIN(SP)
        IOMEGA=INT(OMEGA*180./PI+0.5)
        ALPHA=OMEGA*180./PI
        IF(IOMEGA.EQ.180.or.iomega.eq.0)GOTO 12
        X=2.*SIN(OMEGA)

        !  ********************************************************
        !  BESTIMMUNG DER ROTATIONSACHSE, INCL SONDERFALL OMEGA = 180
        !  ********************************************************
        A(1)=100.*(DRR(2,3,IZ)-DRR(3,2,IZ))/X
        A(2)=100.*(DRR(3,1,IZ)-DRR(1,3,IZ))/X
        A(3)=100.*(DRR(1,2,IZ)-DRR(2,1,IZ))/X
        GOTO 13
12      A(1)=SQRT((DRR(1,1,IZ)+1.)/2.)*100.
        A(2)=SQRT((DRR(2,2,IZ)+1.)/2.)*100.
        A(3)=SQRT((DRR(3,3,IZ)+1.)/2.)*100.
13      DO 11 I=1,3
           IF(A(I).EQ.0.)GOTO 19
           IA(I)=INT(A(I)+0.5*A(I)/ABS(A(I)))
           GOTO 11
19         IA(I)=0
11      END DO

        if(abs(alpha).lt.al_min) then
           al_min=abs(alpha)
           ac_min(1)=ia(1)
           ac_min(2)=ia(2)
           ac_min(3)=ia(3)
        end if

7    END DO
     IF(I31.EQ.1)P1(2)=180.+P1(2)
     IF(I31.EQ.2)P1(2)=360.-P1(2)
     IF(I31.EQ.2)P2(2)=90.-P2(2)
     IF(I31.EQ.3)P1(2)=180.+P1(2)


31 END DO
  i_min=100
  if(abs(ac_min(1)).lt.i_min.and.abs(ac_min(1)).gt.0) i_min=abs(ac_min(1))
  if(abs(ac_min(2)).lt.i_min.and.abs(ac_min(2)).gt.0) i_min=abs(ac_min(2))
  if(abs(ac_min(3)).lt.i_min.and.abs(ac_min(3)).gt.0) i_min=abs(ac_min(3))

  do 5501 ii=1,3
     ac_min(ii)=nint(float(ac_min(ii))/float(i_min))
5501 end do


  return
end subroutine misorinv

!  ______________________

subroutine qrvec(qq,din,dout)

  !  uses quaternion to rotate a vector from DIN to DOUT
  !  note the minus sign in the line with the permutation tensor
  !  we use this version to generate pole figures, which means that
  !  we are effectively performing the reverse transformation
  !  from crystal to sample axes (frames)

  !  cleaned up code  vi 08, ADR
  !  note that -q gives the same result as it should (represents same rotation)

  double precision qq(4),din(3),dout(3)
  double precision t1 , tmp

  !  CODE::
  !        print*,'qrvec: qq: ',qq
  t1 = qq(4)**2-qq(1)**2-qq(2)**2-qq(3)**2
  do i = 1,3
     dout(i) = t1*din(i)

     tmp = 0.
     do j = 1,3
        tmp = tmp + (qq(j)*din(j))
     enddo
     dout(i) = dout(i) + 2. * qq(i) * tmp

     do j = 1,3
        if(i.ne.j) then
           do 80, ijk = 1,3
              if(i.ne.ijk.and.j.ne.ijk) k=ijk
80         end do
           kl = j-i
           if( kl .eq.  2 ) kl = -1
           if( kl .eq. -2 ) kl =  1
           !  poor man's permutation tensor!
           !                 print*,'i,j,k,kl: ',i,j,k,kl
           dout(i) = dout(i) - (2.*qq(k)*qq(4)*float(kl)*din(j))
        endif   !  if(i.ne.j) then
     enddo

  enddo

  return
end subroutine qrvec

! _____________________________
!  }}}}
subroutine invqrvec( qq , din , dout )

  implicit none

  !  uses quaternion to rotate a vector from DIN to DOUT;
  !  but in the inverse sense so as to give sample to crystal axes
  !  derived from qrvec, above routine
  !  note that -q gives the same result as it should (represents same rotation)

  double precision :: qq(4),din(3),dout(3)
  double precision :: t1, tmp
  integer :: i , j , ijk, k , kl
  !  CODE::

  t1 = qq(4)**2-qq(1)**2-qq(2)**2-qq(3)**2
  do i = 1,3
     dout(i) = t1*din(i)

     tmp = 0.
     do j = 1,3
        tmp = tmp + (qq(j)*din(j))
     enddo
     dout(i) = dout(i) + 2. * qq(i) * tmp

     do j = 1,3
        if(i.ne.j) then
           do 80, ijk = 1,3
              if(i.ne.ijk.and.j.ne.ijk) k=ijk
80         end do
           kl = j-i
           if( kl .eq.  2 ) kl = -1
           if( kl .eq. -2 ) kl =  1
           !  poor man's permutation tensor!
           !                 print*,'i,j,k,kl: ',i,j,k,kl
           dout(i) = dout(i) + (2.*qq(k)*qq(4)*float(kl)*din(j))
           !  the ONLY difference between this and QRVE!  is the + instead of - above
        endif   !  if(i.ne.j) then
     enddo

  enddo
  return
end subroutine invqrvec

!  ----------------


subroutine axisang2rot(axis,ang,rot)

  implicit none

  !  AXIS  is the normalized axis vector
  !  ANG is the rotation angle in radians
  !  ROT is the 3x3 matrix
  double precision cc,ss,rnorm,rsum,axis(3),axis2(3),ang,rot(3,3)
  integer i

  !      logical verbose
  !      common /general/ verbose
  !  CODE::
  rsum=0.
  do 100, i=1,3
     axis2(i)=axis(i)**2
     rsum=rsum+axis2(i)
100 end do
  rnorm = sqrt(rsum)
  do 200, i=1,3
     axis(i)=axis(i)/rnorm
     axis2(i)=axis(i)**2
200 end do
  !  normalize anyway!
  cc = cos(ang)
  ss = sin(ang)
  rot(1,1) = axis2(1)+cc*(1-axis2(1))
  rot(2,2) = axis2(2)+cc*(1-axis2(2))
  rot(3,3) = axis2(3)+cc*(1-axis2(3))
  rot(1,2) = (1.-cc)*axis(1)*axis(2)+ss*axis(3)
  rot(2,1) = (1.-cc)*axis(1)*axis(2)-ss*axis(3)
  rot(2,3) = (1.-cc)*axis(2)*axis(3)+ss*axis(1)
  rot(3,2) = (1.-cc)*axis(2)*axis(3)-ss*axis(1)
  rot(3,1) = (1.-cc)*axis(1)*axis(3)+ss*axis(2)
  rot(1,3) = (1.-cc)*axis(1)*axis(3)-ss*axis(2)

  return
end subroutine axisang2rot

! ________________________________________

subroutine vsort( b , a , n )

  implicit none

  !     sorts components of a vector, 
  !  edited to separate input, b(3) & output, a(3)

  double precision :: a(4) , b(4) , rmax , rmin , rhigh , rlow
  integer :: n , i , imax , imin , ihigh , ilow
  !     write(*,*) 're-sorting components of vector'

  !  CODE::
  !      n = size ( b , 1 )
  if ( n < 3 ) stop 'illegal request, too small'
  if ( n > 4 ) stop 'illegal request, too large'
  a(1) = abs(b(1))
  a(2) = abs(b(2))
  a(3) = abs(b(3))
  if ( n == 4 ) a(4) = abs(b(4))
  rmin = a(1)
  rmax = a(1)
  imax = 1
  imin = 1
  do 10, i = 2 , n
     if(a(i).ge.rmax) then
        rmax = a(i)
        imax = i
     endif
     if(a(i).le.rmin) then
        rmin = a(i)
        imin = i
     endif
10 end do

  do 15, i = 1 , n
     if( i.ne.imax .and. i.ne.imin ) then
        rhigh = abs(a(i))
        ihigh = i
     end if
15 end do

  if ( n == 4 ) then
     do i = 1 , n
        if( i.ne.imax .and. i.ne.imin .and. i/=ihigh ) then
           rlow = abs(a(i))
           ilow = i
           if ( rhigh < rlow ) then  !  replace the values
              rlow = rhigh
              ilow = ihigh
              rhigh = abs(a(i))
              ihigh = i
           end if
        end if
     end do
  end if

  a(1) = rmax
  if ( n == 3 ) then
     a(2) = rhigh
     a(3) = rmin
  else
     a(2) = rhigh
     a(3) = rlow
     a(4) = rmin
  end if
  !     write(*,*) 'sorted a= ',a
  !     above code sorts a into descending order

  return
end subroutine vsort

!     
!     ________________________________________
subroutine qnorm(q)

  implicit none
  double precision :: q(4),norm
  integer :: i

  !     CODE::
  norm = sqrt( q(1)**2 + q(2)**2 + q(3)**2 + q(4)**2 )
  do i = 1,4
     q(i) = q(i) / norm
  end do
  return
end subroutine qnorm

! ________________________________________

!     include 'misc.f'

!..........................
subroutine load_cubic_crystal_symmetry_quat( verbose )

  use var

  implicit none

  logical :: verbose

!!$  real quatsymm(4,48) , quatsamp(4,4)
!!$  integer numsymm, numsamp
!!$  common/a1symm/ quatsymm , numsymm , quatsamp , numsamp
  integer :: i , j
  !data quatsymm(4,1) / 0. /

  !  CODE::
  numsymm = 24
  !      Number_Crystal_Symms = numsymm
  do i = 1 , 48
     quatsymm(1,i) = 0. 
     quatsymm(2,i) = 0. 
     quatsymm(3,i) = 0. 
     quatsymm(4,i) = 0. 
  end do
  quatsymm(4,1) = 1 
  quatsymm(1,2) = 1 
  quatsymm(2,3) = 1 
  quatsymm(3,4) = 1   !  180 about Z
  quatsymm(1,5) = 0.707107    !  90 about X
  quatsymm(4,5) = 0.707107 
  quatsymm(2,6) = 0.707107    !  90 about Y 
  quatsymm(4,6) = 0.707107 
  quatsymm(3,7) = 0.707107    !  90 about Z
  quatsymm(4,7) = 0.707107 
  quatsymm(1,8) = -0.707107 
  quatsymm(4,8) =  0.707107 
  quatsymm(2,9) = -0.707107 
  quatsymm(4,9) =  0.707107 
  quatsymm(3,10) = -0.707107 
  quatsymm(4,10) =  0.707107 
  quatsymm(1,11) =  0.707107 
  quatsymm(2,11) =  0.707107 
  quatsymm(1,12) = -0.707107 
  quatsymm(2,12) =  0.707107 
  quatsymm(2,13) =  0.707107 
  quatsymm(3,13) =  0.707107 
  quatsymm(2,14) = -0.707107 
  quatsymm(3,14) =  0.707107 
  quatsymm(1,15) =  0.707107 
  quatsymm(3,15) =  0.707107 
  quatsymm(1,16) = -0.707107 
  quatsymm(3,16) =  0.707107 
  quatsymm(1,17) =  0.5      !  120 about 111
  quatsymm(2,17) =  0.5      !  120 about 111
  quatsymm(3,17) =  0.5      !  120 about 111
  quatsymm(4,17) =  0.5      !  120 about 111
  quatsymm(1,18) = -0.5 
  quatsymm(2,18) = -0.5 
  quatsymm(3,18) = -0.5 
  quatsymm(4,18) =  0.5 
  quatsymm(1,19) =  0.5 
  quatsymm(2,19) = -0.5 
  quatsymm(3,19) =  0.5 
  quatsymm(4,19) =  0.5 
  quatsymm(1,20) = -0.5 
  quatsymm(2,20) =  0.5 
  quatsymm(3,20) = -0.5 
  quatsymm(4,20) =  0.5 
  quatsymm(1,21) = -0.5 
  quatsymm(2,21) =  0.5 
  quatsymm(3,21) =  0.5 
  quatsymm(4,21) =  0.5 
  quatsymm(1,22) =  0.5 
  quatsymm(2,22) = -0.5 
  quatsymm(3,22) = -0.5 
  quatsymm(4,22) =  0.5 
  quatsymm(1,23) = -0.5 
  quatsymm(2,23) = -0.5 
  quatsymm(3,23) =  0.5 
  quatsymm(4,23) =  0.5 
  quatsymm(1,24) =  0.5 
  quatsymm(2,24) =  0.5 
  quatsymm(3,24) = -0.5 
  quatsymm(4,24) =  0.5

  if (verbose) then
     print * , 'Quats for cubic symmetry operators'
     print* , 'No. quatsymm(1-4)'
     do i = 1 , numsymm
        print* , i , (quatsymm( j , i ) , j = 1, 4)
     end do
  end if

  return
end subroutine load_cubic_crystal_symmetry_quat

!// .......................................................
subroutine load_hcp_crystal_symmetry_quat( verbose )

  use var

  implicit none

  logical verbose

  !  /*  from quat.symm.hex
  !   0    0    0    1
  !   0   0    0.5   0.8660254
  !   0   0   0.8660254   0.5
  !   0    0    1    0
  !   0   0   0.8660254  -0.5
  !   0   0    0.5  -0.8660254
  !   1    0    0    0
  !   0.8660254   0.5    0    0
  !   0.5   0.8660254    0    0
  !   0    1    0    0
  !  -0.5   0.8660254    0    0
  !  -0.8660254   0.5    0    0
  !   */

!!$  real quatsymm(4,48) , quatsamp(4,4)
!!$  integer numsymm, numsamp
!!$  common/a1symm/ quatsymm , numsymm , quatsamp , numsamp
  integer :: i , j

  !  CODE::
  numsymm = 12
  !      Number_Crystal_Symms = numsymm
  do i = 1 , numsymm
     quatsymm(1,i) = 0.
     quatsymm(2,i) = 0.
     quatsymm(3,i) = 0.
     quatsymm(4,i) = 0.
  end do
  quatsymm(4,1) = 1
  quatsymm(3,2) =  0.5
  quatsymm(4,2) =  0.8660254
  quatsymm(3,3) =  0.8660254
  quatsymm(4,3) =  0.5
  quatsymm(3,4) = 1
  quatsymm(3,5) =  0.8660254
  quatsymm(4,5) =  -0.5
  quatsymm(3,6) =  0.5
  quatsymm(4,6) =  -0.8660254
  quatsymm(1,7) = 1
  quatsymm(1,8) =  0.8660254
  quatsymm(2,8) =  0.5
  quatsymm(1,9) =  0.5
  quatsymm(2,9) =  0.8660254
  quatsymm(2,10) = 1
  quatsymm(1,11) =  -0.5
  quatsymm(2,11) =  0.8660254
  quatsymm(1,12) =  -0.8660254
  quatsymm(2,12) =  0.5


  if (verbose) then
     print * , 'Quats for cubic symmetry operators'
     print* , 'No. quatsymm(1-4)'
     do i = 1 , numsymm
        print* , i , (quatsymm( j , i ) , j = 1, 4)
     end do
  end if

  return
end subroutine load_hcp_crystal_symmetry_quat

!..........................
subroutine load_orthorhombic_symmetry_quat( verbose )

  use var

  implicit none

  logical :: verbose

  integer :: i , j

  !  CODE::
  numsymm = 4
  do i = 1 , numsymm
     quatsymm(1,i) = 0. 
     quatsymm(2,i) = 0. 
     quatsymm(3,i) = 0. 
     quatsymm(4,i) = 0. 
  end do
  quatsymm(4,1) = 1   !  identity
  quatsymm(1,2) = 1   !  180 about X
  quatsymm(2,3) = 1   !  180 about Y
  quatsymm(3,4) = 1   !  180 about Z

  if (verbose) then
     print * , 'Quats for ortho symmetry operators'
     print* , 'No. quatsymm(1-4)'
     do i = 1 , numsymm
        print* , i , (quatsymm( j , i ) , j = 1, 4)
     end do
  end if

  return
end subroutine load_orthorhombic_symmetry_quat

!..........................
subroutine load_orthorhombic_sample_symmetry_quat( verbose )

  use var

  implicit none

  logical :: verbose

  integer :: i , j

  !  CODE::
  numsamp = 4
  do i = 1 , numsamp
     quatsamp(1,i) = 0. 
     quatsamp(2,i) = 0. 
     quatsamp(3,i) = 0. 
     quatsamp(4,i) = 0. 
  end do
  quatsamp(4,1) = 1   !  identity
  quatsamp(1,2) = 1   !  180 about X
  quatsamp(2,3) = 1   !  180 about Y
  quatsamp(3,4) = 1   !  180 about Z

  if (verbose) then
     print * , 'Quats for ortho SAMPLE symmetry operators'
     print* , 'No. quatsymm(1-4)'
     do i = 1 , numsamp
        print* , i , (quatsamp( j , i ) , j = 1, 4)
     end do
  end if

  return
end subroutine load_orthorhombic_sample_symmetry_quat


!..........................
subroutine load_CSL_quat( verbose )

  use var

  implicit none

  integer :: i
  logical :: verbose

!quatcsl(4,numsigmas)
!Sigma	q1	q2	q3	q4
  numcsl = 28
  i = 0
!     read(3,*) nsig(i),theta,lmn,r1,r2,r3,(quatcsl(ijk,i),ijk=1,4)
  i = i + 1
  nsig(i) = 3
  i = i + 1
  nsig(i) = 5
  i = i + 1
  nsig(i) = 7
  i = i + 1
  nsig(i) = 9
  i = i + 1
  nsig(i) = 11
  i = i + 1
  nsig(i) = 13
  i = i + 1
  nsig(i) = 13
  i = i + 1
  nsig(i) = 15
  i = i + 1
  nsig(i) = 17
  i = i + 1
  nsig(i) = 17
  i = i + 1
  nsig(i) = 19
  i = i + 1
  nsig(i) = 19
  i = i + 1
  nsig(i) = 21
  i = i + 1
  nsig(i) = 21
  i = i + 1
  nsig(i) = 23
  i = i + 1
  nsig(i) = 25
  i = i + 1
  nsig(i) = 25
  i = i + 1
  nsig(i) = 27
  i = i + 1
  nsig(i) = 27
  i = i + 1
  nsig(i) = 29
  i = i + 1
  nsig(i) = 29
  i = i + 1
  nsig(i) = 31
  i = i + 1
  nsig(i) = 31
  i = i + 1
  nsig(i) = 33
  i = i + 1
  nsig(i) = 33
  i = i + 1
  nsig(i) = 33
  i = i + 1
  nsig(i) = 35
  i = i + 1
  nsig(i) = 35
!  now for the quat values
  i = 0
  i = i + 1
  quatcsl( 1 , i ) = 0.288675135
  i = i + 1
  quatcsl( 1 , i ) = 0.000
  i = i + 1
  quatcsl( 1 , i ) = 0.188
  i = i + 1
  quatcsl( 1 , i ) = 0.000
  i = i + 1
  quatcsl( 1 , i ) = 0.000
  i = i + 1
  quatcsl( 1 , i ) = 0.000
  i = i + 1
  quatcsl( 1 , i ) = 0.139
  i = i + 1
  quatcsl( 1 , i ) = 0.000
  i = i + 1
  quatcsl( 1 , i ) = 0.000
  i = i + 1
  quatcsl( 1 , i ) = 0.171
  i = i + 1
  quatcsl( 1 , i ) = 0.000
  i = i + 1
  quatcsl( 1 , i ) = 0.229
  i = i + 1
  quatcsl( 1 , i ) = 0.109
  i = i + 1
  quatcsl( 1 , i ) = 0.154
  i = i + 1
  quatcsl( 1 , i ) = 0.104
  i = i + 1
  quatcsl( 1 , i ) = 0.000
  i = i + 1
  quatcsl( 1 , i ) = 0.100
  i = i + 1
  quatcsl( 1 , i ) = 0.000
  i = i + 1
  quatcsl( 1 , i ) = 0.000
  i = i + 1
  quatcsl( 1 , i ) = 0.000
  i = i + 1
  quatcsl( 1 , i ) = 0.131
  i = i + 1
  quatcsl( 1 , i ) = 0.09
  i = i + 1
  quatcsl( 1 , i ) = 0.180
  i = i + 1
  quatcsl( 1 , i ) = 0.000
  i = i + 1
  quatcsl( 1 , i ) = 0.087
  i = i + 1
  quatcsl( 1 , i ) = 0.000
  i = i + 1
  quatcsl( 1 , i ) = 0.119
  i = i + 1
  quatcsl( 1 , i ) = 0.083
  i = 0
  i = i + 1
  quatcsl( 2 , i ) = 0.288675135
  i = i + 1
  quatcsl( 2 , i ) = 0.000
  i = i + 1
  quatcsl( 2 , i ) = 0.188
  i = i + 1
  quatcsl( 2 , i ) = 0.236  !  S9, corrected
  i = i + 1
  quatcsl( 2 , i ) = 0.302
  i = i + 1
  quatcsl( 2 , i ) = 0.000
  i = i + 1
  quatcsl( 2 , i ) = 0.139
  i = i + 1
  quatcsl( 2 , i ) = 0.183
  i = i + 1
  quatcsl( 2 , i ) = 0.000
  i = i + 1
  quatcsl( 2 , i ) = 0.343
  i = i + 1
  quatcsl( 2 , i ) = 0.162
  i = i + 1
  quatcsl( 2 , i ) = 0.229
  i = i + 1
  quatcsl( 2 , i ) = 0.109
  i = i + 1
  quatcsl( 2 , i ) = 0.154
  i = i + 1
  quatcsl( 2 , i ) = 0.104
  i = i + 1
  quatcsl( 2 , i ) = 0.000
  i = i + 1
  quatcsl( 2 , i ) = 0.300
  i = i + 1
  quatcsl( 2 , i ) = 0.193
  i = i + 1
  quatcsl( 2 , i ) = 0.136
  i = i + 1
  quatcsl( 2 , i ) = 0.000
  i = i + 1
  quatcsl( 2 , i ) = 0.263
  i = i + 1
  quatcsl( 2 , i ) = 0.09
  i = i + 1
  quatcsl( 2 , i ) = 0.18
  i = i + 1
  quatcsl( 2 , i ) = 0.123
  i = i + 1
  quatcsl( 2 , i ) = 0.087
  i = i + 1
  quatcsl( 2 , i ) = 0.348
  i = i + 1
  quatcsl( 2 , i ) = 0.119
  i = i + 1
  quatcsl( 2 , i ) = 0.253
  i = 0
  i = i + 1
  quatcsl( 3 , i ) = 0.288675135
  i = i + 1
  quatcsl( 3 , i ) = 0.316
  i = i + 1
  quatcsl( 3 , i ) = 0.188
  i = i + 1
  quatcsl( 3 , i ) = 0.236  !  S9, corrected
  i = i + 1
  quatcsl( 3 , i ) = 0.302
  i = i + 1
  quatcsl( 3 , i ) = 0.196
  i = i + 1
  quatcsl( 3 , i ) = 0.139
  i = i + 1
  quatcsl( 3 , i ) = 0.365
  i = i + 1
  quatcsl( 3 , i ) = 0.243
  i = i + 1
  quatcsl( 3 , i ) = 0.343
  i = i + 1
  quatcsl( 3 , i ) = 0.162
  i = i + 1
  quatcsl( 3 , i ) = 0.229
  i = i + 1
  quatcsl( 3 , i ) = 0.109
  i = i + 1
  quatcsl( 3 , i ) = 0.308
  i = i + 1
  quatcsl( 3 , i ) = 0.313
  i = i + 1
  quatcsl( 3 , i ) = 0.142
  i = i + 1
  quatcsl( 3 , i ) = 0.300
  i = i + 1
  quatcsl( 3 , i ) = 0.193
  i = i + 1
  quatcsl( 3 , i ) = 0.272
  i = i + 1
  quatcsl( 3 , i ) = 0.393
  i = i + 1
  quatcsl( 3 , i ) = 0.263
  i = i + 1
  quatcsl( 3 , i ) = 0.09
  i = i + 1
  quatcsl( 3 , i ) = 0.359
  i = i + 1
  quatcsl( 3 , i ) = 0.123
  i = i + 1
  quatcsl( 3 , i ) = 0.261
  i = i + 1
  quatcsl( 3 , i ) = 0.348
  i = i + 1
  quatcsl( 3 , i ) = 0.239
  i = i + 1
  quatcsl( 3 , i ) = 0.253
  i = i + 1
  i = 0
  i = i + 1
  quatcsl( 4 , i ) = 0.866025404
  i = i + 1
  quatcsl( 4 , i ) = 0.948
  i = i + 1
  quatcsl( 4 , i ) = 0.944
  i = i + 1
  quatcsl( 4 , i ) = 0.943
  i = i + 1
  quatcsl( 4 , i ) = 0.904
  i = i + 1
  quatcsl( 4 , i ) = 0.981
  i = i + 1
  quatcsl( 4 , i ) = 0.971
  i = i + 1
  quatcsl( 4 , i ) = 0.913
  i = i + 1
  quatcsl( 4 , i ) = 0.970
  i = i + 1
  quatcsl( 4 , i ) = 0.858
  i = i + 1
  quatcsl( 4 , i ) = 0.973
  i = i + 1
  quatcsl( 4 , i ) = 0.918
  i = i + 1
  quatcsl( 4 , i ) = 0.982
  i = i + 1
  quatcsl( 4 , i ) = 0.926
  i = i + 1
  quatcsl( 4 , i ) = 0.938
  i = i + 1
  quatcsl( 4 , i ) = 0.99
  i = i + 1
  quatcsl( 4 , i ) = 0.9
  i = i + 1
  quatcsl( 4 , i ) = 0.962
  i = i + 1
  quatcsl( 4 , i ) = 0.953
  i = i + 1
  quatcsl( 4 , i ) = 0.919
  i = i + 1
  quatcsl( 4 , i ) = 0.919
  i = i + 1
  quatcsl( 4 , i ) = 0.988
  i = i + 1
  quatcsl( 4 , i ) = 0.898
  i = i + 1
  quatcsl( 4 , i ) = 0.985
  i = i + 1
  quatcsl( 4 , i ) = 0.957
  i = i + 1
  quatcsl( 4 , i ) = 0.870
  i = i + 1
  quatcsl( 4 , i ) = 0.956
  i = i + 1
  quatcsl( 4 , i ) = 0.93

return

end subroutine load_CSL_quat

