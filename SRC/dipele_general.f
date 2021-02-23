       subroutine dipele(r1,r2,cgam,dx,dy,dz,nelec,nelecmax)
       implicit real*8(a-h,o-z)

       ! r1 is the distance 01 (small r in Jacobi coordinates)
       ! r2 is the distance between the c-o-m of 01 to atom 2 (big R in
       !                                        Jacobi coordinates
       ! cgam is cos(gamma), gamma being the angle between r1 and r2
       !                      vectors
       ! distances and dipole are in atomic units

       dimension dx(nelecmax),dy(nelecmax),dz(nelecmax)

        return
        end

        subroutine setdipini
        implicit real*8(a-h,o-z)

        return
        end
