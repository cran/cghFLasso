       subroutine avesmooth (n, p, size, matri, result)
       implicit double precision (a-h,o-z)
       integer n,p, size
       double precision matri(n,p), result(n,p)

       integer i,j,k, t, count
       double precision tempsum
       k=(size-1)/2 
       do 100 j=1, p
          do 200 i=1,n
            if(matri(i,j) .EQ. 999) then
             tempsum=0
             count=0
             do 300 t=i-k, i+k
              if((t .GT. 0) .AND. (t .LE. n)) then
                    if(matri(t,j) .NE. 999) then
                        tempsum=tempsum+matri(t,j)
                        count=count+1
                    endif
               endif
 300           continue
             if(count .GT. 0) then
                   result(i,j)=tempsum/count
             else
                   result(i,j)=999
             endif
            endif
 200       continue
 100       continue

       return
       end
