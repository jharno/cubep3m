!! C  (C) Copr. 1986-92 Numerical Recipes Software .-35?421.1-9.
!! Hack of sort2 to be an indexed sort
!! Change brr,b from real to integers and use integer itemp

      SUBROUTINE indexedsort(n,arr,brr)
      IMPLICIT NONE
      INTEGER n,M,NSTACK
      REAL arr(n)
      INTEGER brr(n)
      PARAMETER (M=7,NSTACK=50)
      INTEGER i,ir,j,jstack,k,l,istack(NSTACK)
      INTEGER b,itemp
      REAL a,temp
      jstack=0
      l=1
      ir=n
      do i=1,n
         brr(i)=i
      enddo

1     if(ir-l.lt.M)then
        do 12 j=l+1,ir
          a=arr(j)
          b=brr(j)
          do 11 i=j-1,1,-1
            if(arr(i).le.a)goto 2
            arr(i+1)=arr(i)
            brr(i+1)=brr(i)
11        continue
          i=0
2         arr(i+1)=a
          brr(i+1)=b
12      continue
        if(jstack.eq.0)return
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
        k=(l+ir)/2
        temp=arr(k)
        arr(k)=arr(l+1)
        arr(l+1)=temp
        itemp=brr(k)
        brr(k)=brr(l+1)
        brr(l+1)=itemp
        if(arr(l+1).gt.arr(ir))then
          temp=arr(l+1)
          arr(l+1)=arr(ir)
          arr(ir)=temp
          itemp=brr(l+1)
          brr(l+1)=brr(ir)
          brr(ir)=itemp
        endif
        if(arr(l).gt.arr(ir))then
          temp=arr(l)
          arr(l)=arr(ir)
          arr(ir)=temp
          itemp=brr(l)
          brr(l)=brr(ir)
          brr(ir)=itemp
        endif
        if(arr(l+1).gt.arr(l))then
          temp=arr(l+1)
          arr(l+1)=arr(l)
          arr(l)=temp
          itemp=brr(l+1)
          brr(l+1)=brr(l)
          brr(l)=itemp
        endif
        i=l+1
        j=ir
        a=arr(l)
        b=brr(l)
3       continue
          i=i+1
        if(arr(i).lt.a)goto 3
4       continue
          j=j-1
        if(arr(j).gt.a)goto 4
        if(j.lt.i)goto 5
        temp=arr(i)
        arr(i)=arr(j)
        arr(j)=temp
        itemp=brr(i)
        brr(i)=brr(j)
        brr(j)=itemp
        goto 3
5       arr(l)=arr(j)
        arr(j)=a
        brr(l)=brr(j)
        brr(j)=b
        jstack=jstack+2
        if(jstack.gt.NSTACK)pause 'NSTACK too small in sort2'
        if(ir-i+1.ge.j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        endif
      endif
      goto 1
      END

