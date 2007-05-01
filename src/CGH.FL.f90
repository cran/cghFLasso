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

      subroutine dflas (al1,n,y,mxl2,dlm,thr,nl2,a,al2,nlp,jerr)            
      parameter(eps=1.0e-6)                                                 
      real y(n),a(n,mxl2),al2(mxl2)                                         
      real, dimension (:), allocatable :: b,w,z                              
      integer, dimension (:), allocatable :: m                               
      allocate(b(1:n),stat=jerr)                                            
      allocate(w(1:n),stat=ierr)                                            
      jerr=jerr+ierr                                                        
      allocate(z(1:n),stat=ierr)                                            
      jerr=jerr+ierr                                                        
      allocate(m(1:n),stat=ierr)                                            
      jerr=jerr+ierr                                                        
      if(jerr.ne.0) return                                                  
      sthr=sum(y)/n                                                         
      sthr=sum((y-sthr)**2)/n                                               
      rthr=eps*sthr                                                         
      sthr=thr*sqrt(sthr)                                                   
      nlp=0                                                                 
      z=y                                                                   
      b=sign(max(0.0,abs(z)-al1),z)                                         
      nl2=1                                                                 
      al2(1)=0.0                                                            
      a(:,1)=b                                                              
      if(maxval(abs(b)).le.0.0) return                                      
10010 do 10011 i=1,n-1                                                      
      m(i)=i+1                                                              
10011 continue                                                              
10012 continue                                                              
      m(n)=0                                                                
      w=1.0                                                                 
10020 do 10021 i=2,mxl2                                                     
      al2(i)=al2(i-1)+dlm                                                   
      call dfl(al1,al2(i),z,w,m,rthr,sthr,b,mlp)                            
      call fuse(z,w,b,m)                                                    
      nlp=nlp+mlp                                                           
      nl2=nl2+1                                                             
      call coeffs(m,b,n,a(:,nl2))                                           
      if(m(1).eq.0)goto 10022                                               
      if(m(m(1)).eq.0)goto 10022                                            
10021 continue                                                              
10022 continue                                                              
      return                                                                
      end                                                                   
      subroutine dfl (al1,al2,y,w,m,rthr,thr,a,nlp)                         
      real y(*),w(*),a(*)                                                   
      integer m(*)                                                          
      nlp=0                                                                 
10030 continue                                                              
10031 continue                                                              
      nlp=nlp+1                                                             
      dlx=0.0                                                               
      call dfl2(1,y,w,a,m,al1,al2,rthr,v)                                   
      dlx=max(abs(v-a(1)),dlx)                                              
      a(1)=v                                                                
      i=1                                                                   
10040 continue                                                              
10041 continue                                                              
      i=m(i)                                                                
      call dfl2(i,y,w,a,m,al1,al2,rthr,v)                                   
      dlx=max(abs(v-a(i)),dlx)                                              
      a(i)=v                                                                
      if(m(i).eq.0)goto 10042                                               
      goto 10041                                                            
10042 continue                                                              
      if(dlx.lt.thr)goto 10032                                              
      goto 10031                                                            
10032 continue                                                              
      return                                                                
      end                                                                   
      subroutine dfl2(i,y,w,a,m,al1,al2,eps,ai)                             
      real y(*),w(*),a(*)                                                   
      integer m(*)                                                          
      if(i .ne. 1)goto 10061                                                
      v=edge(y(1),a(m(1)),al1,al2/w(1),al1+al2/w(1))                        
      goto 10051                                                            
10061 if(m(i) .le. 0)goto 10071                                             
      v=soln(y(i),ai,a(m(i)),al1,al2/w(i),al1+2.0*al2/w(i))                 
      goto 10081                                                            
10071 continue                                                              
      ai=edge(y(i),ai,al1,al2/w(i),al1+al2/w(i))                            
      return                                                                
10081 continue                                                              
10051 continue                                                              
      if(v .eq. a(i))goto 10101                                             
      ai=v                                                                  
      return                                                                
10101 continue                                                              
      j=m(i)                                                                
      k=m(j)                                                                
      ww=w(i)+w(j)                                                          
      yy=(w(i)*y(i)+w(j)*y(j))/ww                                           
      if(i .ne. 1)goto 10121                                                
      if(k .le. 0)goto 10141                                                
      vv=edge(yy,a(k),al1,al2/ww,al1+al2/ww)                                
      goto 10151                                                            
10141 continue                                                              
      vv=yy                                                                 
10151 continue                                                              
10131 continue                                                              
      di=dcri(1,vv,a(1),a(j),y(1),y(j),w(1),w(j),dum,a(k),al1,al2)          
      goto 10111                                                            
10121 if(k .eq. 0)goto 10161                                                
      vv=soln(yy,ai,a(k),al1,al2/ww,al1+2.0*al2/ww)                         
      di=dcri(2,vv,a(i),a(j),y(i),y(j),w(i),w(j),ai,a(k),al1,al2)           
      goto 10171                                                            
10161 continue                                                              
      vv=edge(yy,ai,al1,al2/ww,al1+al2/ww)                                  
      di=dcri(3,vv,a(i),a(j),y(i),y(j),w(i),w(j),ai,dum,al1,al2)            
10171 continue                                                              
10111 continue                                                              
      if(di .gt. eps)goto 10191                                             
      ai=a(i)                                                               
      return                                                                
10191 continue                                                              
      ai=vv                                                                 
      a(j)=ai                                                               
      return                                                                
      end                                                                   
      function dcri (k,af,a,ap,y,yp,w,wp,am,app,al1,al2)                    
      dcri=0.5*(w*((y-a)**2-(y-af)**2)+wp*((yp-ap)**2-(yp-af)**2))  +al1*(w*abs(a)+wp*abs(ap)-(w+wp)*abs(af))+al2*abs(a-ap)
      if(k.ne.1) dcri=dcri+al2*(abs(a-am)-abs(af-am))                       
      if(k.ne.3) dcri=dcri+al2*(abs(ap-app)-abs(af-app))                    
      return                                                                
      end                                                                   
      subroutine fuse (y,w,a,m)                                             
      real y(*),w(*),a(*)                                                   
      integer m(*)                                                          
      i=1                                                                   
10200 continue                                                              
10201 if(i.eq.0)goto 10202                                                  
      if(a(i) .ne. 0.0)goto 10221                                           
      i=m(i)                                                                
      goto 10201                                                            
10221 continue                                                              
      j=m(i)                                                                
      if(j.eq.0) return                                                     
10230 continue                                                              
10231 if(a(i).ne.a(j))goto 10232                                            
      y(i)=w(i)*y(i)+w(j)*y(j)                                              
      w(i)=w(i)+w(j)                                                        
      y(i)=y(i)/w(i)                                                        
      j=m(j)                                                                
      if(j.eq.0)goto 10232                                                  
      goto 10231                                                            
10232 continue                                                              
      m(i)=j                                                                
      i=j                                                                   
      goto 10201                                                            
10202 continue                                                              
      return                                                                
      end                                                                   
      subroutine coeffs(m,a,n,ao)                                           
      integer m(n)                                                          
      real a(n),ao(n)                                                       
      i=1                                                                   
10240 continue                                                              
10241 continue                                                              
      j=m(i)                                                                
      if(j .ne. 0)goto 10261                                                
      ao(i:n)=a(i)                                                          
      goto 10242                                                            
10261 continue                                                              
      ao(i:(j-1))=a(i)                                                      
      i=j                                                                   
      goto 10241                                                            
10242 continue                                                              
      return                                                                
      end                                                                   
      function edge(y,a,btl,ombtl,alm)                                      
      if(a .ne. 0.0)goto 10281                                              
      if(y .ne. 0.0)goto 10301                                              
      edge=0.0                                                              
      return                                                                
10301 continue                                                              
      edge=sign(max(0.0,abs(y)-alm),y)                                      
      return                                                                
10281 continue                                                              
      z1=min(0.0,a)                                                         
      z2=max(0.0,a)                                                         
      x=z1-0.1                                                              
      edge=y-sign(btl,x)-sign(ombtl,x-a)                                    
      if(edge.le.z1) return                                                 
      x=0.5*(z1+z2)                                                         
      edge=y-sign(btl,x)-sign(ombtl,x-a)                                    
      if(edge.gt.z1.and.edge.le.z2) return                                  
      x=z2+0.1                                                              
      edge=y-sign(btl,x)-sign(ombtl,x-a)                                    
      if(edge.gt.z2) return                                                 
      x=z1                                                                  
      c1=0.5*(y-x)**2+btl*abs(x)+ombtl*abs(x-a)                             
      x=z2                                                                  
      c2=0.5*(y-x)**2+btl*abs(x)+ombtl*abs(x-a)                             
      if(c1 .ge. c2)goto 10321                                              
      edge=z1                                                               
      goto 10331                                                            
10321 continue                                                              
      edge=z2                                                               
10331 continue                                                              
10311 continue                                                              
      return                                                                
      end                                                                   
      function soln(y,al,au,btl,ombtl,tmbtl)                                
      if(al .ne. au)goto 10351                                              
      if(al .ne. 0.0)goto 10371                                             
      if(y .ne. 0.0)goto 10391                                              
      soln=0.0                                                              
      return                                                                
10391 continue                                                              
      soln=sign(max(0.0,abs(y)-tmbtl),y)                                    
      return                                                                
10371 continue                                                              
      z1=min(0.0,al)                                                        
      z2=max(0.0,al)                                                        
      x=z1-0.1                                                              
      soln=y-sign(btl,x)-sign(2.0*ombtl,x-al)                               
      if(soln.le.z1) return                                                 
      x=0.5*(z1+z2)                                                         
      soln=y-sign(btl,x)-sign(2.0*ombtl,x-al)                               
      if(soln.gt.z1.and.soln.le.z2) return                                  
      x=z2+0.1                                                              
      soln=y-sign(btl,x)-sign(2.0*ombtl,x-al)                               
      if(soln.gt.z2) return                                                 
      x=z1                                                                  
      c1=0.5*(y-x)**2+btl*abs(x)+2.0*ombtl*abs(x-al)                        
      x=z2                                                                  
      c2=0.5*(y-x)**2+btl*abs(x)+2.0*ombtl*abs(x-al)                        
      if(c1 .ge. c2)goto 10411                                              
      soln=z1                                                               
      goto 10421                                                            
10411 continue                                                              
      soln=z2                                                               
10421 continue                                                              
10401 continue                                                              
      return                                                                
10351 continue                                                              
      x1=0.0                                                                
      x2=al                                                                 
      x3=au                                                                 
      if(x2 .ge. x3)goto 10441                                              
      if(x1 .ge. x2)goto 10461                                              
      z1=x1                                                                 
      z2=x2                                                                 
      z3=x3                                                                 
      goto 10451                                                            
10461 if(x1 .ge. x3)goto 10471                                              
      z1=x2                                                                 
      z2=x1                                                                 
      z3=x3                                                                 
      goto 10481                                                            
10471 continue                                                              
      z1=x2                                                                 
      z2=x3                                                                 
      z3=x1                                                                 
10481 continue                                                              
10451 continue                                                              
      goto 10431                                                            
10441 if(x1 .ge. x3)goto 10491                                              
      z1=x1                                                                 
      z2=x3                                                                 
      z3=x2                                                                 
      goto 10431                                                            
10491 if(x2 .ge. x1)goto 10501                                              
      z1=x3                                                                 
      z2=x2                                                                 
      z3=x1                                                                 
      goto 10511                                                            
10501 continue                                                              
      z1=x3                                                                 
      z2=x1                                                                 
      z3=x2                                                                 
10511 continue                                                              
10431 continue                                                              
      x=z1-0.1                                                              
      soln=y-sign(btl,x)+sign(ombtl,au-x)-sign(ombtl,x-al)                  
      if(soln.le.z1) return                                                 
      x=0.5*(z1+z2)                                                         
      soln=y-sign(btl,x)+sign(ombtl,au-x)-sign(ombtl,x-al)                  
      if(soln.gt.z1.and.soln.le.z2) return                                  
      x=0.5*(z2+z3)                                                         
      soln=y-sign(btl,x)+sign(ombtl,au-x)-sign(ombtl,x-al)                  
      if(soln.gt.z2.and.soln.le.z3) return                                  
      x=z3+0.1                                                              
      soln=y-sign(btl,x)+sign(ombtl,au-x)-sign(ombtl,x-al)                  
      if(soln.gt.z3) return                                                 
      x=z1                                                                  
      c1=0.5*(y-x)**2+btl*abs(x)+ombtl*(abs(x-al)+abs(au-x))                
      x=z2                                                                  
      c2=0.5*(y-x)**2+btl*abs(x)+ombtl*(abs(x-al)+abs(au-x))                
      x=z3                                                                  
      c3=0.5*(y-x)**2+btl*abs(x)+ombtl*(abs(x-al)+abs(au-x))                
      if(c1 .ge. c2)goto 10531                                              
      if(c1 .ge. c3)goto 10551                                              
      soln=z1                                                               
      goto 10561                                                            
10551 continue                                                              
      soln=z3                                                               
10561 continue                                                              
10541 continue                                                              
      goto 10521                                                            
10531 if(c3 .ge. c2)goto 10571                                              
      soln=z3                                                               
      goto 10581                                                            
10571 continue                                                              
      soln=z2                                                               
10581 continue                                                              
10521 continue                                                              
      return                                                                
      end                                                                   
