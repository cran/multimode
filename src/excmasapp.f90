
Subroutine calcdist(dataor,ndata,modtemp,maxiter,posx,lenposx,distCext)

Integer ndata,modtemp,maxiter,i,j,k,l,combi(2*modtemp),weightdata,lenposx
Double Precision dataor(ndata),posx(lenposx),distdata,distCext(ndata)

do i=1,2*modtemp
combi(i)=i
end do
		
do i=1,maxiter-1
						
distdata=-SUM((/ (posx(combi(2*l-1)),l=1,modtemp) /)) + SUM((/ (posx(combi(2*l)),l=1,modtemp) /)) 
			
weightdata=0
do k=1,modtemp
weightdata=weightdata-COUNT((/ (dataor(l)<posx(combi(2*k-1)),l=1,ndata) /))
weightdata=weightdata+COUNT((/ (dataor(l)<=posx(combi(2*k)),l=1,ndata) /))
end do
			
if(weightdata>0) then
distCext(weightdata)=distdata
end if
			
combi(2*modtemp)=combi(2*modtemp)+1
			

if (combi(2*modtemp)==(lenposx+1)) then

do j=2*modtemp,2,-1
				
if (combi(j)>(lenposx-(2*modtemp-j))) then
combi(j-1)=combi(j-1)+1
do k=j,2*modtemp
combi(k)=combi(k-1)+1
end do
end if

end do
end if

			
end do

i=maxiter
								
distdata=-SUM((/ (posx(combi(2*l-1)),l=1,modtemp) /)) + SUM((/ (posx(combi(2*l)),l=1,modtemp) /)) 
			
weightdata=0
do k=1,modtemp
weightdata=weightdata-COUNT((/ (dataor(l)<posx(combi(2*k-1)),l=1,ndata) /))
weightdata=weightdata+COUNT((/ (dataor(l)<=posx(combi(2*k)),l=1,ndata) /))
end do
			
if(weightdata>0) then
distCext(weightdata)=distdata
end if
			
		
		
Return
End