# FluidMechanics
A program written for TI-Nspire CX CAS for calculation of nozzles
- Shock Waves
- Fanno process
- Rayleigh Process
``` 
Define LibPub fluidmech()=
Prgm
:
:Local r,k,p0,t0,d0,m,massflowrate,machpoints,area
:Local p0t,p0aux,poypox,ds,datafixer
:Local t,p,d,lambda,arearatio
:Local acr,areacritical,tcritical,pcritical,dcritical
:Local i,j
:i:=1
:j:=1
:
:r:=287.1
:k:=1.4
:p0:=500000
:t0:=300
:machpoints:=[[1,2.2,1.57][0,2,1][0,1.57,0]]
:massflowrate:=1
:Request "kappa",k
:Request "gas Constant",r
:Request "P0",p0
:Request "T0",t0
:Disp "d0",((p0)/(r*t0))
:Request "M",machpoints
:Request "massflowrate",massflowrate
:
:m:=newMat(2,colDim(machpoints))
:poypox:=m
:p0t:=m
:ds:=m
:p0aux:=p0
:j:=1
:For i,1,colDim(machpoints),1
:m[1,i]:=machpoints[1,i]
:p0t[1,i]:=p0aux
:
:If machpoints[2,i]>0 Then
:
:If machpoints[2,i]=1 Then
: m[2,i]:=√(((1+((k-1)/(2))*machpoints[1,i]^(2))/(k*machpoints[1,i]^(2)-((k-1)/(2)))))
:  poypox[1,i]:=((((k+1)*machpoints[1,i]^(2))/(2+(k-1)*machpoints[1,i]^(2))))^(((k)/(k-1)))*(((k+1)/(2*k*machpoints[1,i]^(2)-(k-1))))^(((1)/(k-1)))
:EndIf
:
:If machpoints[2,i]=2 Then
:poypox[1,i]:=((machpoints[1,i])/(machpoints[3,i]))*(((2+(k-1)*machpoints[3,i]^(2))/(2+(k-1)*machpoints[1,i]^(2))))^(((k+1)/(2*(k-1))))
: m[2,i]:=machpoints[3,i]
:EndIf
:
:If machpoints[2,i]=3 Then
:poypox[1,i]:=((1+k*machpoints[1,i]^(2))/(1+k*machpoints[3,i]^(2)))*(((1+(((k-1)*machpoints[3,i]^(2))/(2)))/(1+(((k-1)*machpoints[1,i]^(2))/(2)))))^(((k+1)/(2*(k-1))))
: m[2,i]:=machpoints[3,i]
:EndIf
:
: p0t[1,i]:=p0aux
: p0aux:=p0aux .* poypox[1,i]
: p0t[2,i]:=p0aux
: ds[1,i]:=−r .* ln(poypox[1,i])
:
:EndIf
:
:EndFor
:
:
:d0:=((p0t)/(r*t0))
:acr:=√(((2*k)/(k+1))*r*t0)
:areacritical:=((massflowrate)/(√(2)*√(((k)/(2))*(((2)/(k+1)))^(((k+1)/(k-1)))))) .* (p0t .* d0) .^ (−0.5)
:
:tcritical:=((t0*2)/(k+1))
:pcritical:=p0t .* (((2)/(k+1)))^(((k)/(k-1)))
:dcritical:=d0 .* (((2)/(k+1)))^(((1)/(k-1)))
:
:t:=(m .^ 2 .* (k-1) ./ 2 .+ 1) .^ (−1)
:p:=(m .^ 2 .* (k-1) ./ 2 .+ 1) .^ ((k)/(1-k))
:d:=(m .^ 2 .* (k-1) ./ 2 .+ 1) .^ ((1)/(1-k))
:lambda:=((k+1) .* m .^ 2 ./ ((k-1) .* m .^ 2 .+ 2)) .^ (0.5)
:arearatio:=(((2)/(k+1)) .* (m .^ 2 .* (k-1) ./ 2 .+ 1)) .^ ((k+1)/(2*(k-1))) .* m .^ (−1)
:area:=arearatio .* areacritical
:For i,1,colDim(machpoints),1
:area[2,i]:=ifFn(string(area[2,i])="undef",0,area[2,i])
:EndFor
:
:Disp "OTHER VALUES"
:Disp "Lambda Max",(((k+1)/(k-1)))^(0.5)
:
:Disp "RATIOS"
:Disp "M",m
:Disp "T",t
:Disp "P",p
:Disp "Density",d
:Disp "λ",lambda
:Disp "A/Acr",arearatio
:
:Disp "CRITICAL"
:Disp "Critical Speed of sound",acr
:Disp "A critical",areacritical
:Disp "t critical",tcritical
:Disp "pressure critical",pcritical
:Disp "density critical",dcritical
:
:Disp "VALUES"
:Disp "M",m
:Disp "P0",p0t
:Disp "T0",t0
:Disp "d0",d0
:Disp "T",t .* t0
:Disp "P",p .* p0t
:Disp "Density",d .* d0
:Disp "wind speed c",lambda .* acr
:Disp "A",area
:Disp "Diameter",((4*area)/(π)) .^ (0.5)
:Disp "P0y/P0x",poypox
:Disp "Δs",ds
:
:EndPrgm
``` 
