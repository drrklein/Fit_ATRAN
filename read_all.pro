;RK 7.10.2016
;Read allmost all files

alt=38+lindgen(8)
pwv=[1,2,3,4,5,6,7,8,10,12,14,15,16,18,20,22,25,27,30,32,35,37,40,45,50]
ZA=30+5*lindgen(9)

atran=hash()
for a=0,n_elements(alt)-1 do for z=0,n_elements(ZA)-1 do for p=0,n_elements(pwv)-1 do begin
   fn=string(alt[a],ZA[z],pwv[p],format='(%"atran_%dK_%ddeg_%dpwv_40-300mum.txt.dat")')
   readcol,fn,n,wave,trans
   good=where(wave ge 50 and wave le 205)
   atran[fn]={wave:wave[good],trans:trans[good]}
endfor

save,atran,filename='atran_fifils.sav'

end
